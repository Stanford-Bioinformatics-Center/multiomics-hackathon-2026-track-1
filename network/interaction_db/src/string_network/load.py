"""Readers for the STRING links, UniProt idmapping and STRING aliases files.

Everything is memory-conscious: the links file is read in chunks with narrow
dtypes, and mapping tables are reduced to a two-column ``ensp -> uniprot`` frame.
"""
from __future__ import annotations

import os
import re
from typing import Iterator, Optional

import numpy as np
import pandas as pd

ISOFORM_RE = re.compile(r"-\d+$")

# Official UniProt accession pattern (6 or 10 chars). Used to keep gene symbols
# and other non-accession alias tokens out of the mapping.
UNIPROT_AC_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

# STRING ids look like "9606.ENSP00000000233"; we keep the bare "ENSP..." form
# as the join key so it matches idmapping / aliases which vary in prefixing.


def strip_taxon(series: pd.Series) -> pd.Series:
    """Turn '9606.ENSP000...' into 'ENSP000...'. Leaves other ids untouched."""
    return series.str.replace(r"^\d+\.", "", regex=True)


def strip_isoform(series: pd.Series) -> pd.Series:
    """Remove a trailing '-<digits>' isoform suffix from UniProt accessions."""
    return series.str.replace(ISOFORM_RE, "", regex=True)


def iter_links(path: str, threshold: Optional[int], chunksize: int,
               apply_threshold: bool = True) -> Iterator[pd.DataFrame]:
    """Yield chunks of the STRING links file.

    Columns returned: protein1, protein2 (bare ENSP, category), combined_score
    (int16). When *apply_threshold* is True, rows below *threshold* are dropped
    inside each chunk to keep memory low.
    """
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        header=0,
        names=["protein1", "protein2", "combined_score"],
        dtype={"protein1": "string", "protein2": "string", "combined_score": "int32"},
        chunksize=chunksize,
        engine="c",
    )
    for chunk in reader:
        if apply_threshold and threshold is not None:
            chunk = chunk[chunk["combined_score"] >= threshold]
        if chunk.empty:
            continue
        chunk["protein1"] = strip_taxon(chunk["protein1"]).astype("string")
        chunk["protein2"] = strip_taxon(chunk["protein2"]).astype("string")
        chunk["combined_score"] = chunk["combined_score"].astype("int16")
        yield chunk


def load_uniprot_idmapping(path: str) -> pd.DataFrame:
    """Load HUMAN_9606_idmapping.dat.gz, keeping only STRING rows.

    Returns a frame ``ensp | uniprot`` where uniprot is the isoform-stripped
    accession and ensp is the bare ENSP id.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["uniprot", "type", "id"],
        dtype="string",
    )
    df = df[df["type"] == "STRING"]
    df = df.rename(columns={"id": "ensp"})[["ensp", "uniprot"]]
    df["ensp"] = strip_taxon(df["ensp"]).astype("string")
    df["uniprot"] = strip_isoform(df["uniprot"]).astype("string")
    df = df.dropna().drop_duplicates()
    return df.reset_index(drop=True)


def load_string_aliases(path: str, source: str) -> pd.DataFrame:
    """Load STRING aliases and keep one alias *source*.

    ``source`` is a value from the aliases 'source' column, e.g. 'UniProt_AC'
    or 'Ensembl_UniProt'. Returns ``ensp | uniprot``.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=["ensp", "alias", "source"],
        dtype="string",
    )
    # The source column may list multiple space-joined sources; match membership.
    mask = df["source"].str.split(r"\s+").apply(lambda xs: source in xs if xs is not None else False)
    df = df[mask]
    df = df.rename(columns={"alias": "uniprot"})[["ensp", "uniprot"]]
    df["ensp"] = strip_taxon(df["ensp"]).astype("string")
    df["uniprot"] = strip_isoform(df["uniprot"]).astype("string")
    df = df.dropna().drop_duplicates()
    return df.reset_index(drop=True)


def load_mapping(cfg: dict, mapping_source: str) -> pd.DataFrame:
    """Dispatch to the configured mapping source. Returns ``ensp | uniprot``."""
    raw = cfg["paths"]["raw"]
    if mapping_source == "uniprot_idmapping":
        path = os.path.join(raw, os.path.basename(cfg["urls"]["uniprot_idmapping"]))
        return load_uniprot_idmapping(path)
    if mapping_source.startswith("string_aliases:"):
        source = mapping_source.split(":", 1)[1]
        path = os.path.join(raw, os.path.basename(cfg["urls"]["string_aliases"]))
        return load_string_aliases(path, source)
    raise ValueError(f"Unknown mapping_source: {mapping_source!r}")


def load_protein_info(path: str) -> pd.DataFrame:
    """Load 9606.protein.info: returns ``ensp | preferred_name`` (bare ENSP)."""
    df = pd.read_csv(
        path, sep="\t", header=0,
        names=["ensp", "preferred_name", "protein_size", "annotation"],
        usecols=["ensp", "preferred_name"], dtype="string",
    )
    df["ensp"] = strip_taxon(df["ensp"]).astype("string")
    return df.dropna().drop_duplicates("ensp").reset_index(drop=True)


def _pick_one_per_ensp(cand: pd.DataFrame, reviewed: set) -> pd.DataFrame:
    """From an ``ensp | uniprot`` candidate frame pick a single UniProt per ENSP.

    Preference: a reviewed (Swiss-Prot) accession, then the lexicographically
    smallest accession (stable, deterministic).
    """
    cand = cand.dropna().drop_duplicates()
    cand = cand.assign(_rev=cand["uniprot"].isin(reviewed).astype(int))
    cand = cand.sort_values(["ensp", "_rev", "uniprot"], ascending=[True, False, True])
    return cand.groupby("ensp", as_index=False).first()[["ensp", "uniprot"]]


DEFAULT_ALIAS_SOURCES = ("UniProt_AC", "Ensembl_UniProt")


def load_mapping_with_fallback(cfg: dict, mapping_source: str, taxon: int,
                               alias_sources=None,
                               return_source: bool = False) -> pd.DataFrame:
    """ENSP -> target id (one per ENSP) using the fallback chain:

    1. UniProt accession, resolved from *mapping_source* first, then from the
       STRING alias sources in *alias_sources* (default UniProt_AC then
       Ensembl_UniProt) for any ENSP still unmapped. When several accessions
       exist for an ENSP, a reviewed (Swiss-Prot) one is preferred.
    2. STRING preferred_name (from protein.info)
    3. the raw STRING id, e.g. '9606.ENSP...'

    ``alias_sources`` lets callers (e.g. the tuner) vary only this step without
    touching the resolution logic. With ``return_source=True`` an extra
    ``source`` column records which step produced each mapping.

    The universe of ENSPs is every protein in protein.info, so no STRING protein
    is dropped for lack of a UniProt accession. Returns ``ensp | uniprot`` (and
    ``source`` if requested).
    """
    if alias_sources is None:
        alias_sources = DEFAULT_ALIAS_SOURCES
    try:
        reviewed = load_reviewed_set(cfg)
    except FileNotFoundError:
        reviewed = set()

    primary = load_mapping(cfg, mapping_source)          # ensp -> uniprot
    uni = _pick_one_per_ensp(primary, reviewed)
    uni["source"] = mapping_source

    # Extend UniProt coverage using the STRING aliases, keeping only real UniProt
    # accessions and preferring a reviewed one. This resolves proteins missing
    # from the primary idmapping (e.g. calmodulin ENSP -> P0DP24) instead of
    # falling back to a gene symbol.
    aliases_path = os.path.join(cfg["paths"]["raw"],
                                os.path.basename(cfg["urls"]["string_aliases"]))
    for src in alias_sources:
        try:
            aliases = load_string_aliases(aliases_path, src)
        except (FileNotFoundError, KeyError):
            continue
        aliases = aliases[aliases["uniprot"].str.match(UNIPROT_AC_RE).fillna(False)]
        extra = aliases[~aliases["ensp"].isin(set(uni["ensp"]))]
        extra = _pick_one_per_ensp(extra, reviewed)
        extra["source"] = f"aliases:{src}"
        uni = pd.concat([uni, extra], ignore_index=True)

    info_path = os.path.join(cfg["paths"]["raw"],
                             os.path.basename(cfg["urls"]["protein_info"]))
    info = load_protein_info(info_path)

    mapped_ensp = set(uni["ensp"])
    fallback = info[~info["ensp"].isin(mapped_ensp)].copy()
    pref = fallback["preferred_name"]
    raw_id = (str(taxon) + "." + fallback["ensp"]).astype("string")
    has_pref = pref.notna() & (pref.str.len() > 0)
    fallback["uniprot"] = pref.where(has_pref, raw_id)
    fallback["source"] = np.where(has_pref, "preferred_name", "raw_string_id")
    fallback = fallback[["ensp", "uniprot", "source"]]

    combined = pd.concat([uni, fallback], ignore_index=True)
    combined = combined.dropna(subset=["ensp", "uniprot"]).drop_duplicates("ensp")
    combined = combined.reset_index(drop=True)
    if return_source:
        return combined[["ensp", "uniprot", "source"]]
    return combined[["ensp", "uniprot"]]


def load_reviewed_set(cfg: dict) -> set:
    """Load the Swiss-Prot reviewed accession list as a set."""
    path = os.path.join(cfg["paths"]["raw"], f"reviewed_{cfg['taxon']}.list")
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_keep_list(path: str) -> set:
    """Load a newline-delimited list of UniProt IDs to keep."""
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}
