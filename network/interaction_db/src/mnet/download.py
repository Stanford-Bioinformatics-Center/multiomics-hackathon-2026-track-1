"""Download + cache Rhea, ChEBI and SwissLipids, and build interim parquet caches.

`make mnet-download` runs :func:`run`, which:
  1. downloads the new sources into data/raw (reusing existing STRING/UniProt
     files), refreshing the shared MANIFEST.txt;
  2. parses them into tidy parquet caches in data/interim/;
  3. fills curation/currency_metabolites.csv from ChEBI + Rhea participants;
  4. writes reports/mnet_phase1_downloads.md.

Reaction participants come from the Rhea SPARQL endpoint (see mnet.rhea).
"""
from __future__ import annotations

import gzip
import os
from typing import Dict, List

import pandas as pd

from string_network.download import download_one, write_manifest
from string_network.download import fetch_uniprot_release

from . import chebi as chebi_mod
from . import proteins as prot_mod
from . import rhea as rhea_mod

# Special (non-derivable) filenames for the manifest / download.
SPECIAL_FILENAMES = {
    "reviewed_list": None,               # handled by string_network naming
    "swisslipids_lipids": "swisslipids_lipids.tsv.gz",
}

CURRENCY = [
    ("water", ["water"], "solvent"),
    ("hydron", ["hydron", "proton"], "proton"),
    ("ATP", ["ATP"], "energy currency"),
    ("ADP", ["ADP"], "energy currency"),
    ("AMP", ["AMP"], "energy currency"),
    ("GTP", ["GTP"], "energy currency"),
    ("GDP", ["GDP"], "energy currency"),
    ("NAD(+)", ["NAD(+)", "NAD"], "redox cofactor"),
    ("NADH", ["NADH"], "redox cofactor"),
    ("NADP(+)", ["NADP(+)", "NADP"], "redox cofactor"),
    ("NADPH", ["NADPH"], "redox cofactor"),
    ("FAD", ["FAD"], "redox cofactor"),
    ("FADH2", ["FADH2"], "redox cofactor"),
    ("CO2", ["carbon dioxide"], "small inorganic"),
    ("O2", ["dioxygen"], "small inorganic"),
    ("H2O2", ["hydrogen peroxide"], "small inorganic"),
    ("phosphate", ["phosphate", "hydrogenphosphate", "dihydrogenphosphate"], "leaving group"),
    ("diphosphate", ["diphosphate"], "leaving group"),
    ("CoA", ["coenzyme A"], "acyl carrier"),
    ("ammonium", ["ammonium"], "small inorganic"),
    ("hydrogencarbonate", ["hydrogencarbonate", "hydrogen carbonate"], "small inorganic"),
    ("S-adenosyl-L-methionine", ["S-adenosyl-L-methionine"], "methyl donor"),
    ("S-adenosyl-L-homocysteine", ["S-adenosyl-L-homocysteine"], "methyl donor"),
    ("UDP", ["UDP"], "nucleotide leaving group"),
    ("CTP", ["CTP"], "nucleotide"),
    ("CDP", ["CDP"], "nucleotide"),
    ("ubiquinone", ["ubiquinone", "ubiquinone-10"], "electron carrier"),
    ("ubiquinol", ["ubiquinol", "ubiquinol-10"], "electron carrier"),
    ("reduced flavodoxin", ["reduced flavodoxin"], "electron carrier"),
    ("oxidized flavodoxin", ["oxidized flavodoxin"], "electron carrier"),
]


def _fname(cfg, key):
    if key in SPECIAL_FILENAMES and SPECIAL_FILENAMES[key]:
        return SPECIAL_FILENAMES[key]
    url = cfg["urls"][key]
    return url.split("?")[0].rstrip("/").split("/")[-1]


def _download_sources(cfg, mcfg) -> List[str]:
    raw = cfg["paths"]["raw"]
    keys = ["rhea2uniprot_sprot", "rhea_directions", "rhea2xrefs",
            "rhea_chebi_ph7_3", "chebi_obo", "chebi_structures",
            "swisslipids_lipids"]
    if mcfg.get("use_trembl"):
        keys.append("rhea2uniprot_trembl")
    fetched = []
    for k in keys:
        url = cfg["urls"].get(k)
        if not url:
            print(f"[mnet.download] no url for {k}, skipping")
            continue
        try:
            fetched.append(download_one(url, raw, filename=SPECIAL_FILENAMES.get(k)))
        except Exception as exc:  # noqa: BLE001
            print(f"[mnet.download] ERROR downloading {k}: {exc}")
    return fetched


def _refresh_manifest(cfg):
    raw = cfg["paths"]["raw"]
    files = {}
    for key, url in cfg["urls"].items():
        if key == "rhea_sparql":
            continue
        fname = (f"reviewed_{cfg['taxon']}.list" if key == "reviewed_list"
                 else _fname(cfg, key))
        path = os.path.join(raw, fname)
        if os.path.exists(path):
            files[url] = path
    versions = {
        "string_version": cfg.get("string_version", "?"),
        "taxon": cfg.get("taxon", "?"),
        "uniprot_release": fetch_uniprot_release(cfg),
        "build_threshold (config default)": cfg.get("threshold", "?"),
    }
    write_manifest(raw, files, versions)


def _norm_name(s: str) -> str:
    """Lowercased name with trailing charge / 'zwitterion' / 'residue' removed."""
    import re
    if not isinstance(s, str):
        return ""
    s = s.strip()
    s = re.sub(r"\(\d*[+-]\)$", "", s).strip()      # drop charge suffix e.g. (4-)
    s = re.sub(r"\s+(zwitterion|residue)$", "", s, flags=re.I).strip()
    return s.lower()


def _build_currency(names_df, participants_chebi, out_path) -> pd.DataFrame:
    """Match currency concepts to ChEBI ids (every charge form) present in Rhea.

    Matches on the normalised name (charge / 'zwitterion' stripped) OR on an
    exact synonym, then keeps only ids that actually occur in rhea_participants.
    """
    names_df = names_df.copy()
    names_df["norm"] = names_df["name"].map(_norm_name)
    # explode synonyms for exact matching
    syn = names_df[["chebi_id", "synonyms"]].dropna(subset=["synonyms"]).copy()
    syn["syn"] = syn["synonyms"].str.split("|")
    syn = syn.explode("syn")
    syn["syn_l"] = syn["syn"].str.strip().str.lower()

    rows = []
    for display, bases, reason in CURRENCY:
        bases_l = {b.lower() for b in bases}
        by_name = set(names_df[names_df["norm"].isin(bases_l)]["chebi_id"])
        by_syn = set(syn[syn["syn_l"].isin(bases_l)]["chebi_id"])
        cands = (by_name | by_syn) & participants_chebi
        if not cands:
            rows.append({"chebi_id": "", "name": display,
                         "reason": f"{reason} (NOT FOUND in Rhea participants)"})
            continue
        nm_map = dict(zip(names_df["chebi_id"], names_df["name"]))
        for cid in sorted(cands):
            rows.append({"chebi_id": cid, "name": nm_map.get(cid, display), "reason": reason})
    df = pd.DataFrame(rows, columns=["chebi_id", "name", "reason"])
    # preserve manually-curated rows (e.g. inorganic ions added in Phase 4) so a
    # re-run of mnet-download does not wipe them.
    if os.path.exists(out_path):
        prev = pd.read_csv(out_path, dtype=str)
        keep = prev[prev["reason"].fillna("").str.contains("inorganic ion|manual", regex=True)]
        keep = keep[~keep["chebi_id"].isin(set(df["chebi_id"]))]
        if len(keep):
            df = pd.concat([df, keep[["chebi_id", "name", "reason"]]], ignore_index=True)
    df.to_csv(out_path, index=False)
    return df


def run(cfg, args) -> int:
    mcfg = cfg.get("mnet", {}) or {}
    raw = cfg["paths"]["raw"]
    interim = mcfg.get("interim", "data/interim")
    interim = interim if os.path.isabs(interim) else os.path.join(cfg["_project_root"], interim)
    os.makedirs(interim, exist_ok=True)
    reports = cfg["paths"]["reports"]

    print("[mnet.download] downloading sources ...")
    _download_sources(cfg, mcfg)
    _refresh_manifest(cfg)

    def rp(key):
        return os.path.join(raw, _fname(cfg, key))

    # --- Rhea flat files ---
    rhea_enzymes = rhea_mod.parse_enzymes(rp("rhea2uniprot_sprot"))
    rhea_directions = rhea_mod.parse_directions(rp("rhea_directions"))
    rhea_xrefs = rhea_mod.parse_xrefs(rp("rhea2xrefs"))
    print(f"[mnet.download] rhea_enzymes={len(rhea_enzymes):,} "
          f"directions={len(rhea_directions):,} xrefs={len(rhea_xrefs):,}")

    # --- ChEBI ---
    print("[mnet.download] parsing chebi.obo ...")
    chebi_names, chebi_relations, chebi_structures = chebi_mod.parse_obo(rp("chebi_obo"))
    print("[mnet.download] loading InChIKeys from structures.tsv.gz ...")
    ik = chebi_mod.load_inchikeys(rp("chebi_structures"))
    chebi_structures = chebi_structures.merge(ik, on="chebi_id", how="left")
    print(f"[mnet.download] chebi_names={len(chebi_names):,} "
          f"relations={len(chebi_relations):,} structures={len(chebi_structures):,} "
          f"(with inchikey={chebi_structures['inchikey'].notna().sum():,})")

    formula_map: Dict[str, str] = dict(
        zip(chebi_structures["chebi_id"], chebi_structures["formula"]))
    alt2primary = chebi_mod.alt_to_primary(chebi_names)

    # --- Rhea participants (SPARQL) ---
    print("[mnet.download] fetching Rhea participants via SPARQL ...")
    participants = rhea_mod.build_participants(cfg["urls"]["rhea_sparql"], formula_map)
    participants["chebi_id"] = participants["chebi_id"].map(lambda c: alt2primary.get(c, c))
    print(f"[mnet.download] rhea_participants={len(participants):,} "
          f"reactions={participants['master_rhea_id'].nunique():,}")

    # --- UniProt genes ---
    uniprot_genes = prot_mod.parse_uniprot_genes(rp("uniprot_idmapping"))
    print(f"[mnet.download] uniprot_genes={len(uniprot_genes):,}")

    # --- write interim parquet ---
    def wp(df, name):
        p = os.path.join(interim, f"{name}.parquet")
        df.to_parquet(p, index=False)
        return p
    wp(rhea_enzymes, "rhea_enzymes")
    wp(rhea_directions, "rhea_directions")
    wp(rhea_xrefs, "rhea_xrefs")
    wp(participants, "rhea_participants")
    wp(chebi_names, "chebi_names")
    wp(chebi_relations, "chebi_relations")
    wp(chebi_structures, "chebi_structures")
    wp(uniprot_genes, "uniprot_genes")

    # swisslipids: cache a slim table if the file is present + valid gzip
    sl_status = _cache_swisslipids(cfg, raw, interim)

    # --- currency ---
    curation_dir = os.path.join(cfg["_project_root"], "curation")
    cur_df = _build_currency(chebi_names, set(participants["chebi_id"]),
                             os.path.join(curation_dir, "currency_metabolites.csv"))
    print(f"[mnet.download] currency rows written: {len(cur_df):,}")

    # --- report ---
    _write_report(cfg, mcfg, reports, raw, interim,
                  rhea_enzymes, rhea_directions, rhea_xrefs, participants,
                  chebi_names, chebi_structures, uniprot_genes, cur_df, sl_status)
    return 0


def _cache_swisslipids(cfg, raw, interim) -> str:
    path = os.path.join(raw, SPECIAL_FILENAMES["swisslipids_lipids"])
    if not os.path.exists(path):
        return "missing"
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
        opener = gzip.open if magic == b"\x1f\x8b" else open
        df = pd.read_csv(opener(path, "rt", encoding="utf-8", errors="replace"),
                         sep="\t", dtype=str, low_memory=False)
        df.to_parquet(os.path.join(interim, "swisslipids.parquet"), index=False)
        return f"ok ({len(df):,} rows, cols={list(df.columns)[:4]}...)"
    except Exception as exc:  # noqa: BLE001
        return f"invalid: {exc}"


def _write_report(cfg, mcfg, reports, raw, interim, rhea_enzymes, rhea_directions,
                  rhea_xrefs, participants, chebi_names, chebi_structures,
                  uniprot_genes, cur_df, sl_status):
    import numpy as np
    from string_network.config import resolve_output_path

    # human enzyme set: any accession in idmapping OR reviewed list
    idmap_ac = set(pd.read_csv(os.path.join(raw, "HUMAN_9606_idmapping.dat.gz"),
                               sep="\t", header=None, usecols=[0], dtype=str)[0])
    with open(os.path.join(raw, f"reviewed_{cfg['taxon']}.list")) as fh:
        reviewed = {ln.strip() for ln in fh if ln.strip()}
    human = idmap_ac | reviewed

    enz = rhea_enzymes.copy()
    enz["is_human"] = enz["uniprot"].isin(human)
    reac_any = enz["MASTER_ID"].dropna().nunique()
    reac_human = enz[enz.is_human]["MASTER_ID"].dropna().nunique()
    human_enzymes = enz[enz.is_human]["uniprot"].nunique()

    # STRING nodes (Phase 4 peek)
    try:
        snodes = pd.read_parquet(resolve_output_path(cfg))
        string_nodes = set(snodes.protein1) | set(snodes.protein2)
    except Exception:
        string_nodes = set()
    human_enz_in_string = len({u for u in enz[enz.is_human]["uniprot"].unique()
                               if u in string_nodes})

    # participants restricted to human reactions
    human_reacs = set(enz[enz.is_human]["MASTER_ID"].dropna().astype(int))
    ph = participants[participants["master_rhea_id"].isin(human_reacs)]
    chebi_human = ph["chebi_id"].nunique()
    by_type = ph.drop_duplicates(["chebi_id"]).groupby("compound_type").size().to_dict()
    top30 = (ph.groupby("chebi_id")["master_rhea_id"].nunique()
             .sort_values(ascending=False).head(30))
    name_map = dict(zip(chebi_names["chebi_id"], chebi_names["name"]))

    # metabolite csv check
    mpath = mcfg.get("metabolites", "data/input/metabolite_chebi_ids.csv")
    mpath = mpath if os.path.isabs(mpath) else os.path.join(cfg["_project_root"], mpath)
    mdf = pd.read_csv(mpath, dtype=str)

    def sz(key_or_file):
        p = os.path.join(raw, key_or_file)
        return f"{os.path.getsize(p):,}" if os.path.exists(p) else "-"

    L = ["# mnet Phase 1 - downloads & interim caches", ""]
    L.append(f"UniProt release: {fetch_uniprot_release(cfg)} | STRING v{cfg.get('string_version')}")
    L.append("")
    L.append("## Source files (data/raw; sizes; see MANIFEST.txt for SHA-256 + date)")
    L.append("")
    L.append("| file | size (bytes) |")
    L.append("|---|---|")
    src_files = ["rhea2uniprot_sprot.tsv", "rhea-directions.tsv", "rhea2xrefs.tsv",
                 "chebi_pH7_3_mapping.tsv", "chebi.obo.gz", "structures.tsv.gz",
                 "swisslipids_lipids.tsv.gz", "HUMAN_9606_idmapping.dat.gz"]
    for f in src_files:
        L.append(f"| {f} | {sz(f)} |")
    L.append("")
    L.append("Rhea participants: fetched from the SPARQL endpoint "
             f"({cfg['urls'].get('rhea_sparql')}); chose SPARQL over parsing "
             "rhea.rdf.gz with rdflib for speed/memory. ChEBI structures/InChIKey "
             "from the 2025 `flat_files/structures.tsv.gz` (obo lacks InChIKey; "
             "formula/charge/smiles come from the obo `chemrof:` property_values).")
    L.append("")
    L.append("## Interim parquet row counts")
    L.append("")
    L.append("| table | rows |")
    L.append("|---|---|")
    for name, df in [("rhea_enzymes", rhea_enzymes), ("rhea_directions", rhea_directions),
                     ("rhea_xrefs", rhea_xrefs), ("rhea_participants", participants),
                     ("chebi_names", chebi_names), ("chebi_structures", chebi_structures),
                     ("uniprot_genes", uniprot_genes)]:
        L.append(f"| {name} | {len(df):,} |")
    L.append(f"| chebi_structures with inchikey | {chebi_structures['inchikey'].notna().sum():,} |")
    L.append("")
    L.append("## Rhea")
    L.append("")
    L.append(f"- master reactions with participants: **{participants['master_rhea_id'].nunique():,}**")
    L.append(f"- reactions with >=1 enzyme (Swiss-Prot): **{reac_any:,}**")
    L.append(f"- reactions with >=1 HUMAN enzyme: **{reac_human:,}**")
    L.append(f"- human enzymes (UniProt accessions): **{human_enzymes:,}**")
    L.append(f"- of those already STRING node ids as-is: **{human_enz_in_string:,}** "
             f"(Phase 4 crosswalk peek; STRING nodes={len(string_nodes):,})")
    L.append(f"- transport reactions (same compound both sides): "
             f"{participants[participants.is_transport]['master_rhea_id'].nunique():,}")
    L.append("")
    if mcfg.get("use_trembl"):
        L.append("- TrEMBL: enabled.")
    else:
        L.append("- TrEMBL: disabled (`mnet.use_trembl=false`). Swiss-Prot covers most "
                 "human enzymes; enabling it would add TrEMBL-only human enzymes "
                 "(rhea2uniprot_trembl.tsv.gz, not downloaded).")
    L.append("")
    L.append("## ChEBI participants in human reactions")
    L.append("")
    L.append(f"- unique ChEBI ids in human reactions: **{chebi_human:,}**")
    L.append(f"- by compound_type: {by_type}")
    L.append("")
    L.append("### Top 30 ChEBI ids by number of human reactions")
    L.append("")
    L.append("| chebi_id | name | #human_reactions |")
    L.append("|---|---|---|")
    for cid, n in top30.items():
        L.append(f"| {cid} | {name_map.get(cid,'?')} | {n} |")
    L.append("")
    L.append("## Currency metabolites")
    L.append("")
    nf = cur_df[cur_df.chebi_id == ""]
    L.append(f"- rows written to curation/currency_metabolites.csv: {len(cur_df):,} "
             f"(charge forms present in Rhea participants)")
    if len(nf):
        L.append(f"- concepts NOT found in Rhea participants: "
                 f"{sorted(set(nf['name']))}")
    L.append("")
    L.append("## Metabolite CSV check")
    L.append("")
    L.append(f"- rows: {len(mdf):,} | columns ok: "
             f"{list(mdf.columns)==['metabolite','refmet_name','refmet_id','super_class','main_class','pubchem_cid','inchi_key','chebi_id','chebi_all','chebi_method']}")
    L.append(f"- with ChEBI: {mdf['chebi_id'].notna().sum():,}")
    L.append("")
    L.append("## SwissLipids")
    L.append("")
    L.append(f"- status: {sl_status}")
    L.append("")

    os.makedirs(reports, exist_ok=True)
    outp = os.path.join(reports, "mnet_phase1_downloads.md")
    with open(outp, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"[mnet.download] wrote {outp}")


# config url keys this step fetches (informational)
SOURCE_KEYS = ("rhea2uniprot_sprot", "rhea_directions", "rhea2xrefs",
               "rhea_chebi_ph7_3", "chebi_obo", "chebi_structures",
               "swisslipids_lipids")
