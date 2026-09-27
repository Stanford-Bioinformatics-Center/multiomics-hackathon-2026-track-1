"""Protein helpers: gene symbols / HGNC ids from the existing UniProt idmapping.

Parses the already-downloaded HUMAN_9606_idmapping.dat.gz (types Gene_Name and
HGNC) into uniprot -> gene_symbol, hgnc_id. The by-organism idmapping file does
not contain UniProt secondary->primary accessions (no such type is present), so
that crosswalk is left for Phase 4; this is reported, not guessed.
"""
from __future__ import annotations

import os

import pandas as pd

from string_network.load import (load_mapping_with_fallback, strip_isoform,
                                  strip_taxon)

from . import network as net


def map_accessions_to_nodes(cfg, accessions) -> pd.DataFrame:
    """Reusable crosswalk: UniProt accession -> STRING (>=ppi_threshold) node id.

    Chain (same mapping options as the STRING build, read from config.yaml):
      1. isoform-stripped accession already a network node id     -> string_direct
      2. accession -> ENSP (idmapping STRING) -> node via
         load_mapping_with_fallback                                -> via_ensp
      3. same gene (Gene_Name/HGNC) with a unique network node     -> via_gene
      4. otherwise keep the accession as a new node (in_string=False) -> new_node

    Returns accession, node_id, method, in_string.
    """
    string_nodes = net.ppi_nodes(cfg)
    taxon = int(cfg.get("taxon", 9606))
    m = load_mapping_with_fallback(
        cfg, cfg.get("mapping_source", "uniprot_idmapping"), taxon,
        alias_sources=cfg.get("mapping_alias_sources", ["UniProt_AC", "Ensembl_UniProt"]))
    ensp2node = dict(zip(m["ensp"], m["uniprot"]))

    raw = cfg["paths"]["raw"]
    idm = pd.read_csv(os.path.join(raw, "HUMAN_9606_idmapping.dat.gz"),
                      sep="\t", header=None, names=["uniprot", "type", "id"], dtype=str)
    st = idm[idm["type"] == "STRING"].copy()
    st["ensp"] = strip_taxon(st["id"].astype("string")).astype(str)
    uni2ensp = st.groupby("uniprot")["ensp"].apply(list).to_dict()

    genes = parse_uniprot_genes(os.path.join(raw, "HUMAN_9606_idmapping.dat.gz"))
    uni2gene = dict(zip(genes["uniprot"], genes["gene_symbol"]))
    gene2nodes: dict = {}
    for uni, g in uni2gene.items():
        if isinstance(g, str) and uni in string_nodes:
            gene2nodes.setdefault(g, set()).add(uni)

    rows = []
    for acc in accessions:
        a = strip_isoform(pd.Series([acc], dtype="string")).iloc[0]
        if a in string_nodes:
            rows.append((acc, a, "string_direct", True)); continue
        node = None
        for e in uni2ensp.get(a, []):
            n = ensp2node.get(e)
            if n:
                node = n; break
        if node is not None:
            rows.append((acc, node, "via_ensp", node in string_nodes)); continue
        g = uni2gene.get(a)
        if isinstance(g, str) and len(gene2nodes.get(g, set())) == 1:
            rows.append((acc, next(iter(gene2nodes[g])), "via_gene", True)); continue
        rows.append((acc, a, "new_node", False))
    return pd.DataFrame(rows, columns=["accession", "node_id", "method", "in_string"])


def parse_uniprot_genes(idmapping_path: str) -> pd.DataFrame:
    """uniprot | gene_symbol | hgnc_id from the idmapping .dat.gz."""
    df = pd.read_csv(idmapping_path, sep="\t", header=None,
                     names=["uniprot", "type", "id"], dtype=str)
    df["uniprot"] = strip_isoform(df["uniprot"].astype("string")).astype(str)

    gene = (df[df["type"] == "Gene_Name"][["uniprot", "id"]]
            .rename(columns={"id": "gene_symbol"})
            .drop_duplicates("uniprot"))
    hgnc = (df[df["type"] == "HGNC"][["uniprot", "id"]]
            .rename(columns={"id": "hgnc_id"})
            .drop_duplicates("uniprot"))
    out = gene.merge(hgnc, on="uniprot", how="outer")
    return out.reset_index(drop=True)


def secondary_to_primary(idmapping_path: str) -> pd.DataFrame:
    """Return an empty table + reason: no secondary accession type is present."""
    return pd.DataFrame(columns=["secondary", "primary"])


def run(cfg, args) -> int:
    print("[mnet.proteins] use the mnet-download step (download.run).")
    return 0
