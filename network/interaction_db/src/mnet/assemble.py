"""Assemble the final edges.csv and nodes.csv.

Will union ppi + catalysis + lipid_is_a edges (undirected, deduplicated, no
self-loops), attach score/source/edge_type/evidence from curation/edge_scores.csv,
and build nodes.csv with the documented columns and per-node degree.

Final edges.csv columns:
  node1, node2, node1_type, node2_type, edge_type, score, source, evidence,
  n_evidence, confidence
Final nodes.csv columns:
  node_id, node_type, label, is_measured, in_string, refmet_id, chebi_id,
  chebi_rhea, uniprot, gene_symbol, super_class, main_class, lipid_class,
  lipid_level, chains, swisslipids_id, mapping_method, confidence, degree

Phase 0: scaffold only.
"""
from __future__ import annotations

EDGE_COLUMNS = [
    "node1", "node2", "node1_type", "node2_type", "edge_type",
    "score", "source", "evidence", "n_evidence", "confidence",
]
NODE_COLUMNS = [
    "node_id", "node_type", "label", "is_measured", "in_string", "refmet_id",
    "chebi_id", "chebi_rhea", "uniprot", "gene_symbol", "super_class",
    "main_class", "lipid_class", "lipid_level", "chains", "swisslipids_id",
    "mapping_method", "confidence", "degree",
]
EDGE_TYPES = ("ppi", "catalysis", "lipid_is_a")


def _canon(df):
    import numpy as np
    a = np.minimum(df["node1"].values.astype(object), df["node2"].values.astype(object))
    b = np.maximum(df["node1"].values.astype(object), df["node2"].values.astype(object))
    return [f"{x}\t{y}" for x, y in zip(a, b)]


def run(cfg, args) -> int:
    import os
    import numpy as np
    import pandas as pd
    from . import metabolites as MB
    from . import network as net

    interim = MB._interim(cfg)
    outdir = (cfg.get("mnet", {}) or {}).get("output_dir", "data/output")
    outdir = outdir if os.path.isabs(outdir) else os.path.join(cfg["_project_root"], outdir)
    os.makedirs(outdir, exist_ok=True)
    curation = os.path.join(cfg["_project_root"], "curation")

    # ---- edge scores ----
    es = pd.read_csv(os.path.join(curation, "edge_scores.csv"), dtype=str)
    escore = {r["edge_type"]: r for _, r in es.iterrows()}

    # ---- ppi ----
    ppi_path = net.ppi_path(cfg)
    ppi = net.load_ppi(cfg)
    assert ppi["combined_score"].min() >= 500, "ppi min score < 500!"
    print(f"[mnet.assemble] ppi {ppi_path} rows={len(ppi):,} min={ppi['combined_score'].min():.1f}")
    ppi_e = pd.DataFrame({
        "node1": ppi["protein1"].astype(str), "node2": ppi["protein2"].astype(str),
        "node1_type": "protein", "node2_type": "protein", "edge_type": "ppi",
        "score": ppi["combined_score"].astype("float32"),
        "source": "STRING v12 (physical else 0.9x full)",
        "evidence": None, "n_evidence": pd.NA, "confidence": "from_string",
        "shared_generic": False,
    })

    # ---- node attributes (loaded early to reconcile duplicate measured nodes) ----
    p2 = pd.read_parquet(os.path.join(interim, "metabolites_nonlipid.parquet"))
    lip = pd.read_parquet(os.path.join(interim, "lipids.parquet"))
    cross = pd.read_parquet(os.path.join(interim, "protein_crosswalk.parquet"))
    names = pd.read_parquet(os.path.join(interim, "chebi_names.parquet")).set_index("chebi_id")["name"].to_dict()
    genes = pd.read_parquet(os.path.join(interim, "uniprot_genes.parquet"))
    uni2gene = dict(zip(genes["uniprot"], genes["gene_symbol"]))
    currency = set(pd.read_csv(os.path.join(curation, "currency_metabolites.csv")).chebi_id.dropna())
    instr = dict(zip(cross["node_id"], cross["in_string"]))
    p2i = p2.set_index("node_id")
    lipi = lip.drop_duplicates("node_id").set_index("node_id")

    # reconcile duplicate measured nodes by refmet_id: prefer the Phase-3 lipid
    # node, REDIRECT the dropped node's edges to it, fold the class ChEBI into an
    # attribute.
    p2_ref = (p2.dropna(subset=["refmet_id"]).drop_duplicates("refmet_id")
              .set_index("refmet_id")["node_id"].to_dict())
    lip_ref = (lip.dropna(subset=["refmet_id"]).drop_duplicates("refmet_id")
               .set_index("refmet_id")["node_id"].to_dict())
    node_remap = {}
    fold_chebi = {}
    for rm, lipnode in lip_ref.items():
        p2node = p2_ref.get(rm)
        if p2node and p2node != lipnode:
            node_remap[p2node] = lipnode
            cc = p2i.loc[p2node].get("chebi_rhea") if p2node in p2i.index else None
            fold_chebi[lipnode] = cc
    drop_nodes = set(node_remap)

    # canonical node -> all refmet ids landing on it (two input rows for the same
    # molecule, e.g. EPA and "Eicosapentaenoic acid", collapse onto one node and
    # keep BOTH refmet ids so no input row is lost).
    node_refmets: dict = {}
    for _, r in p2.iterrows():
        if pd.notna(r["refmet_id"]):
            nid = node_remap.get(r["node_id"], r["node_id"])
            node_refmets.setdefault(nid, set()).add(r["refmet_id"])
    for _, r in lip.iterrows():
        if pd.notna(r["refmet_id"]):
            node_refmets.setdefault(r["node_id"], set()).add(r["refmet_id"])

    def remap(df):
        df = df.copy()
        df["node1"] = df["node1"].map(lambda n: node_remap.get(n, n))
        df["node2"] = df["node2"].map(lambda n: node_remap.get(n, n))
        return df

    # ---- rhea (catalysis/transport) ----
    er = pd.read_parquet(os.path.join(interim, "edges_rhea.parquet"))
    er["score"] = er["edge_type"].map(lambda t: float(escore[t]["score"]))
    er["source"] = er["edge_type"].map(lambda t: f"Rhea ({escore[t]['source']})")
    er["confidence"] = er["edge_type"].map(lambda t: escore[t]["confidence"])
    er = remap(er)[["node1", "node2", "node1_type", "node2_type", "edge_type", "score",
                    "source", "evidence", "n_evidence", "confidence", "shared_generic"]]

    # ---- lipid_is_a ----
    el = pd.read_parquet(os.path.join(interim, "edges_lipid.parquet"))
    el["score"] = float(escore["lipid_is_a"]["score"])
    el["source"] = f"lipid hierarchy ({escore['lipid_is_a']['source']})"
    el["confidence"] = escore["lipid_is_a"]["confidence"]
    el["evidence"] = None
    el["n_evidence"] = pd.NA
    el["shared_generic"] = False
    el = remap(el)[["node1", "node2", "node1_type", "node2_type", "edge_type", "score",
                    "source", "evidence", "n_evidence", "confidence", "shared_generic"]]

    edges = pd.concat([ppi_e, er, el], ignore_index=True)
    edges = edges[edges["node1"] != edges["node2"]]                 # no self-loops
    edges["_k"] = _canon(edges)
    # for duplicate (pair, edge_type) keep the max score
    edges = (edges.sort_values("score", ascending=False)
             .drop_duplicates(["_k", "edge_type"]).drop(columns="_k").reset_index(drop=True))
    assert edges["score"].notna().all(), "edges without a score exist!"

    # all node ids = edge endpoints + every measured metabolite (even 0-degree)
    all_nodes = set(edges["node1"]) | set(edges["node2"])
    all_nodes |= (set(p2["node_id"]) | set(lip["node_id"])) - drop_nodes

    deg = pd.concat([edges["node1"], edges["node2"]]).value_counts().to_dict()

    def ntype(nid):
        if nid.startswith("LIPIDCLASS:"):
            return "lipid_class"
        if nid.startswith(("CHEBI:", "REFMET:", "MEAS:")):
            return "metabolite"
        return "protein"

    rows = []
    for nid in sorted(all_nodes):
        t = ntype(nid)
        rec = {"node_id": nid, "node_type": t, "label": None, "is_measured": False,
               "in_string": None, "refmet_id": None, "chebi_id": None, "chebi_rhea": None,
               "uniprot": None, "gene_symbol": None, "super_class": None, "main_class": None,
               "lipid_class": None, "lipid_level": None, "chains": None, "swisslipids_id": None,
               "mapping_method": None, "confidence": None, "is_currency": False,
               "degree": int(deg.get(nid, 0))}
        if t == "protein":
            rec["in_string"] = bool(instr.get(nid, True))
            rec["uniprot"] = nid if not nid.startswith("9606.") and "-" not in nid else None
            g = uni2gene.get(nid)
            rec["gene_symbol"] = g
            rec["label"] = g if isinstance(g, str) else nid
        elif t == "lipid_class":
            rec["label"] = nid.split(":", 1)[1]
        else:  # metabolite
            if nid in p2i.index:
                r = p2i.loc[nid]
                r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
                rec.update(is_measured=True, refmet_id=r.get("refmet_id"),
                           chebi_id=r.get("chebi_id"), chebi_rhea=r.get("chebi_rhea"),
                           super_class=r.get("super_class"), main_class=r.get("main_class"),
                           mapping_method=r.get("mapping_method"), confidence=r.get("confidence"),
                           is_currency=bool(r.get("is_currency")), label=r.get("metabolite"))
            elif nid in lipi.index:
                r = lipi.loc[nid]
                cid = r.get("exact_chebi")
                if (cid is None or pd.isna(cid)) and nid in fold_chebi:
                    cid = fold_chebi[nid]              # folded class ChEBI from the dropped dup
                rec.update(is_measured=True, refmet_id=r.get("refmet_id"),
                           chebi_id=cid, chebi_rhea=cid, lipid_class=r.get("lipid_class"),
                           lipid_level=r.get("lipid_level"), chains=r.get("chains"),
                           swisslipids_id=r.get("swisslipids_id"),
                           mapping_method=r.get("mapping_method"), confidence=r.get("confidence"),
                           label=r.get("metabolite"))
            if nid.startswith("CHEBI:"):
                rec["chebi_id"] = rec["chebi_id"] or nid
                if not rec["label"]:
                    rec["label"] = names.get(nid, nid)
                rec["is_currency"] = rec["is_currency"] or (nid in currency)
            if nid in node_refmets:                    # keep every refmet on the node
                rec["refmet_id"] = ";".join(sorted(node_refmets[nid]))
                rec["is_measured"] = True
            if not rec["label"]:
                rec["label"] = nid
        rows.append(rec)
    nodes = pd.DataFrame(rows)

    # ---- write ----
    edges.to_csv(os.path.join(outdir, "edges.csv"), index=False)
    # compact parquet copy of the full edge list (easier to share than the CSV)
    edges.to_parquet(os.path.join(outdir, "edges.parquet"), index=False, compression="zstd")
    nodes.to_csv(os.path.join(outdir, "nodes.csv"), index=False)
    # edges_3col: ONE row per undirected pair (collapse catalysis+transport dup
    # pairs, keep the max score). edges.csv keeps one row per (pair, edge_type).
    e3 = edges[["node1", "node2", "score"]].copy()
    e3["_k"] = _canon(e3)
    e3 = (e3.sort_values("score", ascending=False)
          .drop_duplicates("_k").drop(columns="_k"))
    e3["node1"] = e3["node1"].astype("string")
    e3["node2"] = e3["node2"].astype("string")
    e3["score"] = e3["score"].astype("float32")
    e3 = e3.sort_values("score", ascending=False, kind="mergesort").reset_index(drop=True)
    e3.to_parquet(os.path.join(outdir, "edges_3col.parquet"), index=False)
    _write_output_readme(cfg, outdir, edges, nodes)

    print(f"[mnet.assemble] edges={len(edges):,} nodes={len(nodes):,} "
          f"by_type={edges['edge_type'].value_counts().to_dict()}")
    return 0


def _write_output_readme(cfg, outdir, edges, nodes):
    import os
    from string_network.download import fetch_uniprot_release
    L = ["# Metabolite-extended STRING network (output)", ""]
    L.append("Protein-protein (STRING) + protein-metabolite (Rhea catalysis/transport) "
             "+ lipid->class (lipid_is_a) undirected edge list.")
    L.append("")
    L.append("## The three output files")
    L.append("- `edges.csv` - **full detail**: one row per (undirected pair, edge_type) "
             "with source/evidence/confidence/shared_generic.")
    L.append("- `edges_3col.parquet` - **the same network in the STRING/watershed layout**: "
             "node1, node2 (string), score (float32); ONE row per undirected pair "
             "(catalysis+transport duplicates collapsed to the max score); sorted by "
             "score desc; no index.")
    L.append("- `edges.parquet` - same columns as edges.csv, zstd-compressed (easier to "
             "share than the large CSV).")
    L.append("- `nodes.csv` - **node attributes** (every measured metabolite listed, even "
             "degree 0).")
    L.append("")
    L.append("> `degree` in nodes.csv = number of rows in edges.csv touching the node "
             "(a pair with both catalysis and transport counts twice), so it can be "
             "slightly higher than the neighbour count in edges_3col.parquet.")
    L.append("")
    L.append("> **Scores:** only PPI scores are real STRING combined scores. The non-PPI "
             "scores are **placeholders** to be revised: catalysis = 900, transport = 800, "
             "lipid_is_a = 500 (see curation/edge_scores.csv).")
    L.append("")
    L.append("## edges.csv columns")
    L.append("node1, node2, node1_type, node2_type, edge_type, score, source, evidence, "
             "n_evidence, confidence, shared_generic")
    L.append("")
    L.append("## edge_type / source / score placeholders")
    L.append("- ppi: STRING v12 combined score (physical else 0.9x full), >= 500")
    L.append("- catalysis / transport: Rhea reactions (evidence = RHEA master IDs); "
             "placeholder scores 900 / 800 (curation/edge_scores.csv)")
    L.append("- lipid_is_a: lipid -> lipid-class hierarchy; placeholder score 500")
    L.append("- scores for the metabolite layers are curation placeholders, to be revised.")
    L.append("")
    L.append("## node types")
    L.append("protein | metabolite | lipid_class (never parse the id - use node_type).")
    L.append("")
    L.append(f"## versions / dates")
    L.append(f"- STRING v{cfg.get('string_version')} (taxon {cfg.get('taxon')})")
    L.append(f"- UniProt: {fetch_uniprot_release(cfg)}")
    L.append("- Rhea, ChEBI, SwissLipids: current releases (see data/raw/MANIFEST.txt)")
    L.append("")
    L.append("## counts")
    L.append(f"- edges: {len(edges):,} {edges['edge_type'].value_counts().to_dict()}")
    L.append(f"- nodes: {len(nodes):,} {nodes['node_type'].value_counts().to_dict()}")
    L.append("")
    L.append("## Licence / citation")
    L.append("STRING, Rhea, ChEBI, SwissLipids and UniProt are all released under CC BY 4.0.")
    L.append("Cite: STRING - Szklarczyk et al., NAR 2023; Rhea - Bansal et al., NAR 2022; "
             "ChEBI - Hastings et al., NAR 2016; SwissLipids - Aimo et al., Bioinformatics 2015; "
             "UniProt - The UniProt Consortium, NAR 2023.")
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
