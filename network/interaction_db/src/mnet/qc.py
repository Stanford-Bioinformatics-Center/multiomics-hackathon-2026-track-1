"""Quality-control checks and report for the assembled network.

Will verify the ppi sub-network is identical to the STRING output, that edges are
undirected/deduplicated/self-loop-free, that every node id carries a node_type,
and will summarise node/edge counts by type into reports/mnet_*.md.

Phase 0: scaffold only.
"""
from __future__ import annotations


def run(cfg, args) -> int:
    import os
    import sys
    import numpy as np
    import pandas as pd
    import networkx as nx
    from . import metabolites as MB
    from . import network as net

    outdir = (cfg.get("mnet", {}) or {}).get("output_dir", "data/output")
    outdir = outdir if os.path.isabs(outdir) else os.path.join(cfg["_project_root"], outdir)
    reports = cfg["paths"]["reports"]
    curation = os.path.join(cfg["_project_root"], "curation")

    edges = pd.read_csv(os.path.join(outdir, "edges.csv"), low_memory=False)
    nodes = pd.read_csv(os.path.join(outdir, "nodes.csv"), low_memory=False)
    node_ids = set(nodes["node_id"])

    fails = []

    def check(cond, msg):
        if not cond:
            fails.append(msg)

    # 1. ppi identical to STRING >=500
    ppi = net.load_ppi(cfg)
    def canon_set(df, a, b):
        x = np.minimum(df[a].astype(str).values.astype(object), df[b].astype(str).values.astype(object))
        y = np.maximum(df[a].astype(str).values.astype(object), df[b].astype(str).values.astype(object))
        return {f"{p}\t{q}": s for p, q, s in zip(x, y, df["combined_score"] if "combined_score" in df else df["score"])}
    ppi_ref = canon_set(ppi, "protein1", "protein2")
    ppi_e = edges[edges["edge_type"] == "ppi"]
    ppi_got = canon_set(ppi_e, "node1", "node2")
    check(len(ppi_ref) == len(ppi_got) == len(ppi), f"ppi count mismatch {len(ppi_got)} vs {len(ppi)}")
    check(set(ppi_ref) == set(ppi_got), "ppi undirected pairs differ from STRING >=500")
    check(all(abs(ppi_ref[k] - ppi_got[k]) < 1e-3 for k in ppi_ref if k in ppi_got),
          "ppi scores differ from STRING >=500")
    check(ppi_e["score"].min() >= 500, "ppi min score < 500")

    # 2. structural
    check(not (edges["node1"] == edges["node2"]).any(), "self-loops present")
    dup = edges.assign(k=[f"{min(a,b)}\t{max(a,b)}" for a, b in zip(edges.node1, edges.node2)]) \
        .duplicated(["k", "edge_type"]).sum()
    check(dup == 0, f"{dup} duplicate (pair, edge_type) rows")
    missing = (set(edges["node1"]) | set(edges["node2"])) - node_ids
    check(not missing, f"{len(missing)} edge endpoints not in nodes.csv")
    # id prefixes
    mt = nodes[nodes.node_type == "metabolite"]["node_id"]
    check(mt.str.startswith(("CHEBI:", "REFMET:", "MEAS:")).all(), "metabolite id bad prefix")
    lcn = nodes[nodes.node_type == "lipid_class"]["node_id"]
    check(lcn.str.startswith("LIPIDCLASS:").all(), "lipid_class id bad prefix")

    # 3. no two protein nodes share an ENSP
    raw = cfg["paths"]["raw"]
    idm = pd.read_csv(os.path.join(raw, "HUMAN_9606_idmapping.dat.gz"), sep="\t",
                      header=None, names=["u", "t", "i"], dtype=str)
    st = idm[idm.t == "STRING"].copy()
    st["ensp"] = st["i"].str.replace(r"^\d+\.", "", regex=True)
    prot = set(nodes[nodes.node_type == "protein"]["node_id"])
    u2e = st[st.u.isin(prot)]
    ensp_to_nodes = u2e.groupby("ensp")["u"].nunique()
    check((ensp_to_nodes <= 1).all(),
          f"{int((ensp_to_nodes>1).sum())} ENSP map to >1 protein node")

    # 4. scores + currency
    check(edges["score"].notna().all(), "null scores")
    check(edges["score"].between(0, 1000).all(), "scores out of 0-1000")
    currency = set(pd.read_csv(os.path.join(curation, "currency_metabolites.csv")).chebi_id.dropna())
    measured_chebi_nodes = set(nodes[(nodes.node_type == "metabolite") & (nodes.is_measured)]["node_id"])
    pm = edges[edges.edge_type.isin(["catalysis", "transport"])]
    unmeasured_currency = [n for n in set(pm["node2"])
                           if n in currency and n not in measured_chebi_nodes]
    check(not unmeasured_currency,
          f"{len(unmeasured_currency)} unmeasured currency metabolites in catalysis/transport")

    # 5. no metabolite-metabolite edges except lipid_is_a
    tmap = dict(zip(nodes.node_id, nodes.node_type))
    def is_mm(r):
        return tmap.get(r.node1) == "metabolite" and tmap.get(r.node2) == "metabolite"
    mm = edges[edges.edge_type != "lipid_is_a"].apply(is_mm, axis=1)
    check(not mm.any(), f"{int(mm.sum())} metabolite-metabolite edges outside lipid_is_a")

    # 6. every input row represented on exactly one node; no refmet on two nodes.
    # (A handful of input rows are the same molecule under two names, e.g. EPA /
    # Eicosapentaenoic acid; these collapse onto one node carrying both refmet ids,
    # so the node count is slightly < 450 rather than duplicating a ChEBI node.)
    meas = nodes[(nodes.node_type == "metabolite") & (nodes.is_measured == True)]  # noqa: E712
    input_refmets = set(pd.read_csv(
        os.path.join(cfg["_project_root"], (cfg.get("mnet", {}) or {}).get("metabolites")),
        dtype=str)["refmet_id"].dropna())
    node_refmets = set()
    dup = []
    seen = set()
    for v in meas["refmet_id"].dropna():
        for r in str(v).split(";"):
            if r in seen:
                dup.append(r)
            seen.add(r); node_refmets.add(r)
    check(input_refmets <= node_refmets,
          f"{len(input_refmets - node_refmets)} input refmet ids missing a node: "
          f"{sorted(input_refmets - node_refmets)[:5]}")
    check(not dup, f"refmet_id on >1 node: {dup[:5]}")
    print(f"[mnet.qc] measured metabolite nodes: {len(meas)} "
          f"(input rows 450; same-molecule collapses: {450-len(meas)})")

    # 7. edges_3col: no duplicate undirected pairs
    e3 = pd.read_parquet(os.path.join(outdir, "edges_3col.parquet"))
    k3 = [f"{min(a,b)}\t{max(a,b)}" for a, b in zip(e3.node1.astype(str), e3.node2.astype(str))]
    check(len(k3) == len(set(k3)), f"{len(k3)-len(set(k3))} duplicate pairs in edges_3col")

    if fails:
        print("[mnet.qc] QC FAILED:")
        for f in fails:
            print("   -", f)
        sys.exit(2)
    print("[mnet.qc] all QC checks passed.")

    _coverage(cfg, edges, nodes, reports)
    return 0


def _coverage(cfg, edges, nodes, reports):
    import os
    import networkx as nx
    import pandas as pd

    tmap = dict(zip(nodes.node_id, nodes.node_type))
    G = nx.Graph()
    G.add_nodes_from(nodes.node_id)
    G.add_edges_from(zip(edges.node1, edges.node2))
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    largest = comps[0] if comps else set()

    # protein neighbour (direct or via lipid class)
    def protein_reach(n):
        for nb in G.neighbors(n):
            if tmap.get(nb) == "protein":
                return True
            if tmap.get(nb) == "lipid_class":
                if any(tmap.get(x) == "protein" for x in G.neighbors(nb)):
                    return True
        return False

    meas = nodes[nodes.is_measured == True].copy()  # noqa: E712
    meas["deg"] = meas["node_id"].map(dict(pd.concat([edges.node1, edges.node2]).value_counts()))
    meas["deg"] = meas["deg"].fillna(0).astype(int)
    meas["has_edge"] = meas["deg"] > 0
    meas["protein_reach"] = meas["node_id"].apply(lambda n: protein_reach(n) if n in G else False)
    meas["in_largest"] = meas["node_id"].isin(largest)

    L = ["# mnet final coverage", ""]
    L.append(f"- nodes: {len(nodes):,} {nodes.node_type.value_counts().to_dict()}")
    L.append(f"- edges: {len(edges):,} {edges.edge_type.value_counts().to_dict()}")
    L.append(f"- components: {len(comps):,} | largest component: {len(largest):,} nodes")
    L.append("")
    L.append("## Measured metabolites by super_class")
    L.append("")
    L.append("| super_class | total | >=1 edge | protein-reachable | in largest comp |")
    L.append("|---|---|---|---|---|")
    g = meas.groupby(meas.super_class.fillna("(lipid/none)"))
    for sc, sub in g:
        L.append(f"| {sc} | {len(sub)} | {int(sub.has_edge.sum())} | "
                 f"{int(sub.protein_reach.sum())} | {int(sub.in_largest.sum())} |")
    L.append(f"| **TOTAL** | {len(meas)} | {int(meas.has_edge.sum())} | "
             f"{int(meas.protein_reach.sum())} | {int(meas.in_largest.sum())} |")
    L.append("")
    L.append("## Proteins")
    L.append(f"- protein nodes: {int((nodes.node_type=='protein').sum()):,}; "
             f"Rhea enzymes new (in_string=False): "
             f"{int(((nodes.node_type=='protein') & (nodes.in_string==False)).sum()):,}")
    L.append("")
    L.append("## Lipids by lipid_class (measured)")
    L.append(f"{meas[meas.lipid_class.notna()].lipid_class.value_counts().to_dict()}")
    L.append("")
    # hubs
    deg_all = pd.concat([edges.node1, edges.node2]).value_counts()
    L.append("## Top 10 hubs per node type")
    for t in ["protein", "metabolite", "lipid_class"]:
        ids = [n for n in deg_all.index if tmap.get(n) == t][:10]
        L.append(f"- {t}: " + ", ".join(f"{n}({deg_all[n]})" for n in ids))
    L.append("")
    # measured with no protein connection + reason
    # membership in human Rhea participants (for reason categorisation)
    import os as _os
    parts = pd.read_parquet(_os.path.join(_os.path.dirname(_os.path.dirname(
        __file__)), "..", "data", "interim", "rhea_participants.parquet")) \
        if False else None
    try:
        from . import metabolites as MB
        parts = pd.read_parquet(_os.path.join(MB._interim(cfg), "rhea_participants.parquet"))
        rhea_chebi = set(parts["chebi_id"])
    except Exception:
        rhea_chebi = set()

    # is_a relations for parent/sibling suggestions
    from . import metabolites as MBc
    rel = pd.read_parquet(os.path.join(MBc._interim(cfg), "chebi_relations.parquet"))
    isa = rel[rel.relation == "is_a"]
    up, down = {}, {}
    for s, o in zip(isa.subject, isa.object):
        up.setdefault(s, set()).add(o)
        down.setdefault(o, set()).add(s)
    nm = pd.read_parquet(os.path.join(MBc._interim(cfg), "chebi_names.parquet")).set_index("chebi_id")["name"].to_dict()
    deg_all = pd.concat([edges.node1, edges.node2]).value_counts().to_dict()
    edged_chebi = {n for n in deg_all if str(n).startswith("CHEBI:") and deg_all[n] > 0}

    noconn = meas[~meas.protein_reach]
    # suggestions: degree-0 measured ChEBI that is a parent/sibling of an edged node
    suggestions = []
    for _, r in noconn.iterrows():
        cr = str(r.get("chebi_rhea") or r.get("chebi_id") or "")
        for c in cr.split(";"):
            if not c.startswith("CHEBI:"):
                continue
            cands = set(down.get(c, set()))                 # children (c is parent)
            for p in up.get(c, set()):
                cands |= down.get(p, set())                  # siblings
            edged = sorted((x for x in cands if x in edged_chebi and x != c),
                           key=lambda x: -deg_all.get(x, 0))
            if edged:
                t = edged[0]
                suggestions.append((r["node_id"], r.get("label"), t, nm.get(t, "?"), deg_all.get(t, 0)))
                break
    L.append(f"## Suggested targets for degree-0 metabolites ({len(suggestions)}) - NOT auto-applied")
    for nid, lab, t, tn, d in suggestions:
        L.append(f"- {nid} ({lab}) -> {t} ({tn}, degree {d})")
    L.append("")

    L.append(f"## Measured metabolites with NO protein connection ({len(noconn)})")
    for _, r in noconn.head(60).iterrows():
        cr = str(r.get("chebi_rhea") or r.get("chebi_id") or "")
        in_rhea = any(c in rhea_chebi for c in cr.split(";") if c)
        if r.get("mapping_method") == "unmapped":
            reason = "unmapped"
        elif not cr or cr == "nan":
            reason = "no ChEBI mapping"
        elif in_rhea:
            reason = "in Rhea but only as currency/cofactor (edges dropped)"
        else:
            reason = "not in human Rhea (diet/xenobiotic or no human enzyme)"
        L.append(f"- {r['node_id']} ({r.get('label')}): {reason}")
    L.append("")
    # unmeasured evidence: largest-comp coverage with vs without unmeasured metabolites
    unmeas = set(nodes[(nodes.node_type == "metabolite") & (nodes.is_measured == False)]["node_id"])  # noqa: E712
    G2 = G.copy()
    G2.remove_nodes_from(unmeas)
    comps2 = sorted(nx.connected_components(G2), key=len, reverse=True)
    largest2 = comps2[0] if comps2 else set()
    cov_with = meas.in_largest.mean()
    cov_without = meas.node_id.isin(largest2).mean()
    L.append("## Evidence for including unmeasured metabolites")
    L.append(f"- measured metabolites in largest component WITH unmeasured: {100*cov_with:.1f}%")
    L.append(f"- WITHOUT unmeasured metabolites: {100*cov_without:.1f}%")
    L.append("")
    with open(os.path.join(reports, "mnet_final_coverage.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"[mnet.qc] wrote {os.path.join(reports, 'mnet_final_coverage.md')}")
