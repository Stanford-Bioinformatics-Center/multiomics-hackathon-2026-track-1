"""Rhea reactions: enzymes, directions, cross-refs and ChEBI participants.

Participants are fetched from the Rhea SPARQL endpoint rather than parsed from
rhea.rdf.gz: the full RDF/XML is large and rdflib parsing is slow and
memory-heavy, whereas the endpoint returns the ~83k participant rows in a few
paged queries. (Documented choice, per the phase brief.)

compound_type (small molecule | generic | polymer | other) is not carried on the
Rhea compound node, so it is derived from the ChEBI molecular formula: R-group /
wildcard placeholders -> generic, polymer repeat notation -> polymer, a concrete
formula -> small molecule, missing -> other.
"""
from __future__ import annotations

import re
import time
from typing import Dict, Optional

import pandas as pd

RHEA_NS = "http://rdf.rhea-db.org/"
CHEBI_PURL = "http://purl.obolibrary.org/obo/CHEBI_"


# --- flat TSV parsers -------------------------------------------------------

def parse_enzymes(path: str) -> pd.DataFrame:
    """rhea2uniprot_sprot.tsv -> RHEA_ID, DIRECTION, MASTER_ID, uniprot."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    ren = {"ID": "uniprot"}
    df = df.rename(columns=ren)
    keep = ["RHEA_ID", "DIRECTION", "MASTER_ID", "uniprot"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df["MASTER_ID"] = pd.to_numeric(df["MASTER_ID"], errors="coerce").astype("Int64")
    return df.dropna(subset=["uniprot"]).drop_duplicates().reset_index(drop=True)


def parse_directions(path: str) -> pd.DataFrame:
    """rhea-directions.tsv -> master / LR / RL / BI reaction ids."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df.reset_index(drop=True)


def parse_xrefs(path: str) -> pd.DataFrame:
    """rhea2xrefs.tsv -> RHEA_ID, DIRECTION, MASTER_ID, ID, DB (EC/KEGG/...)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df.reset_index(drop=True)


# --- ChEBI pH 7.3 mapping ---------------------------------------------------

def parse_ph7_3_mapping(path: str) -> pd.DataFrame:
    """chebi_pH7_3_mapping.tsv -> CHEBI, CHEBI_PH7_3 (Rhea's charged form)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip().upper() for c in df.columns]
    return df.reset_index(drop=True)


# --- participants via SPARQL ------------------------------------------------

def _short_chebi(uri: str) -> str:
    return "CHEBI:" + uri.rsplit("CHEBI_", 1)[-1] if "CHEBI_" in uri else uri


def _master_id(uri: str) -> Optional[int]:
    m = re.search(r"/(\d+)$", uri)
    return int(m.group(1)) if m else None


def _side(uri: str) -> str:
    return "R" if uri.endswith("_R") else ("L" if uri.endswith("_L") else "?")


# Rhea compound class -> our compound_type (Rhea is the source of truth).
RHEA_CLASS_MAP = {
    "SmallMolecule": "small molecule",
    "Polymer": "polymer",
    "GenericPolypeptide": "generic",
    "GenericPolynucleotide": "generic",
    "GenericCompound": "generic",
    "ReactivePart": "generic",
}


def _class_name(uri: str) -> str:
    return uri.rsplit("/", 1)[-1]


def fetch_participants(sparql_url: str, page: int = 20000,
                       max_pages: int = 100) -> pd.DataFrame:
    """Fetch participants with Rhea compound class + reactive/underlying ChEBI.

    Returns master_rhea_id, side, chebi_id (the associated ChEBI: the small
    molecule, or the reactive part / underlying unit for generics / polymers),
    compound_type (from the Rhea class) and underlying_chebi (the reactive-part /
    underlying-unit ChEBI, null for small molecules).
    """
    from SPARQLWrapper import JSON, SPARQLWrapper

    s = SPARQLWrapper(sparql_url, agent="mnet/0.1 (research)")
    s.setReturnFormat(JSON)
    q = """PREFIX rh: <http://rdf.rhea-db.org/>
SELECT ?rhea ?side ?class ?chebi ?rpchebi ?uchebi WHERE {
  ?rhea rdfs:subClassOf rh:Reaction .
  ?rhea rh:side ?side .
  ?side rh:contains ?p .
  ?p rh:compound ?c .
  ?c rdfs:subClassOf ?class . FILTER(STRSTARTS(STR(?class), "http://rdf.rhea-db.org/"))
  OPTIONAL { ?c rh:chebi ?chebi }
  OPTIONAL { ?c rh:reactivePart ?rp . ?rp rh:chebi ?rpchebi }
  OPTIONAL { ?c rh:underlyingChebi ?uchebi }
} ORDER BY ?rhea ?side ?class ?chebi LIMIT %d OFFSET %d"""

    rows = []
    for pg in range(max_pages):
        s.setQuery(q % (page, pg * page))
        res = s.query().convert()["results"]["bindings"]
        if not res:
            break
        for r in res:
            cls = _class_name(r["class"]["value"])
            chebi = _short_chebi(r["chebi"]["value"]) if "chebi" in r else None
            rp = _short_chebi(r["rpchebi"]["value"]) if "rpchebi" in r else None
            uc = _short_chebi(r["uchebi"]["value"]) if "uchebi" in r else None
            assoc = chebi or rp or uc
            underlying = rp or uc
            if assoc is None:
                continue
            rows.append((_master_id(r["rhea"]["value"]), _side(r["side"]["value"]),
                         assoc, RHEA_CLASS_MAP.get(cls, "other"), underlying))
        if len(res) < page:
            break
        time.sleep(0.2)
    df = pd.DataFrame(rows, columns=["master_rhea_id", "side", "chebi_id",
                                     "compound_type", "underlying_chebi"])
    # a compound may match >1 class row; keep the most specific (non-small) type
    order = {"polymer": 0, "generic": 1, "small molecule": 2, "other": 3}
    df["_o"] = df["compound_type"].map(order).fillna(3)
    df = (df.sort_values("_o").drop_duplicates(["master_rhea_id", "side", "chebi_id"])
          .drop(columns="_o"))
    return df.reset_index(drop=True)


def classify_compound_type(formula: Optional[str]) -> str:
    """Legacy formula-based heuristic (kept for the Phase 1 vs Phase 2 comparison)."""
    if formula is None or (isinstance(formula, float)) or str(formula).strip() in ("", "nan", "."):
        return "other"
    f = str(formula)
    if re.search(r"\)n\b", f) or re.search(r"n$", f) or "*" in f:
        return "polymer" if ("n" in f and "*" not in f) else "generic"
    if re.search(r"\bR\b", f) or "R" in re.sub(r"[A-QS-Za-qs-z]", "", f) or "X" in f:
        return "generic"
    return "small molecule"


def build_participants(sparql_url: str, chebi_formula: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Full participants table (Rhea-class compound_type) with is_transport flag."""
    df = fetch_participants(sparql_url)
    sides = df.groupby(["master_rhea_id", "chebi_id"])["side"].nunique()
    trans = {mid for (mid, _c) in set(sides[sides > 1].index)}
    df["is_transport"] = df["master_rhea_id"].isin(trans)
    return df


def get_protein_generic_chebis(sparql_url: str) -> set:
    """ChEBI ids that are GenericPolypeptide / GenericPolynucleotide reactive parts
    (e.g. 'L-seryl-[protein]') - these are proteins, not metabolites, and are
    excluded from protein-metabolite edges."""
    from SPARQLWrapper import JSON, SPARQLWrapper
    s = SPARQLWrapper(sparql_url, agent="mnet/0.1")
    s.setReturnFormat(JSON)
    s.setQuery("""PREFIX rh: <http://rdf.rhea-db.org/>
SELECT DISTINCT ?chebi WHERE {
  ?c rdfs:subClassOf ?cls .
  FILTER(?cls IN (rh:GenericPolypeptide, rh:GenericPolynucleotide))
  ?c rh:reactivePart ?rp . ?rp rh:chebi ?chebi .
}""")
    return {_short_chebi(r["chebi"]["value"]) for r in s.query().convert()["results"]["bindings"]}


def run(cfg, args) -> int:
    import os
    import pandas as pd
    from . import metabolites as MB
    from . import network as net
    from . import proteins as P

    mcfg = cfg.get("mnet", {}) or {}
    interim = MB._interim(cfg)
    reports = cfg["paths"]["reports"]
    curation = os.path.join(cfg["_project_root"], "curation")
    ctx = MB.Ctx(cfg)

    parts = MB._read(cfg, "rhea_participants")
    enz = MB._read(cfg, "rhea_enzymes")
    p2 = MB._read(cfg, "metabolites_nonlipid")
    lip = MB._read(cfg, "lipids")

    # ---- human enzymes + reactions ----
    human = ctx._human_set()
    enz_h = enz[enz["uniprot"].isin(human)].copy()
    enz_h["master_rhea_id"] = enz_h["MASTER_ID"].astype("Int64")
    human_reacs = set(enz_h["master_rhea_id"].dropna().astype(int))

    # ---- 4A crosswalk ----
    accs = sorted(set(enz_h["uniprot"].dropna()))
    cross = P.map_accessions_to_nodes(cfg, accs)
    cross.to_parquet(os.path.join(interim, "protein_crosswalk.parquet"), index=False)
    acc2node = dict(zip(cross["accession"], cross["node_id"]))
    print(f"[mnet.rhea] crosswalk: {len(cross):,} accessions, "
          f"new nodes (in_string=False): {int((~cross['in_string']).sum()):,}")

    # ---- 4B participants for human reactions ----
    # keep the FULL participant list (incl GenericPolypeptide/GenericPolynucleotide)
    # for the side/currency analysis; protein-generics count as non-currency real
    # participants but never become metabolite edges.
    ph = parts[parts["master_rhea_id"].isin(human_reacs)].copy()
    prot_generic = get_protein_generic_chebis(cfg["urls"]["rhea_sparql"])
    # a ChEBI that also occurs as a FREE small molecule (e.g. FAD, which is also a
    # reactive part of flavoproteins) must NOT be excluded as a protein-generic.
    small_mol = set(parts.loc[parts["compound_type"] == "small molecule", "chebi_id"])
    prot_generic = prot_generic - small_mol
    chebi_ctype = dict(zip(ph["chebi_id"], ph["compound_type"]))

    measured_chebi = set()
    for v in p2["chebi_rhea"].dropna():
        measured_chebi |= {c for c in str(v).split(";") if c}
    measured_chebi |= {str(c) for c in lip["exact_chebi"].dropna()}

    def eligible_metabolite(c):
        if c in prot_generic:
            return False
        ct = chebi_ctype.get(c)
        if ct in ("small molecule", "generic"):
            return True
        if ct == "polymer" and c in measured_chebi:
            return True
        return False

    # ---- 4C currency ----
    currency = set(ctx.currency)
    hub_ids = {c for c in ctx.human_rhea_chebi
               if ctx.human_rxn_count.get(c, 0) > mcfg.get("currency_max_reactions", 150)}
    pd.DataFrame([{"chebi_id": c, "name": ctx.chebi_name.get(c, "?"),
                   "n_human_reactions": ctx.human_rxn_count.get(c, 0)}
                  for c in sorted(hub_ids, key=lambda c: -ctx.human_rxn_count.get(c, 0))]
                 ).to_csv(os.path.join(reports, "mnet_phase4_auto_currency.csv"), index=False)
    currency_all = currency | hub_ids

    pairs = pd.read_csv(os.path.join(curation, "currency_pairs.csv"), dtype=str)
    partner: dict = {}
    for _, r in pairs.iterrows():
        partner.setdefault(r["chebi_a"], set()).add(r["chebi_b"])
        partner.setdefault(r["chebi_b"], set()).add(r["chebi_a"])

    keep_rxns = set(pd.read_csv(os.path.join(curation, "currency_keep_reactions.csv"),
                                dtype=str)["rhea_id"].astype(int))

    def kept_participants(policy: str) -> pd.DataFrame:
        rows = []
        for rxn, sub in ph.groupby("master_rhea_id"):
            Lall = set(sub[sub.side == "L"]["chebi_id"])
            Rall = set(sub[sub.side == "R"]["chebi_id"])
            real_L = {c for c in Lall if c not in currency_all}  # protein-generics count as real
            real_R = {c for c in Rall if c not in currency_all}
            has_real_each = bool(real_L) and bool(real_R)
            no_real_any = (not real_L) and (not real_R)        # pure hydrolysis / cofactor-only
            keep_rxn = int(rxn) in keep_rxns
            for c, s in zip(sub["chebi_id"], sub["side"]):
                if not eligible_metabolite(c):
                    continue                            # protein-generics / ineligible: no edge
                if c in currency_all:
                    measured = c in measured_chebi
                    if not measured:
                        continue                        # unmeasured currency: always drop
                    if keep_rxn or policy == "keep_all":
                        pass                            # explicit keep
                    elif policy == "drop":
                        continue
                    else:                               # cofactor_rule
                        opp = Rall if s == "L" else Lall
                        cofactor = bool(partner.get(c, set()) & opp) and has_real_each
                        if cofactor or no_real_any:
                            continue
                rows.append((rxn, c))
        return pd.DataFrame(rows, columns=["master_rhea_id", "chebi_id"]).drop_duplicates()

    policy = mcfg.get("measured_currency_policy", "cofactor_rule")

    # ---- metabolite node resolution ----
    lc = pd.read_csv(os.path.join(curation, "lipid_class_to_chebi.csv"), dtype=str)
    chebi2class: dict = {}
    for _, r in lc.iterrows():
        if isinstance(r["chebi_id"], str) and r["chebi_id"]:
            shared = str(r.get("shared_generic", "")).lower() in ("true", "1")
            chebi2class.setdefault(r["chebi_id"], []).append(
                (f"LIPIDCLASS:{r['class']}", shared))
    # chebi -> set of measured node ids claiming it (a chebi can back several
    # measured metabolites, e.g. Leucine, and the combined Leucine/Isoleucine).
    chebi2meas: dict = {}
    for _, r in p2.iterrows():
        if pd.notna(r["chebi_rhea"]):
            for c in str(r["chebi_rhea"]).split(";"):
                if c:
                    chebi2meas.setdefault(c, set()).add(r["node_id"])
    for _, r in lip.iterrows():
        if pd.notna(r["exact_chebi"]):
            chebi2meas.setdefault(str(r["exact_chebi"]), set()).add(r["node_id"])

    def resolve(chebi):
        if chebi in chebi2class:
            return [(cn, "lipid_class", False, sh) for cn, sh in chebi2class[chebi]]
        if chebi in chebi2meas:
            return [(nid, "metabolite", True, False) for nid in sorted(chebi2meas[chebi])]
        return [(chebi, "metabolite", False, False)]

    rxn_transport = ph.groupby("master_rhea_id")["is_transport"].any().to_dict()
    enz_h["node"] = enz_h["uniprot"].map(acc2node)
    enz_map = enz_h[["master_rhea_id", "node"]].dropna().drop_duplicates()

    def build_edges(kept: pd.DataFrame):
        metab_rows = []
        for rxn, c in zip(kept["master_rhea_id"], kept["chebi_id"]):
            for node, ntype, ismeas, shared in resolve(c):
                metab_rows.append((rxn, node, ntype, ismeas, shared))
        md = pd.DataFrame(metab_rows, columns=["master_rhea_id", "mnode",
                                               "mtype", "is_measured", "shared_generic"]).drop_duplicates()
        e = enz_map.merge(md, on="master_rhea_id")
        e = e[e["node"] != e["mnode"]]
        e["edge_type"] = e["master_rhea_id"].map(rxn_transport).map(
            lambda t: "transport" if t else "catalysis")
        return e, md

    # counts under all 3 policies (for report)
    policy_counts = {}
    for pol in ("cofactor_rule", "drop", "keep_all"):
        k = kept_participants(pol)
        e, _ = build_edges(k)
        agg = e.groupby(["node", "mnode", "edge_type"]).size().reset_index()
        policy_counts[pol] = len(agg)

    # build the configured policy for output
    kept = kept_participants(policy)
    e, md = build_edges(kept)
    e["rhea"] = "RHEA:" + e["master_rhea_id"].astype(str)
    edges = (e.groupby(["node", "mnode", "edge_type"])
             .agg(evidence=("rhea", lambda s: ";".join(sorted(set(s)))),
                  n_evidence=("rhea", lambda s: s.nunique()),
                  node2_type=("mtype", "first"),
                  is_measured=("is_measured", "first"),
                  shared_generic=("shared_generic", "any"))
             .reset_index()
             .rename(columns={"node": "node1", "mnode": "node2"}))
    edges["node1_type"] = "protein"

    edges.to_parquet(os.path.join(interim, "edges_rhea.parquet"), index=False)

    # ---- currency degree report (before/after) + ATP>300 stop check ----
    import sys
    cur_ids = {"ATP": "CHEBI:30616", "ADP": "CHEBI:456216", "AMP": "CHEBI:456215",
               "GTP": "CHEBI:37565", "NAD+": "CHEBI:57540", "FAD": "CHEBI:57692"}
    before = {"ATP": 871, "ADP": 782, "AMP": 120, "GTP": 244, "NAD+": 107, "FAD": 0}
    cur_deg = {n: int((edges["node2"] == cid).sum()) for n, cid in cur_ids.items()}
    cur_ex = {}
    for n, cid in cur_ids.items():
        ev = set()
        for v in edges[edges["node2"] == cid]["evidence"].dropna():
            ev |= set(str(v).split(";"))
        cur_ex[n] = sorted(ev)[:5]
    print(f"[mnet.rhea] currency degrees after cofactor rule: {cur_deg}")
    if cur_deg["ATP"] > 300:
        top = (edges[edges["node2"] == cur_ids["ATP"]]
               .assign(d=1).groupby("node1").size().sort_values(ascending=False).head(15))
        print("[mnet.rhea] ATP still > 300 protein edges; top enzymes:")
        print(top.to_string())
        sys.exit(2)

    # ---- nodes_rhea ----
    prot_nodes = cross[["node_id", "in_string"]].drop_duplicates().copy()
    prot_nodes["node_type"] = "protein"
    metab_nodes = md[["mnode", "mtype", "is_measured"]].drop_duplicates().rename(
        columns={"mnode": "node_id", "mtype": "node_type"})
    metab_nodes["chebi_id"] = metab_nodes["node_id"].where(
        metab_nodes["node_id"].str.startswith("CHEBI:"))
    metab_nodes["label"] = metab_nodes["chebi_id"].map(ctx.chebi_name)
    nodes_rhea = pd.concat([prot_nodes, metab_nodes], ignore_index=True)
    nodes_rhea.to_parquet(os.path.join(interim, "nodes_rhea.parquet"), index=False)

    # ---- report + checkpoint ----
    _rhea_report(cfg, edges, cross, p2, md, policy_counts, policy, reports,
                 before, cur_deg, cur_ex)

    print(f"[mnet.rhea] edges={len(edges):,} policy={policy} "
          f"counts={policy_counts}")
    _checkpoint(edges, cross, p2, currency_all)
    return 0


def _checkpoint(edges, cross, p2, currency_all):
    import sys
    # measured non-lipid metabolites with >=1 protein edge
    meas_nodes = set(p2[p2["mapping_method"] != "unmapped"]["node_id"])
    prot_edges = edges  # all are protein-metabolite
    connected = set(prot_edges["node2"]) & meas_nodes
    frac = len(connected) / max(len(meas_nodes), 1)
    # the >500 hub check targets a *missed* common metabolite, so known currency
    # (incl. measured-currency pairs like ATP/ADP handled by the cofactor rule)
    # is excluded - it is expected to be high-degree.
    deg = prot_edges.groupby("node2").size()
    deg_noncur = deg[~deg.index.isin(currency_all)]
    max_deg = int(deg_noncur.max()) if len(deg_noncur) else 0
    new_frac = (~cross["in_string"]).mean()
    fails = []
    if frac < 0.50:
        fails.append(f"only {100*frac:.1f}% of measured non-lipid metabolites have a protein edge (<50%)")
    if max_deg > 500:
        top = deg_noncur.idxmax()
        fails.append(f"non-currency metabolite {top} has {max_deg} protein edges (>500)")
    if new_frac > 0.05:
        fails.append(f"{100*new_frac:.1f}% of Rhea enzymes became new nodes (>5%)")
    print(f"[mnet.rhea] CHECKPOINT: measured-connected={100*frac:.1f}% "
          f"max_metab_degree={max_deg} new_node_frac={100*new_frac:.1f}%")
    if fails:
        print("[mnet.rhea] CHECKPOINT FAILED:")
        for f in fails:
            print("   -", f)
        sys.exit(2)


def _rhea_report(cfg, edges, cross, p2, md, policy_counts, policy, reports,
                 before=None, cur_deg=None, cur_ex=None):
    import os
    L = ["# mnet Phase 4 - Rhea protein-metabolite layer", ""]
    L.append(f"- edges (configured policy `{policy}`): {len(edges):,}")
    L.append(f"  - catalysis: {int((edges.edge_type=='catalysis').sum()):,} | "
             f"transport: {int((edges.edge_type=='transport').sum()):,}")
    L.append(f"- edge counts by measured_currency_policy: {policy_counts}")
    if cur_deg:
        L.append("")
        L.append("## Currency protein-edge counts: before -> after (cofactor rule)")
        L.append("")
        L.append("| metabolite | before | after |")
        L.append("|---|---|---|")
        for n in ["ATP", "ADP", "AMP", "GTP", "NAD+", "FAD"]:
            L.append(f"| {n} | {before.get(n,'-')} | {cur_deg.get(n,0)} |")
        L.append("")
        L.append("### Example kept reactions")
        for n in ["ATP", "ADP", "AMP", "GTP", "NAD+", "FAD"]:
            L.append(f"- {n}: {cur_ex.get(n, [])}")
    L.append(f"- protein nodes: {cross['node_id'].nunique():,} "
             f"(new/in_string=False: {int((~cross['in_string']).sum()):,})")
    L.append(f"- crosswalk methods: {cross['method'].value_counts().to_dict()}")
    L.append(f"- metabolite/class nodes with an edge: {edges['node2'].nunique():,} "
             f"(measured: {int(edges['is_measured'].sum()):,} edges)")
    L.append(f"- shared_generic edges: {int(edges['shared_generic'].sum()):,}")
    with open(os.path.join(reports, "mnet_phase4_rhea.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
