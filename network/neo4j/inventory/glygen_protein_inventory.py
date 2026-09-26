#!/usr/bin/env python3
# =====================================================================================================
# network/neo4j/inventory/glygen_protein_inventory.py — WHAT DOES GLYGEN HOLD FOR EACH OF OUR 471 PROTEINS?
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Before choosing which GlyGen data to add to the Neo4j graph, count, for each of our 471 proteins, how
#   many records GlyGen (https://www.glygen.org, Homo sapiens) has in every data section: glycosylation
#   (reported, with glycans, predicted), phosphorylation, glycation, mutations, expression, disease,
#   pathways, reactions, structures, and so on. Nothing is bulk-downloaded: one API call per protein.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Maps our genes (Entrez ID, symbol, UniProt) to GlyGen's canonical UniProt accession, using GlyGen's
#      small Entrez cross-reference and masterlist files.
#   2. Calls https://api.glygen.org/protein/detail/<accession>/ for each protein (6 at a time) and keeps only
#      counts: the length of every data section, GlyGen's own section statistics (e.g. N-linked sites with
#      glycans), and a few site-level details (unique phosphosites, sites with a kinase, glycosylation sites
#      by type and evidence). The full responses are not stored.
#   3. Writes one row per protein and a long list of GlyGen phosphosites (for matching MoTrPAC sites).
#
# HOW TO RUN
#   After the pipeline (steps 1-2):   python3 network/neo4j/inventory/glygen_protein_inventory.py
#   (about 5 minutes; standard library only; re-running skips proteins already fetched)
#
# DATA AND PROVENANCE
#   GlyGen data release 2.11.1 (API https://api.glygen.org; mapping files from
#   https://data.glygen.org/ln2data/releases/data/current/reviewed/). Our genes: $HACK_OUT/02_nodes_string.csv.
#
# INPUTS:  $HACK_OUT/02_nodes_string.csv (entrez_gene, gene_symbol, uniprot)
# OUTPUTS: $HACK_OUT/inventory/glygen_protein_counts.csv   one row per protein, one column per count
#          $HACK_OUT/inventory/glygen_phosphosites.csv      entrez_gene, glygen_ac, position, residue, kinase
#          $HACK_OUT/inventory/glygen_cache/                small per-protein count files (resume support)
# =====================================================================================================

import csv, io, json, os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

# Where the pipeline's tables are, and where the inventory goes (outside the repo).
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
INV = os.path.join(OUT, "inventory"); CACHE = os.path.join(INV, "glygen_cache")
os.makedirs(CACHE, exist_ok=True)
# GlyGen endpoints.
DATA = "https://data.glygen.org/ln2data/releases/data/current/reviewed/"
API = "https://api.glygen.org/protein/detail/{}/"

def get(url, tries=4):
    """Download a URL (with a few retries), returning the text."""
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "motrpac-hackathon-inventory"}), timeout=120) as r:
                return r.read().decode("utf-8")
        except Exception as e:
            if i == tries - 1: raise
            time.sleep(3 * (i + 1))

def read_csv_text(text):
    """Parse CSV text into a list of dicts."""
    return list(csv.DictReader(io.StringIO(text)))

# ---- 1. our genes -> GlyGen canonical accession --------------------------------------------------------
genes = read_csv_text(open(os.path.join(OUT, "02_nodes_string.csv")).read())
# GlyGen's Entrez cross-references (canonical accession <-> Entrez) and masterlist (all isoforms per canonical).
xref = read_csv_text(get(DATA + "human_protein_xref_geneid.csv"))
master = read_csv_text(get(DATA + "human_protein_masterlist.csv"))
by_entrez = {}
for r in xref: by_entrez.setdefault(r["xref_id"], []).append(r["uniprotkb_canonical_ac"])
by_acc = {}
for r in master:
    for iso in (r["reviewed_isoforms"] + "|" + r["unreviewed_isoforms"] + "|" + r["uniprotkb_canonical_ac"]).replace(";", "|").split("|"):
        if iso: by_acc[iso.split("-")[0]] = r["uniprotkb_canonical_ac"]
mapping = []
for g in genes:
    cands = by_entrez.get(g["entrez_gene"], [])
    # prefer the canonical whose accession matches our UniProt; else the first Entrez match; else our UniProt
    acc_ours = (g.get("uniprot") or "").split("-")[0]
    ac = next((c for c in cands if c.split("-")[0] == acc_ours), cands[0] if cands else by_acc.get(acc_ours))
    how = "entrez+uniprot" if cands and ac.split("-")[0] == acc_ours else ("entrez" if cands else ("uniprot" if ac else "unmapped"))
    mapping.append({"entrez_gene": g["entrez_gene"], "gene_symbol": g["gene_symbol"], "uniprot": g.get("uniprot", ""), "glygen_ac": ac or "", "mapped_via": how})
print(f"mapped {sum(1 for m in mapping if m['glygen_ac'])} of {len(mapping)} genes to GlyGen accessions", file=sys.stderr)

# ---- 2. per-protein counts from the API ------------------------------------------------------------------
def summarise(d):
    """Reduce one GlyGen protein record to counts (plus the phosphosite list)."""
    c = {}
    # every data section: number of records
    for k, v in d.items():
        if isinstance(v, list): c[f"n_{k}"] = len(v)
        elif isinstance(v, dict) and k in ("go_annotation",):
            c["n_go_terms"] = sum(len(cat.get("go_terms", [])) for cat in v.get("categories", []))
    # GlyGen's own statistics per table (glycosylation_reported, _predicted, ...)
    for t in d.get("section_stats", []):
        for s in t.get("table_stats", []):
            c[f"{t['table_id']}__{s['field']}"] = s["count"]
    # glycosylation. Records without a position are protein-level evidence (mostly the O-GlcNAc Database:
    # "this protein is O-GlcNAcylated, site not known"); they are counted separately from sites.
    allgly = d.get("glycosylation", [])
    nosite = [g for g in allgly if not g.get("start_pos")]
    c["gly_protein_level_no_site"] = int(len(nosite) > 0)
    c["gly_protein_level_oglcnac"] = int(any((g.get("subtype") or "") == "O-GlcNAcylation" for g in nosite))
    c["gly_oglcnac_sites"] = len({g.get("start_pos") for g in allgly if g.get("start_pos") and (g.get("subtype") or "") == "O-GlcNAcylation"})
    gly = [g for g in allgly if g.get("start_pos")]
    sites = {(g.get("type"), g.get("start_pos")) for g in gly if g.get("start_pos")}
    c["gly_sites_unique"] = len(sites)
    for ty in ("N-linked", "O-linked", "C-linked"):
        c[f"gly_sites_{ty.split('-')[0]}"] = len({p for t, p in sites if t == ty})
    cats = {}
    for g in gly:
        for cat in (g.get("site_category_dict") or {g.get("site_category"): True}):
            cats.setdefault(cat, set()).add(g.get("start_pos"))
    for cat, s in cats.items(): c[f"gly_sites_{cat}"] = len(s)
    # glycan structures observed at a site (excluding the generic O-GlcNAc monosaccharide on protein-level records)
    c["glycans_at_sites"] = len({g["glytoucan_ac"] for g in gly if g.get("glytoucan_ac")})
    # phosphorylation: unique sites and sites with a kinase
    ph = d.get("phosphorylation", [])
    c["phospho_sites_unique"] = len({p.get("start_pos") for p in ph})
    c["phospho_sites_with_kinase"] = len({p.get("start_pos") for p in ph if p.get("kinase_uniprot_canonical_ac") or p.get("kinase_gene_name")})
    sites_list = [{"position": p.get("start_pos"), "residue": p.get("residue"),
                   "kinase": p.get("kinase_gene_name") or p.get("kinase_uniprot_canonical_ac") or "",
                   "source": ";".join(sorted({e.get("database", "") for e in p.get("evidence", []) if e.get("database") != "PubMed"}))} for p in ph]
    # glycation and other modifications, if present
    c["glycation_sites_unique"] = len({p.get("start_pos") for p in d.get("glycation", [])})
    return c, sites_list

def fetch(m):
    """Fetch (or read from cache) one protein's counts."""
    ac = m["glygen_ac"]
    if not ac: return m, None, []
    f = os.path.join(CACHE, ac + ".json")
    if os.path.exists(f):
        j = json.load(open(f)); return m, j["counts"], j["phospho"]
    try:
        c, ps = summarise(json.loads(get(API.format(ac))))
    except Exception as e:
        print(f"  failed {ac}: {e}", file=sys.stderr); return m, None, []
    json.dump({"counts": c, "phospho": ps}, open(f, "w"))
    return m, c, ps

with ThreadPoolExecutor(max_workers=6) as ex:
    results = list(ex.map(fetch, mapping))

# ---- 3. write the tables ---------------------------------------------------------------------------------
cols = sorted({k for _, c, _ in results if c for k in c})
with open(os.path.join(INV, "glygen_protein_counts.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["entrez_gene", "gene_symbol", "uniprot", "glygen_ac", "mapped_via", "fetched"] + cols)
    for m, c, _ in results:
        w.writerow([m["entrez_gene"], m["gene_symbol"], m["uniprot"], m["glygen_ac"], m["mapped_via"], int(c is not None)] + [(c or {}).get(k, 0) for k in cols])
with open(os.path.join(INV, "glygen_phosphosites.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["entrez_gene", "gene_symbol", "glygen_ac", "position", "residue", "kinase", "source"])
    for m, _, ps in results:
        for p in ps: w.writerow([m["entrez_gene"], m["gene_symbol"], m["glygen_ac"], p["position"], p["residue"], p["kinase"], p["source"]])
print(f"fetched {sum(1 for _, c, _ in results if c)} proteins -> {INV}", file=sys.stderr)
