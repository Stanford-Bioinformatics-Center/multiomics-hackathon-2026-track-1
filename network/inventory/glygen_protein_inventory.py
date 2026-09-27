#!/usr/bin/env python3
# =====================================================================================================
# network/inventory/glygen_protein_inventory.py — WHAT DOES GLYGEN HOLD FOR EACH OF OUR 471 PROTEINS?
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
#   After the pipeline (steps 1-2):   python3 network/inventory/glygen_protein_inventory.py
#   (about 5 minutes; standard library only; re-running skips proteins already fetched)
#
# DATA AND PROVENANCE
#   GlyGen data release 2.11.1 (API https://api.glygen.org; mapping files from
#   https://data.glygen.org/ln2data/releases/data/current/reviewed/). Our genes: $HACK_OUT/02_nodes_string.csv.
#
# INPUTS:  $HACK_OUT/02_nodes_string.csv (entrez_gene, gene_symbol, uniprot)
# OUTPUTS: $HACK_OUT/inventory/glygen_protein_counts.csv   one row per protein, one column per count
#          $HACK_OUT/inventory/glygen_phosphosites.csv      entrez_gene, glygen_ac, position, residue, kinase
#          $HACK_OUT/inventory/glygen_glycosites.csv        entrez_gene, glygen_ac, position, residue, type, subtype,
#                                                           category, glytoucan_ac, source (one row per site x glycan)
#          $HACK_OUT/inventory/glygen_cache/                small per-protein count files (resume support)
# =====================================================================================================

# Postpone evaluation of type hints (so they are only documentation and never run as code).
from __future__ import annotations

# Built-in modules: csv (tables), io (read text as a file), json (web answers and cache files), os (paths, settings),
# sys (messages to the error stream), time (pauses between retries), urllib.request (web requests).
import csv, io, json, os, sys, time, urllib.request
# ThreadPoolExecutor runs several web requests at the same time.
from concurrent.futures import ThreadPoolExecutor
# Type names used only in the function type hints (no effect on results).
from typing import Any, Dict, List, Optional, Tuple

# Where the pipeline's tables are, and where the inventory goes (outside the repo).
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
# The inventory folder and, inside it, the cache of per-protein results (so a rerun can resume).
INV = os.path.join(OUT, "inventory"); CACHE = os.path.join(INV, "glygen_cache")
# Create both folders if they do not exist yet.
os.makedirs(CACHE, exist_ok=True)
# GlyGen endpoints.
DATA = "https://data.glygen.org/ln2data/releases/data/current/reviewed/"
API = "https://api.glygen.org/protein/detail/{}/"

def get(url: str, tries: int = 4) -> str:
    """Download a URL (with a few retries), returning the text."""
    # Try up to `tries` times.
    for i in range(tries):
        try:
            # download (with a name for our script, and a 2-minute time limit) and return the text
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "motrpac-hackathon-inventory"}), timeout=120) as r:
                return r.read().decode("utf-8")
        # on any failure: give up (raise the error) on the last try, otherwise wait 3 s, 6 s, 9 s ... and retry
        except Exception as e:
            if i == tries - 1: raise
            time.sleep(3 * (i + 1))

def cached(name: str) -> str:
    """A GlyGen data file, downloaded once into the inventory folder and reused (reproducible reruns)."""
    # where the local copy lives
    f = os.path.join(INV, "glygen_files", name)
    # not downloaded yet: create the folder, download and save it
    if not os.path.exists(f):
        os.makedirs(os.path.dirname(f), exist_ok=True)
        open(f, "w").write(get(DATA + name))
    # return the file's contents
    return open(f).read()

def read_csv_text(text: str) -> List[Dict[str, str]]:
    """Parse CSV text into a list of dicts."""
    # each row becomes a dictionary keyed by the column names
    return list(csv.DictReader(io.StringIO(text)))

# ---- 1. our genes -> GlyGen canonical accession --------------------------------------------------------
# Our 471 genes (Entrez ID, symbol, UniProt) from step 2.
genes = read_csv_text(open(os.path.join(OUT, "02_nodes_string.csv")).read())
# GlyGen's Entrez cross-references (canonical accession <-> Entrez) and masterlist (all isoforms per canonical).
xref = read_csv_text(cached("human_protein_xref_geneid.csv"))
master = read_csv_text(cached("human_protein_masterlist.csv"))
# Lookup: Entrez gene ID -> GlyGen canonical accession(s).
by_entrez = {}
for r in xref: by_entrez.setdefault(r["xref_id"], []).append(r["uniprotkb_canonical_ac"])
# Lookup: any UniProt accession (isoform suffix removed) -> its GlyGen canonical accession.
by_acc = {}
for r in master:
    # every reviewed isoform, unreviewed isoform and the canonical accession itself (separated by | or ;)
    for iso in (r["reviewed_isoforms"] + "|" + r["unreviewed_isoforms"] + "|" + r["uniprotkb_canonical_ac"]).replace(";", "|").split("|"):
        if iso: by_acc[iso.split("-")[0]] = r["uniprotkb_canonical_ac"]
# One mapping row per gene.
mapping = []
for g in genes:
    # GlyGen canonical accessions listed for this gene's Entrez ID
    cands = by_entrez.get(g["entrez_gene"], [])
    # prefer the canonical whose accession matches our UniProt; else the first Entrez match; else our UniProt
    acc_ours = (g.get("uniprot") or "").split("-")[0]
    ac = next((c for c in cands if c.split("-")[0] == acc_ours), cands[0] if cands else by_acc.get(acc_ours))
    # record how the accession was found: Entrez and UniProt agree, Entrez only, UniProt only, or not mapped
    how = "entrez+uniprot" if cands and ac.split("-")[0] == acc_ours else ("entrez" if cands else ("uniprot" if ac else "unmapped"))
    # save the gene's IDs, the chosen accession (blank if unmapped) and the mapping route
    mapping.append({"entrez_gene": g["entrez_gene"], "gene_symbol": g["gene_symbol"], "uniprot": g.get("uniprot", ""), "glygen_ac": ac or "", "mapped_via": how})
# Report how many genes got a GlyGen accession.
print(f"mapped {sum(1 for m in mapping if m['glygen_ac'])} of {len(mapping)} genes to GlyGen accessions", file=sys.stderr)

# ---- 2. per-protein counts from the API ------------------------------------------------------------------
def summarise(d: Dict[str, Any]) -> Tuple[Dict[str, int], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Reduce one GlyGen protein record to counts (plus the phosphosite list)."""
    # the counts for this protein, keyed by column name
    c = {}
    # every data section: number of records
    # (lists: count their entries; the GO annotation is a nested dictionary, so count the GO terms inside it)
    for k, v in d.items():
        if isinstance(v, list): c[f"n_{k}"] = len(v)
        elif isinstance(v, dict) and k in ("go_annotation",):
            c["n_go_terms"] = sum(len(cat.get("go_terms", [])) for cat in v.get("categories", []))
    # GlyGen's own statistics per table (glycosylation_reported, _predicted, ...)
    for t in d.get("section_stats", []):
        for s in t.get("table_stats", []):
            # column name = <table>__<field>, value = GlyGen's count
            c[f"{t['table_id']}__{s['field']}"] = s["count"]
    # glycosylation. Records without a position are protein-level evidence (mostly the O-GlcNAc Database:
    # "this protein is O-GlcNAcylated, site not known"); they are counted separately from sites.
    allgly = d.get("glycosylation", [])
    # glycosylation records without a position
    nosite = [g for g in allgly if not g.get("start_pos")]
    # 1 if there is any protein-level (no position) glycosylation evidence, else 0
    c["gly_protein_level_no_site"] = int(len(nosite) > 0)
    # 1 if any of that protein-level evidence is O-GlcNAcylation, else 0
    c["gly_protein_level_oglcnac"] = int(any((g.get("subtype") or "") == "O-GlcNAcylation" for g in nosite))
    # number of distinct positions with O-GlcNAcylation
    c["gly_oglcnac_sites"] = len({g.get("start_pos") for g in allgly if g.get("start_pos") and (g.get("subtype") or "") == "O-GlcNAcylation"})
    # glycosylation records WITH a position, and the distinct (type, position) pairs
    gly = [g for g in allgly if g.get("start_pos")]
    sites = {(g.get("type"), g.get("start_pos")) for g in gly if g.get("start_pos")}
    # number of distinct (type, position) pairs, and distinct positions per linkage type (N, O, C)
    c["gly_sites_unique"] = len(sites)
    for ty in ("N-linked", "O-linked", "C-linked"):
        c[f"gly_sites_{ty.split('-')[0]}"] = len({p for t, p in sites if t == ty})
    # distinct positions per evidence category (e.g. reported, predicted, literature mining)
    cats = {}
    for g in gly:
        # a record can carry several categories (site_category_dict) or just one (site_category)
        for cat in (g.get("site_category_dict") or {g.get("site_category"): True}):
            cats.setdefault(cat, set()).add(g.get("start_pos"))
    # one column per category
    for cat, s in cats.items(): c[f"gly_sites_{cat}"] = len(s)
    # glycan structures observed at a site (excluding the generic O-GlcNAc monosaccharide on protein-level records)
    c["glycans_at_sites"] = len({g["glytoucan_ac"] for g in gly if g.get("glytoucan_ac")})
    # phosphorylation: unique sites and sites with a kinase
    # the phosphorylation records
    ph = d.get("phosphorylation", [])
    c["phospho_sites_unique"] = len({p.get("start_pos") for p in ph})
    c["phospho_sites_with_kinase"] = len({p.get("start_pos") for p in ph if p.get("kinase_uniprot_canonical_ac") or p.get("kinase_gene_name")})
    # one row per phosphorylation record: position, residue, kinase (gene name, else accession) and source databases
    # (PubMed entries are left out of the source list)
    sites_list = [{"position": p.get("start_pos"), "residue": p.get("residue"),
                   "kinase": p.get("kinase_gene_name") or p.get("kinase_uniprot_canonical_ac") or "",
                   "source": ";".join(sorted({e.get("database", "") for e in p.get("evidence", []) if e.get("database") != "PubMed"}))} for p in ph]
    # glycation and other modifications, if present
    c["glycation_sites_unique"] = len({p.get("start_pos") for p in d.get("glycation", [])})
    # glycosylation sites with a position (for site-level matching to phosphosites)
    # one row per positioned glycosylation record: position, residue, type, subtype, category, glycan and sources
    gsites = [{"position": g.get("start_pos"), "residue": g.get("residue") or g.get("start_aa") or "", "type": g.get("type") or "",
               "subtype": g.get("subtype") or "", "category": g.get("site_category") or "", "glytoucan_ac": g.get("glytoucan_ac") or "",
               "source": ";".join(sorted({e.get("database", "") for e in g.get("evidence", []) if e.get("database") != "PubMed"}))} for g in gly]
    # return the counts, the phosphosite rows and the glycosite rows
    return c, sites_list, gsites

def fetch(m: Dict[str, str]) -> Tuple[Dict[str, str], Optional[Dict[str, int]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Fetch (or read from cache) one protein's counts."""
    # the protein's GlyGen accession; unmapped genes get no counts
    ac = m["glygen_ac"]
    if not ac: return m, None, [], []
    # this protein's cache file
    f = os.path.join(CACHE, ac + ".json")
    # a complete cache file (it has the glycosite list) is reused instead of calling the API again
    if os.path.exists(f):
        j = json.load(open(f))
        if "glyco" in j: return m, j["counts"], j["phospho"], j["glyco"]
    # call the API and summarise the answer; a failure is reported and the protein gets no counts
    try:
        c, ps, gs = summarise(json.loads(get(API.format(ac))))
    except Exception as e:
        print(f"  failed {ac}: {e}", file=sys.stderr); return m, None, [], []
    # save the result to the cache and return it
    json.dump({"counts": c, "phospho": ps, "glyco": gs}, open(f, "w"))
    return m, c, ps, gs

# Fetch all proteins, 6 at a time; results come back in the same order as the genes.
with ThreadPoolExecutor(max_workers=6) as ex:
    results = list(ex.map(fetch, mapping))

# ---- 3. write the tables ---------------------------------------------------------------------------------
# Count columns = every count name seen in any fetched protein, alphabetically.
cols = sorted({k for _, c, _, _ in results if c for k in c})
# Per-protein counts table: IDs, mapping route, fetched (1/0), then one column per count (0 if absent).
with open(os.path.join(INV, "glygen_protein_counts.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["entrez_gene", "gene_symbol", "uniprot", "glygen_ac", "mapped_via", "fetched"] + cols)
    for m, c, _, _ in results:
        w.writerow([m["entrez_gene"], m["gene_symbol"], m["uniprot"], m["glygen_ac"], m["mapped_via"], int(c is not None)] + [(c or {}).get(k, 0) for k in cols])
# GlyGen phosphosites table: one row per phosphorylation record.
with open(os.path.join(INV, "glygen_phosphosites.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["entrez_gene", "gene_symbol", "glygen_ac", "position", "residue", "kinase", "source"])
    for m, _, ps, _ in results:
        for p in ps: w.writerow([m["entrez_gene"], m["gene_symbol"], m["glygen_ac"], p["position"], p["residue"], p["kinase"], p["source"]])
# GlyGen glycosites table: one row per positioned glycosylation record.
with open(os.path.join(INV, "glygen_glycosites.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["entrez_gene", "gene_symbol", "glygen_ac", "position", "residue", "type", "subtype", "category", "glytoucan_ac", "source"])
    for m, _, _, gs in results:
        for g in gs: w.writerow([m["entrez_gene"], m["gene_symbol"], m["glygen_ac"], g["position"], g["residue"], g["type"], g["subtype"], g["category"], g["glytoucan_ac"], g["source"]])
# Report how many proteins were fetched.
print(f"fetched {sum(1 for _, c, _, _ in results if c)} proteins -> {INV}", file=sys.stderr)
