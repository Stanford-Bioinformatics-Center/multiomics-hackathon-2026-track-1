#!/usr/bin/env python3
# =====================================================================================================
# 01c_metabolite_ids.py — STEP 1c: ChEBI (AND OTHER DATABASE) IDs FOR THE 450 METABOLITE NODES
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Step 1b made 450 metabolite nodes, named with RefMet names (e.g. "Lactic acid", "PC 16:0_18:1").
#   Other databases identify molecules by their own IDs. This script looks up, for every metabolite, its
#   ChEBI ID (Chemical Entities of Biological Interest, the EBI's reference dictionary of small
#   molecules), plus its RefMet ID, chemical class, PubChem CID and InChIKey, so teammates and other
#   tools can connect our metabolites to outside resources.
#
# HOW THE LOOKUP WORKS (three public web services, queried by name or structure; no data values are sent)
#   1. RefMet (Metabolomics Workbench REST API): for each RefMet name, fetch its RefMet ID, class,
#      PubChem CID and InChIKey. The InChIKey is a fixed-length fingerprint of the exact chemical
#      structure. Names containing "/" are sent with "_" instead because "/" breaks the web address; if
#      the exact-name lookup fails, RefMet's own name matcher is tried.
#   2. UniChem (EBI's cross-reference service between chemical databases): InChIKey -> ChEBI ID(s). This
#      is a structure-level match, the most reliable route. If PubChem is reachable it is used as a
#      fallback for any structure UniChem does not cover (ChEBI IDs are listed among PubChem synonyms).
#   3. For metabolites with no single structure (most lipids are defined only as a "species", e.g.
#      "PC 16:0_18:1" = a phosphatidylcholine with a 16:0 and an 18:1 chain in unknown positions), there
#      is no InChIKey. For those only, the name is reformatted in one fixed way to ChEBI's style,
#      "PC 16:0_18:1" -> "PC(16:0_18:1)", and accepted ONLY if a ChEBI entry has exactly that name
#      (via EBI's Ontology Lookup Service). No fuzzy or approximate name matching is done: a blank is
#      better than a wrong ID.
#
# HOW TO READ THE OUTPUT
#   chebi_id      the first ChEBI ID (lowest number) — use this if you need exactly one
#   chebi_all     every ChEBI ID found; several can share one structure (ChEBI sometimes keeps separate
#                 entries for the same molecule, e.g. an acid and its charged form)
#   chebi_method  unichem_inchikey (structure match), pubchem_synonym (structure via PubChem),
#                 chebi_exact_label (name match, lipids only), or none
#   Most lipid species have no ChEBI entry at all; they are left blank on purpose.
#
# TECH STACK
#   Python 3 standard library only (urllib, json, csv, re, concurrent.futures); no packages to install.
#   Needs internet access. Web services change over time, so a rerun can differ slightly; the output
#   records the method used for every ID.
#
# INPUT:   $HACK_OUT/01b_metab_nodes_EE.csv  (the 450 metabolite names, from step 1b)
# OUTPUT:  $HACK_OUT/01c_metabolite_ids.csv   (one row per metabolite)
#          optional second copy: path in $EXPORT_COPY (e.g. ~/Desktop/metabolite_chebi_ids.csv)
# =====================================================================================================

# Standard-library modules: file paths, CSV, JSON, regular expressions, time, web requests, threads.
import csv, json, os, re, time, urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

# Results folder (same default as the R steps; override with HACK_OUT).
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
# Optional second copy of the output (e.g. on the Desktop for a teammate); empty = no copy.
EXPORT_COPY = os.path.expanduser(os.environ.get("EXPORT_COPY", ""))
# Base addresses of the three web services.
REFMET = "https://www.metabolomicsworkbench.org/rest/refmet"
UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON"
OLS = "https://www.ebi.ac.uk/ols4/api/search?q=%s&ontology=chebi&exact=true&queryFields=label&rows=5"
# UniChem's internal number for the ChEBI database.
UNICHEM_CHEBI = 7


def fetch(url, data=None, headers=None):
    """Download a web address and decode its JSON answer. Retry when the server is busy; None on failure."""
    # up to 6 attempts, waiting a little longer each time the server says it is busy
    for attempt in range(6):
        try:
            # send the request (a POST if `data` is given, otherwise a GET) and read the answer
            req = urllib.request.Request(url, data=data, headers=headers or {})
            body = urllib.request.urlopen(req, timeout=60).read()
            # an empty answer means "nothing found"
            return json.loads(body) if body else {}
        except urllib.error.HTTPError as e:
            # 404 = the service has no entry for this query
            if e.code == 404:
                return {}
            # busy / temporary server errors: wait and try again
            if e.code in (429, 500, 502, 503):
                time.sleep(2 * (attempt + 1))
                continue
            return None
        except Exception:
            # network hiccup: wait and try again
            time.sleep(2 * (attempt + 1))
    return None


def first_record(d):
    """RefMet sometimes answers with a list; return its first entry as a dictionary (or empty)."""
    if isinstance(d, list):
        d = d[0] if d else {}
    return d if isinstance(d, dict) else {}


def refmet_lookup(name):
    """Step 1: RefMet record (ID, class, PubChem CID, InChIKey) for one metabolite name."""
    # "/" cannot appear inside a web address, so it is sent as "_" (the chemical class is unaffected)
    q = urllib.parse.quote(name.replace("/", "_"), safe="")
    # first try the exact-name lookup
    rec = first_record(fetch(f"{REFMET}/name/{q}/all"))
    if rec.get("refmet_id"):
        return rec
    # otherwise ask RefMet's name matcher for the standard name, then look that up
    m = first_record(fetch(f"{REFMET}/match/{q}"))
    std = m.get("refmet_name")
    if std and std != "-":
        rec = first_record(fetch(f"{REFMET}/name/{urllib.parse.quote(std, safe='')}/all"))
        if not rec.get("refmet_id"):
            # the matcher knows the name but the detail lookup does not: keep what the matcher gave
            rec = {"refmet_id": m.get("refmet_id"), "name": std, "super_class": m.get("super_class"),
                   "main_class": m.get("main_class")}
        return rec
    return {}


def chebi_sort(ids):
    """Sort ChEBI IDs by their number (CHEBI:422 before CHEBI:16651)."""
    return sorted(set(ids), key=lambda x: int(x.split(":")[1]))


def chebi_from_inchikey(inchikey):
    """Step 2: ChEBI IDs for an exact structure, via UniChem."""
    d = fetch(UNICHEM, json.dumps({"type": "inchikey", "compound": inchikey}).encode(),
              {"Content-Type": "application/json"})
    if not d:
        return []
    # collect every ChEBI entry UniChem lists for this structure
    ids = ["CHEBI:" + str(s["compoundId"]).upper().replace("CHEBI:", "")
           for comp in d.get("compounds", []) for s in comp.get("sources", []) if s.get("id") == UNICHEM_CHEBI]
    return chebi_sort(ids)


def chebi_from_pubchem(cid):
    """Fallback for step 2: ChEBI IDs listed among a PubChem compound's synonyms."""
    d = fetch(PUBCHEM % cid)
    if not d:
        return []
    syn = [s for info in d.get("InformationList", {}).get("Information", []) for s in info.get("Synonym", [])]
    return chebi_sort(s.upper() for s in syn if re.fullmatch(r"CHEBI:\d+", s, re.I))


def chebi_from_exact_label(name):
    """Step 3 (no-structure metabolites only): one fixed reformat, then an exact ChEBI name match."""
    # "PC 16:0_18:1" -> "PC(16:0_18:1)"; names without a space are not reformatted
    m = re.fullmatch(r"(\S+) (.+)", name)
    if not m:
        return []
    label = f"{m.group(1)}({m.group(2)})"
    d = fetch(OLS % urllib.parse.quote(label))
    if not d:
        return []
    # accept only entries whose name is exactly the reformatted label
    return chebi_sort(x["obo_id"] for x in d.get("response", {}).get("docs", []) if x.get("label") == label)


def lookup(name):
    """All steps for one metabolite; returns one output row."""
    rec = refmet_lookup(name)
    cid, ik = rec.get("pubchem_cid") or "", rec.get("inchi_key") or ""
    ids, method = [], "none"
    if ik:
        # structure known: structure-level match
        ids, method = chebi_from_inchikey(ik), "unichem_inchikey"
        if not ids and cid:
            ids, method = chebi_from_pubchem(cid), "pubchem_synonym"
    else:
        # no single structure (mostly lipid species): exact-name match only
        ids, method = chebi_from_exact_label(name), "chebi_exact_label"
    if not ids:
        method = "none"
    return {"metabolite": name, "refmet_name": rec.get("name", ""), "refmet_id": rec.get("refmet_id", ""),
            "super_class": rec.get("super_class", ""), "main_class": rec.get("main_class", ""),
            "pubchem_cid": cid, "inchi_key": ik,
            "chebi_id": ids[0] if ids else "", "chebi_all": ";".join(ids), "chebi_method": method}


def main():
    # read the 450 metabolite names from step 1b's node table (first column)
    with open(os.path.join(OUT, "01b_metab_nodes_EE.csv"), newline="") as f:
        names = [row["metabolite"] for row in csv.DictReader(f)]
    assert len(names) == 450, "expected 450 metabolites from step 1b"
    # look them up 4 at a time (polite to the servers), keeping the input order
    with ThreadPoolExecutor(4) as pool:
        rows = list(pool.map(lookup, names))
    # write the table, and the optional second copy
    cols = list(rows[0].keys())
    for path in [os.path.join(OUT, "01c_metabolite_ids.csv")] + ([EXPORT_COPY] if EXPORT_COPY else []):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print("->", path)
    # summary: how many got a ChEBI ID, by method
    n = {m: sum(r["chebi_method"] == m for r in rows)
         for m in ("unichem_inchikey", "pubchem_synonym", "chebi_exact_label", "none")}
    print(f"ChEBI found for {len(rows) - n['none']} of {len(rows)}: {n}")
    print(f"with PubChem CID: {sum(bool(r['pubchem_cid']) for r in rows)}; "
          f"no RefMet record: {sum(not r['refmet_id'] for r in rows)}")


# run main() when the file is executed as a script
if __name__ == "__main__":
    main()
