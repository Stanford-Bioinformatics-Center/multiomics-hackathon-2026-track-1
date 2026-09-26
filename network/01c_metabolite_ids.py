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
#      is a structure-level match, the most reliable route. PubChem is used as a fallback for any
#      structure UniChem does not cover (ChEBI IDs are listed among PubChem synonyms).
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

# Built-in modules: csv (read/write tables), json (decode web answers), os (file paths, settings),
# re (text patterns), time (pauses), urllib (web requests).
import csv, json, os, re, time, urllib.error, urllib.parse, urllib.request
# ThreadPoolExecutor runs several lookups at the same time (4 at once, to stay polite to the servers).
from concurrent.futures import ThreadPoolExecutor

# Results folder: the HACK_OUT setting if given, otherwise the same default folder the R steps use.
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
# Optional second copy of the output (e.g. on the Desktop for a teammate); empty means no copy.
EXPORT_COPY = os.path.expanduser(os.environ.get("EXPORT_COPY", ""))
# Web address of the RefMet service (metabolite names -> IDs, classes, structures).
REFMET = "https://www.metabolomicsworkbench.org/rest/refmet"
# Web address of UniChem (structure fingerprint -> IDs in other chemical databases, incl. ChEBI).
UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
# Web address of PubChem's synonym list for one compound number (%s is replaced by the number).
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON"
# Web address of EBI's Ontology Lookup Service, searching ChEBI for an exact name (%s = the name).
OLS = "https://www.ebi.ac.uk/ols4/api/search?q=%s&ontology=chebi&exact=true&queryFields=label&rows=5"
# UniChem numbers each database it knows; ChEBI is database number 7.
UNICHEM_CHEBI = 7


# Helper: download a web address and decode its JSON answer, retrying when the server is busy.
def fetch(url, data=None, headers=None):
    """Download a web address and decode its JSON answer. Retry when the server is busy; None on failure."""
    # Try up to 6 times.
    for attempt in range(6):
        # Attempt the download; problems are handled in the "except" branches below.
        try:
            # Build the request (a POST if `data` is given, otherwise a plain GET).
            req = urllib.request.Request(url, data=data, headers=headers or {})
            # Send it, wait up to 60 seconds, and read the answer.
            body = urllib.request.urlopen(req, timeout=60).read()
            # Decode the JSON answer; an empty answer means "nothing found".
            return json.loads(body) if body else {}
        # The server answered with an error code.
        except urllib.error.HTTPError as e:
            # 404 means the service has no entry for this query: return "nothing found".
            if e.code == 404:
                return {}  # nothing found
            # 429/500/502/503 mean "busy" or a temporary fault: wait a little longer each time...
            if e.code in (429, 500, 502, 503):
                time.sleep(2 * (attempt + 1))  # pause before retrying
                # ...and try again.
                continue
            # Any other error code: give up on this query.
            return None
        # A network hiccup (timeout, dropped connection).
        except Exception:
            # Wait a little longer each time, then loop round and try again.
            time.sleep(2 * (attempt + 1))
    # All attempts failed.
    return None


# Helper: RefMet sometimes answers with a list of records instead of one record.
def first_record(d):
    """RefMet sometimes answers with a list; return its first entry as a dictionary (or empty)."""
    # If the answer is a list, take its first entry (or nothing if the list is empty).
    if isinstance(d, list):
        d = d[0] if d else {}  # first entry, or empty if the list is empty
    # Return it if it is a proper record, otherwise an empty record.
    return d if isinstance(d, dict) else {}


# Step 1: look up one metabolite name in RefMet.
def refmet_lookup(name):
    """Step 1: RefMet record (ID, class, PubChem CID, InChIKey) for one metabolite name."""
    # "/" cannot appear inside a web address, so send it as "_"; then make the name web-safe.
    q = urllib.parse.quote(name.replace("/", "_"), safe="")
    # First try RefMet's exact-name lookup.
    rec = first_record(fetch(f"{REFMET}/name/{q}/all"))
    # If it returned a RefMet ID, we are done.
    if rec.get("refmet_id"):
        return rec  # found: return the record
    # Otherwise ask RefMet's name matcher which standard name this corresponds to.
    m = first_record(fetch(f"{REFMET}/match/{q}"))
    # The standard name it suggests ("-" means it found none).
    std = m.get("refmet_name")
    # If it suggested a name...
    if std and std != "-":
        # ...look that standard name up in detail.
        rec = first_record(fetch(f"{REFMET}/name/{urllib.parse.quote(std, safe='')}/all"))
        # If the detail lookup has nothing, keep the basic facts the matcher already gave us.
        if not rec.get("refmet_id"):
            # keep the matcher's basic facts
            rec = {"refmet_id": m.get("refmet_id"), "name": std, "super_class": m.get("super_class"),
                   "main_class": m.get("main_class")}
        # Return what we found.
        return rec
    # RefMet does not know this name at all.
    return {}


# Helper: put ChEBI IDs in numeric order.
def chebi_sort(ids):
    """Sort ChEBI IDs by their number (CHEBI:422 before CHEBI:16651)."""
    # Remove duplicates, then sort by the number after "CHEBI:".
    return sorted(set(ids), key=lambda x: int(x.split(":")[1]))


# Step 2: exact structure -> ChEBI, via UniChem.
def chebi_from_inchikey(inchikey):
    """Step 2: ChEBI IDs for an exact structure, via UniChem."""
    # Ask UniChem which database entries share this structure fingerprint.
    d = fetch(UNICHEM, json.dumps({"type": "inchikey", "compound": inchikey}).encode(),
              {"Content-Type": "application/json"})
    # No answer: no ChEBI IDs.
    if not d:
        return []  # nothing found
    # Keep only the ChEBI entries (database number 7), written as "CHEBI:<number>".
    ids = ["CHEBI:" + str(s["compoundId"]).upper().replace("CHEBI:", "")
           for comp in d.get("compounds", []) for s in comp.get("sources", []) if s.get("id") == UNICHEM_CHEBI]
    # Return them in numeric order.
    return chebi_sort(ids)


# Fallback for step 2: ChEBI IDs listed among a PubChem compound's alternative names.
def chebi_from_pubchem(cid):
    """Fallback for step 2: ChEBI IDs listed among a PubChem compound's synonyms."""
    # Download the compound's list of synonyms (alternative names and IDs).
    d = fetch(PUBCHEM % cid)
    # No answer: no ChEBI IDs.
    if not d:
        return []  # nothing found
    # Flatten the answer into one list of synonyms.
    syn = [s for info in d.get("InformationList", {}).get("Information", []) for s in info.get("Synonym", [])]
    # Keep only synonyms that are exactly a ChEBI ID ("CHEBI:" followed by digits), in numeric order.
    return chebi_sort(s.upper() for s in syn if re.fullmatch(r"CHEBI:\d+", s, re.I))


# Step 3 (only for metabolites without a structure): one fixed reformat, then an exact ChEBI name match.
def chebi_from_exact_label(name):
    """Step 3 (no-structure metabolites only): one fixed reformat, then an exact ChEBI name match."""
    # Split the name at its first space into class and chains, e.g. "PC" and "16:0_18:1".
    m = re.fullmatch(r"(\S+) (.+)", name)
    # Names without a space are not reformatted: no match attempted.
    if not m:
        return []  # name has no space: skip
    # Rebuild in ChEBI's style: "PC(16:0_18:1)".
    label = f"{m.group(1)}({m.group(2)})"
    # Search ChEBI for that name.
    d = fetch(OLS % urllib.parse.quote(label))
    # No answer: no ChEBI IDs.
    if not d:
        return []  # nothing found
    # Accept only entries whose name is EXACTLY the reformatted label, in numeric order.
    return chebi_sort(x["obo_id"] for x in d.get("response", {}).get("docs", []) if x.get("label") == label)


# All steps for one metabolite, producing one row of the output table.
def lookup(name):
    """All steps for one metabolite; returns one output row."""
    # Step 1: RefMet record.
    rec = refmet_lookup(name)
    # Its PubChem compound number and structure fingerprint (blank if RefMet has none).
    cid, ik = rec.get("pubchem_cid") or "", rec.get("inchi_key") or ""
    # Start with "no ChEBI found".
    ids, method = [], "none"
    # If the structure is known...
    if ik:
        # ...match by structure through UniChem.
        ids, method = chebi_from_inchikey(ik), "unichem_inchikey"
        # If UniChem has nothing and there is a PubChem number, try PubChem's synonyms.
        if not ids and cid:
            ids, method = chebi_from_pubchem(cid), "pubchem_synonym"  # PubChem fallback
    # If there is no single structure (mostly lipid species)...
    else:
        # ...try the exact-name match only.
        ids, method = chebi_from_exact_label(name), "chebi_exact_label"
    # If nothing was found by any route, say so.
    if not ids:
        method = "none"  # record that nothing was found
    # Return one output row: the name, RefMet facts, structure IDs and ChEBI result.
    return {"metabolite": name, "refmet_name": rec.get("name", ""), "refmet_id": rec.get("refmet_id", ""),
            "super_class": rec.get("super_class", ""), "main_class": rec.get("main_class", ""),
            "pubchem_cid": cid, "inchi_key": ik,
            "chebi_id": ids[0] if ids else "", "chebi_all": ";".join(ids), "chebi_method": method}


# The main program.
def main():
    # Open step 1b's node table...
    with open(os.path.join(OUT, "01b_metab_nodes_EE.csv"), newline="") as f:
        # ...and read the metabolite name from every row.
        names = [row["metabolite"] for row in csv.DictReader(f)]
    # Safety check: there must be exactly 450 metabolites.
    assert len(names) == 450, "expected 450 metabolites from step 1b"
    # Look the names up 4 at a time; results come back in the same order as the names.
    with ThreadPoolExecutor(4) as pool:
        rows = list(pool.map(lookup, names))  # run all 450 lookups
    # The column names of the output table.
    cols = list(rows[0].keys())
    # Write the table to the results folder, and to the optional second location.
    for path in [os.path.join(OUT, "01c_metabolite_ids.csv")] + ([EXPORT_COPY] if EXPORT_COPY else []):
        # Open the file for writing...
        with open(path, "w", newline="") as f:
            # ...set up a table writer with our columns...
            w = csv.DictWriter(f, fieldnames=cols)
            # ...write the header row...
            w.writeheader()
            # ...and one row per metabolite.
            w.writerows(rows)
        # Report where it was written.
        print("->", path)
    # Count how many metabolites got a ChEBI ID by each route (and how many got none).
    n = {m: sum(r["chebi_method"] == m for r in rows)
         for m in ("unichem_inchikey", "pubchem_synonym", "chebi_exact_label", "none")}
    # Print that summary.
    print(f"ChEBI found for {len(rows) - n['none']} of {len(rows)}: {n}")
    # Also print how many had a PubChem number and how many RefMet did not know at all.
    print(f"with PubChem CID: {sum(bool(r['pubchem_cid']) for r in rows)}; "
          f"no RefMet record: {sum(not r['refmet_id'] for r in rows)}")


# Run the main program when this file is executed as a script (not when imported by another script).
if __name__ == "__main__":
    main()  # start the program
