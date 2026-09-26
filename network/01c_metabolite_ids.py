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
# UPSTREAM (what the names already are when they reach this script)
#   The consortium standardised every metabolite name to RefMet during its QC (see step 1b), so our names
#   should be RefMet names already. This script adds identifiers only; it does not change the data.
#
# HOW THE LOOKUP WORKS (public web services, queried by name or structure; no data values are sent)
#   1. RefMet (Metabolomics Workbench REST API): find the metabolite's RefMet record, then fetch its
#      details (RefMet ID, class, PubChem CID, InChIKey) BY RefMet ID. The InChIKey is a fixed-length
#      fingerprint of the exact chemical structure.
#      Why names with "/" need care: in lipid names "/" and "_" mean different things ("/" = the position
#      of each fatty-acid chain is known, "_" = it is not). RefMet's server refuses web addresses that
#      contain an encoded "/", so such names cannot be looked up directly. For those, the name is sent
#      to RefMet's name MATCHER with "_" in place of "/", and the matcher returns the proper record
#      (e.g. "TG 18:1/18:1/18:1", RM0134295, not the less specific "TG 18:1_18:1_18:1"). Details are then
#      fetched by RefMet ID, so no "/" ever has to go into a web address.
#      Every row records whether RefMet's name equals ours exactly (column name_match), so any record
#      that is not exactly our molecule is visible.
#   2. UniChem (EBI's cross-reference service between chemical databases): InChIKey -> ChEBI ID(s). This
#      is a structure-level match, the most reliable route. PubChem is used as a fallback for any
#      structure UniChem does not cover (ChEBI IDs are listed among PubChem synonyms).
#   3. Metabolites with no single structure (most lipids are defined only as a "species", e.g.
#      "PC 16:0_18:1" = a phosphatidylcholine with a 16:0 and an 18:1 chain in unknown positions) have
#      no InChIKey. For any metabolite without an InChIKey (in practice almost all are lipids), the name is
#      reformatted in one fixed way to ChEBI's style, "PC 16:0_18:1" -> "PC(16:0_18:1)", and accepted ONLY
#      if a ChEBI entry has exactly that name (via EBI's Ontology Lookup Service). No fuzzy or approximate
#      name matching is done: a blank is better than a wrong ID.
#   Failed web requests (timeouts, outages) are NOT silently treated as "not found": they are counted,
#   marked lookup_status = "error" in the row, and reported at the end, so a rerun can fill them in.
#
# HOW TO READ THE OUTPUT
#   refmet_name   RefMet's name for the record used; name_match = exact / differs / no_record
#   chebi_id      the first ChEBI ID (lowest number) — use this if you need exactly one
#   chebi_all     every ChEBI ID found; several can share one structure (ChEBI sometimes keeps separate
#                 entries for the same molecule, e.g. an acid and its charged form)
#   chebi_method  unichem_inchikey (structure match), pubchem_synonym (structure via PubChem),
#                 chebi_exact_label (exact name match, metabolites without a structure), or none
#   lookup_status ok, or error if any web request for this row failed
#   Stereochemistry: RefMet assigns a specific stereoisomer where it can (e.g. lactic acid -> L-lactic
#   acid), even when the assay does not separate stereoisomers; the ChEBI ID inherits that choice.
#
# TECH STACK
#   Python 3.9+ standard library only (urllib, json, csv, re, concurrent.futures); nothing to install.
#   Needs internet access. Web services change over time, so a rerun can differ slightly; the output
#   records the method used for every ID. Last run: 2026-09-26.
#
# INPUT:   $HACK_OUT/01b_metab_nodes_EE.csv  (the 450 metabolite names, from step 1b)
# OUTPUT:  $HACK_OUT/01c_metabolite_ids.csv   (one row per metabolite)
#          optional second copy: path in $EXPORT_COPY (e.g. ~/Desktop/metabolite_chebi_ids.csv)
# =====================================================================================================

# Built-in modules: csv (read/write tables), json (decode web answers), os (file paths, settings),
# re (text patterns), threading (a lock for the shared error counter), time (pauses), urllib (web requests).
import csv, json, os, re, threading, time, urllib.error, urllib.parse, urllib.request
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
# Number of attempts per web request before giving up.
ATTEMPTS = 6


# Marker returned by fetch() when a request FAILED (as opposed to "answered: nothing found").
class FetchError(Exception):
    """Raised when a web request fails after all retries, so the failure is never mistaken for 'not found'."""


# Helper: download a web address and decode its JSON answer, retrying when the server is busy.
def fetch(url, data=None, headers=None):
    """Return the decoded JSON answer ({} if the service says 'not found'); raise FetchError on failure."""
    # Try up to ATTEMPTS times.
    for attempt in range(ATTEMPTS):
        # Wait before every retry (not before the first try): 2 s, 4 s, 6 s, ...
        if attempt:
            time.sleep(2 * attempt)  # pause before this retry
        # Attempt the download; problems are handled in the "except" branches below.
        try:
            # Build the request (a POST if `data` is given, otherwise a plain GET).
            req = urllib.request.Request(url, data=data, headers=headers or {})
            # Send it, wait up to 60 seconds, and read the answer.
            body = urllib.request.urlopen(req, timeout=60).read()
            # An empty answer means "nothing found".
            if not body:
                return {}  # empty answer: nothing found
            # Decode the JSON answer (a non-JSON answer raises ValueError, handled below).
            return json.loads(body)
        # The server answered with an error code.
        except urllib.error.HTTPError as e:
            # 404 means the service has no entry for this query: a genuine "nothing found".
            if e.code == 404:
                return {}  # 404: nothing found
            # Busy or temporary faults: loop round and retry after the pause.
            if e.code in (429, 500, 502, 503, 504):
                continue  # busy: retry
            # Any other error code is a real failure.
            raise FetchError(f"HTTP {e.code}: {url}")
        # The answer was not valid JSON (e.g. an HTML error page): a real failure, no point retrying.
        except ValueError:
            raise FetchError(f"non-JSON answer: {url}")  # invalid answer: fail
        # A network hiccup (timeout, dropped connection): loop round and retry after the pause.
        except Exception:
            continue  # hiccup: retry
    # All attempts failed.
    raise FetchError(f"gave up after {ATTEMPTS} attempts: {url}")


# Shared counter of failed requests (several lookups run at once, so a lock protects it).
ERRORS = {"n": 0}
# The lock that protects the counter.
ERR_LOCK = threading.Lock()


# Helper: like fetch(), but a failure is counted and noted in the row instead of stopping the whole run.
def safe_fetch(status, url, data=None, headers=None):
    """fetch() that records a failure in `status` (a one-item list) and returns {} so the row can continue."""
    # Try the request.
    try:
        return fetch(url, data, headers)  # normal case: the answer
    # On failure: count it, mark this row as having an error, and carry on with an empty answer.
    except FetchError:
        with ERR_LOCK:  # take the lock (other lookups run at the same time)
            ERRORS["n"] += 1  # count the failure
        status[0] = "error"  # mark this row
        return {}  # carry on with an empty answer


# Helper: RefMet sometimes answers with a list of records instead of one record.
def first_record(d):
    """RefMet sometimes answers with a list; return its first entry as a dictionary (or empty)."""
    # If the answer is a list, take its first entry (or nothing if the list is empty).
    if isinstance(d, list):
        d = d[0] if d else {}  # first entry, or empty if the list is empty
    # Return it if it is a proper record, otherwise an empty record.
    return d if isinstance(d, dict) else {}


# Step 1: find one metabolite's RefMet record.
def refmet_lookup(name, status):
    """Step 1: RefMet record (ID, class, PubChem CID, InChIKey) for one metabolite name."""
    # Names without "/" can be looked up directly by exact name.
    if "/" not in name:
        # Exact-name lookup (the name made web-safe).
        rec = first_record(safe_fetch(status, f"{REFMET}/name/{urllib.parse.quote(name, safe='')}/all"))
        # If it returned a RefMet ID, we are done.
        if rec.get("refmet_id"):
            return rec  # found by exact name
    # Otherwise (a "/" in the name, or no exact hit): ask RefMet's name matcher. "/" is sent as "_"
    # because RefMet refuses an encoded "/"; the matcher returns the record with the proper notation.
    m = first_record(safe_fetch(status, f"{REFMET}/match/{urllib.parse.quote(name.replace('/', '_'), safe='')}"))
    # The RefMet ID the matcher chose ("-" or missing means it found none).
    rid = m.get("refmet_id")
    # No match: RefMet does not know this name.
    if not rid or rid == "-":
        return {}  # RefMet does not know this name
    # Fetch the full record BY RefMet ID (no "/" in the web address).
    rec = first_record(safe_fetch(status, f"{REFMET}/refmet_id/{urllib.parse.quote(rid, safe='')}/all"))
    # If the detail lookup has nothing, keep the basic facts the matcher already gave us.
    if not rec.get("refmet_id"):
        # keep the matcher's basic facts
        rec = {"refmet_id": rid, "name": m.get("refmet_name"), "super_class": m.get("super_class"),
               "main_class": m.get("main_class")}
    # Return the record.
    return rec


# Helper: put ChEBI IDs in numeric order.
def chebi_sort(ids):
    """Sort ChEBI IDs by their number (CHEBI:422 before CHEBI:16651)."""
    # Remove duplicates, then sort by the number after "CHEBI:".
    return sorted(set(ids), key=lambda x: int(x.split(":")[1]))


# Step 2: exact structure -> ChEBI, via UniChem.
def chebi_from_inchikey(inchikey, status):
    """Step 2: ChEBI IDs for an exact structure, via UniChem."""
    # Ask UniChem which database entries share this structure fingerprint.
    d = safe_fetch(status, UNICHEM, json.dumps({"type": "inchikey", "compound": inchikey}).encode(),
                   {"Content-Type": "application/json"})
    # Keep only the ChEBI entries (database number 7), written as "CHEBI:<number>".
    ids = ["CHEBI:" + str(s["compoundId"]).upper().replace("CHEBI:", "")
           for comp in d.get("compounds", []) for s in comp.get("sources", []) if s.get("id") == UNICHEM_CHEBI]
    # Return them in numeric order (empty if none).
    return chebi_sort(ids)


# Fallback for step 2: ChEBI IDs listed among a PubChem compound's alternative names.
def chebi_from_pubchem(cid, status):
    """Fallback for step 2: ChEBI IDs listed among a PubChem compound's synonyms."""
    # Download the compound's list of synonyms (alternative names and IDs).
    d = safe_fetch(status, PUBCHEM % cid)
    # Flatten the answer into one list of synonyms.
    syn = [s for info in d.get("InformationList", {}).get("Information", []) for s in info.get("Synonym", [])]
    # Keep only synonyms that are exactly a ChEBI ID ("CHEBI:" followed by digits), in numeric order.
    return chebi_sort(s.upper() for s in syn if re.fullmatch(r"CHEBI:\d+", s, re.I))


# Step 3 (only for metabolites without a structure): one fixed reformat, then an exact ChEBI name match.
def chebi_from_exact_label(name, status):
    """Step 3 (metabolites without a structure): one fixed reformat, then an exact ChEBI name match."""
    # Split the name at its first space into class and chains, e.g. "PC" and "16:0_18:1".
    m = re.fullmatch(r"(\S+) (.+)", name)
    # Names without a space are not reformatted: no match attempted.
    if not m:
        return []  # name has no space: skip
    # Rebuild in ChEBI's style: "PC(16:0_18:1)".
    label = f"{m.group(1)}({m.group(2)})"
    # Search ChEBI for that name.
    d = safe_fetch(status, OLS % urllib.parse.quote(label))
    # Accept only entries whose name is EXACTLY the reformatted label, in numeric order.
    return chebi_sort(x["obo_id"] for x in d.get("response", {}).get("docs", []) if x.get("label") == label)


# All steps for one metabolite, producing one row of the output table.
def lookup(name):
    """All steps for one metabolite; returns one output row."""
    # This row's status; any failed request below switches it to "error".
    status = ["ok"]
    # Step 1: RefMet record.
    rec = refmet_lookup(name, status)
    # RefMet's name for that record, and whether it is exactly our name.
    rname = rec.get("name") or ""
    name_match = "no_record" if not rec.get("refmet_id") else ("exact" if rname == name else "differs")  # compare RefMet's name with ours
    # Its PubChem compound number and structure fingerprint (blank if RefMet has none).
    cid, ik = rec.get("pubchem_cid") or "", rec.get("inchi_key") or ""
    # Start with "no ChEBI found".
    ids, method = [], "none"
    # If the structure is known...
    if ik:
        # ...match by structure through UniChem.
        ids, method = chebi_from_inchikey(ik, status), "unichem_inchikey"
        # If UniChem has nothing and there is a PubChem number, try PubChem's synonyms.
        if not ids and cid:
            ids, method = chebi_from_pubchem(cid, status), "pubchem_synonym"  # PubChem fallback
    # If there is no single structure (mostly lipid species)...
    else:
        # ...try the exact-name match only.
        ids, method = chebi_from_exact_label(name, status), "chebi_exact_label"
    # If nothing was found by any route, say so.
    if not ids:
        method = "none"  # record that nothing was found
    # Return one output row: the name, RefMet facts, structure IDs, ChEBI result and status.
    return {"metabolite": name, "refmet_name": rname, "name_match": name_match,
            "refmet_id": rec.get("refmet_id", ""),
            "super_class": rec.get("super_class", ""), "main_class": rec.get("main_class", ""),
            "pubchem_cid": cid, "inchi_key": ik,
            "chebi_id": ids[0] if ids else "", "chebi_all": ";".join(ids), "chebi_method": method,
            "lookup_status": status[0]}


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
    # Print how many had a PubChem number, how RefMet's names compare with ours, and any failures.
    print(f"with PubChem CID: {sum(bool(r['pubchem_cid']) for r in rows)}; name_match: "
          f"{ {k: sum(r['name_match'] == k for r in rows) for k in ('exact', 'differs', 'no_record')} }")
    # List the rows whose RefMet name differs from ours, so they can be checked by eye.
    for r in rows:
        # Only the rows that differ.
        if r["name_match"] == "differs":
            # Our name -> RefMet's name.
            print(f"   name differs: {r['metabolite']!r} -> {r['refmet_name']!r}")
    # Report failed requests; a nonzero count means some blanks may be failures, so rerun.
    print(f"failed web requests: {ERRORS['n']} (rows with lookup_status=error: "
          f"{sum(r['lookup_status'] == 'error' for r in rows)})")


# Run the main program when this file is executed as a script (not when imported by another script).
if __name__ == "__main__":
    main()  # start the program
