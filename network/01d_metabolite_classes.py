#!/usr/bin/env python3
# =====================================================================================================
# 01d_metabolite_classes.py — STEP 1d: HOW MANY OF THE 450 METABOLITES FALL IN EACH CHEMICAL CLASS
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Groups the 450 metabolite nodes by chemical class and counts them. The classes come from RefMet
#   (the Metabolomics Workbench reference nomenclature), which files every metabolite in a hierarchy:
#     super class  (broad family, e.g. "Glycerophospholipids", "Organic acids")
#       main class (narrower group inside it, e.g. "Glycerophosphocholines" = PC lipids, "Amino acids")
#   The class of each metabolite was already looked up in step 1c (01c_metabolite_ids.py); this script
#   only counts, so it needs no internet.
#
# LIPID vs NON-LIPID
#   Following the LIPID MAPS categories used by RefMet, these super classes are counted as lipids:
#   Fatty Acyls, Glycerolipids, Glycerophospholipids, Sphingolipids, Sterol Lipids, Prenol Lipids,
#   Saccharolipids, Polyketides. Everything else is non-lipid (mostly polar: amino acids, organic acids,
#   nucleotides, sugars, ...).
#
# MISSING CLASSES
#   A metabolite with no RefMet class is counted as "Unclassified" (never dropped, never assigned by
#   hand), so the counts always add up to 450.
#
# TECH STACK
#   Python 3 standard library only (csv, collections, os).
#
# INPUT:   $HACK_OUT/01c_metabolite_ids.csv (from step 1c)
# OUTPUT:  $HACK_OUT/01d_metabolite_class_counts.csv, one row per class, with columns:
#            level        "super_class" (broad) or "main_class" (narrower, nested inside a super class)
#            is_lipid     yes / no
#            super_class, main_class, n_metabolites, pct_of_450
#            n_with_chebi how many in that class got a ChEBI ID in step 1c
#            examples     up to 5 member names
#          optional second copy at $EXPORT_COPY (e.g. ~/Desktop/metabolite_class_counts.csv)
# =====================================================================================================

# Built-in modules: csv (read/write tables) and os (file paths, settings).
import csv, os
# defaultdict: a dictionary that starts every new key with an empty list (handy for grouping).
from collections import defaultdict

# Results folder: the HACK_OUT setting if given, otherwise the pipeline's default folder.
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
# Optional second copy of the output (e.g. on the Desktop); empty means no copy.
EXPORT_COPY = os.path.expanduser(os.environ.get("EXPORT_COPY", ""))
# The super classes counted as lipids (the LIPID MAPS categories).
LIPID = {"Fatty Acyls", "Glycerolipids", "Glycerophospholipids", "Sphingolipids", "Sterol Lipids",
         "Prenol Lipids", "Saccharolipids", "Polyketides"}


# Helper: tidy one class name.
def clean(x):
    """A blank or '-' class becomes 'Unclassified' so nothing is silently dropped."""
    # Remove surrounding spaces (and treat a missing value as empty text).
    x = (x or "").strip()
    # Keep a real class name; turn blank or "-" into "Unclassified".
    return x if x and x != "-" else "Unclassified"


# The main program.
def main():
    # Open step 1c's table (one row per metabolite)...
    with open(os.path.join(OUT, "01c_metabolite_ids.csv"), newline="") as f:
        # ...and read all its rows.
        rows = list(csv.DictReader(f))
    # The number of metabolites.
    total = len(rows)
    # Safety check: there must be exactly 450.
    assert total == 450, "expected 450 metabolites from step 1c"

    # Groups of metabolites, keyed by (level, super class, main class).
    groups = defaultdict(list)
    # Go through the metabolites one by one.
    for r in rows:
        # This metabolite's tidied super class and main class.
        sup, main_ = clean(r["super_class"]), clean(r["main_class"])
        # Add it to its super-class group...
        groups[("super_class", sup, "")].append(r)
        # ...and to its main-class group (nested under the super class).
        groups[("main_class", sup, main_)].append(r)

    # The output rows.
    out = []
    # One output row per group.
    for (level, sup, main_), members in groups.items():
        # Describe the group.
        out.append({
            # "super_class" or "main_class"
            "level": level,
            # whether the super class is a lipid category
            "is_lipid": "yes" if sup in LIPID else "no",
            # the class names (main class is blank on super-class rows)
            "super_class": sup,
            "main_class": main_,
            # how many metabolites, and what percentage of the 450
            "n_metabolites": len(members),
            "pct_of_450": round(100 * len(members) / total, 1),
            # how many of them got a ChEBI ID in step 1c
            "n_with_chebi": sum(bool(m["chebi_id"]) for m in members),
            # up to five member names, alphabetically
            "examples": "; ".join(sorted(m["metabolite"] for m in members)[:5]),
        })
    # The size of each super class (used to order the table).
    super_size = {o["super_class"]: o["n_metabolites"] for o in out if o["level"] == "super_class"}
    # Order: all super-class rows first (largest first), then main-class rows grouped under their super
    # class in the same order (largest main class first within each).
    out.sort(key=lambda o: (o["level"] != "super_class", -super_size[o["super_class"]], o["super_class"],
                            -o["n_metabolites"], o["main_class"]))

    # Safety check: at each level the counts must add up to 450.
    for lvl in ("super_class", "main_class"):
        # Sum the counts for this level and compare with the total.
        assert sum(o["n_metabolites"] for o in out if o["level"] == lvl) == total

    # Write the table to the results folder, and to the optional second location.
    for path in [os.path.join(OUT, "01d_metabolite_class_counts.csv")] + ([EXPORT_COPY] if EXPORT_COPY else []):
        # Open the file for writing...
        with open(path, "w", newline="") as f:
            # ...set up a table writer with the columns of the output rows...
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            # ...write the header row...
            w.writeheader()
            # ...and all class rows.
            w.writerows(out)
        # Report where it was written.
        print("->", path)

    # Print the super-class summary on screen.
    for o in out:
        # Only the super-class rows.
        if o["level"] == "super_class":
            # Count, percentage, class name and whether it is a lipid class.
            print(f"{o['n_metabolites']:4d}  {o['pct_of_450']:5.1f}%  {o['super_class']}  (lipid: {o['is_lipid']})")
    # The total number of lipids (sum of the lipid super classes).
    n_lipid = sum(o["n_metabolites"] for o in out if o["level"] == "super_class" and o["is_lipid"] == "yes")
    # Print the lipid total and the number of main classes.
    print(f"lipids {n_lipid} / {total}; main classes: {sum(o['level'] == 'main_class' for o in out)}")


# Run the main program when this file is executed as a script (not when imported by another script).
if __name__ == "__main__":
    main()  # start the program
