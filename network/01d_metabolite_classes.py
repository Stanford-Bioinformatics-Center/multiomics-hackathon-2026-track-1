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
#   A metabolite with no RefMet class is counted as "Unclassified" (never dropped), so the counts always
#   add up to 450.
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

# Standard-library modules: CSV files, counting, file paths.
import csv, os
from collections import defaultdict

# Results folder (override with HACK_OUT) and optional second copy.
OUT = os.environ.get("HACK_OUT", os.path.expanduser("~/Desktop/output/hackathon-2026-track1/network"))
EXPORT_COPY = os.path.expanduser(os.environ.get("EXPORT_COPY", ""))
# The super classes counted as lipids (LIPID MAPS categories).
LIPID = {"Fatty Acyls", "Glycerolipids", "Glycerophospholipids", "Sphingolipids", "Sterol Lipids",
         "Prenol Lipids", "Saccharolipids", "Polyketides"}


def clean(x):
    """A blank or '-' class becomes 'Unclassified' so nothing is silently dropped."""
    x = (x or "").strip()
    return x if x and x != "-" else "Unclassified"


def main():
    # read step 1c's table: one row per metabolite
    with open(os.path.join(OUT, "01c_metabolite_ids.csv"), newline="") as f:
        rows = list(csv.DictReader(f))
    total = len(rows)
    assert total == 450, "expected 450 metabolites from step 1c"

    # group the metabolites: once by super class, once by (super class, main class)
    groups = defaultdict(list)
    for r in rows:
        sup, main_ = clean(r["super_class"]), clean(r["main_class"])
        groups[("super_class", sup, "")].append(r)
        groups[("main_class", sup, main_)].append(r)

    # one output row per group
    out = []
    for (level, sup, main_), members in groups.items():
        out.append({
            "level": level,
            "is_lipid": "yes" if sup in LIPID else "no",
            "super_class": sup,
            "main_class": main_,
            "n_metabolites": len(members),
            "pct_of_450": round(100 * len(members) / total, 1),
            "n_with_chebi": sum(bool(m["chebi_id"]) for m in members),
            "examples": "; ".join(sorted(m["metabolite"] for m in members)[:5]),
        })
    # order: super classes first (largest first), then main classes grouped under their super class
    super_size = {o["super_class"]: o["n_metabolites"] for o in out if o["level"] == "super_class"}
    out.sort(key=lambda o: (o["level"] != "super_class", -super_size[o["super_class"]], o["super_class"],
                            -o["n_metabolites"], o["main_class"]))

    # safety check: each level adds up to 450
    for lvl in ("super_class", "main_class"):
        assert sum(o["n_metabolites"] for o in out if o["level"] == lvl) == total

    # write the table, and the optional second copy
    for path in [os.path.join(OUT, "01d_metabolite_class_counts.csv")] + ([EXPORT_COPY] if EXPORT_COPY else []):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print("->", path)

    # print the super-class summary and the lipid total
    for o in out:
        if o["level"] == "super_class":
            print(f"{o['n_metabolites']:4d}  {o['pct_of_450']:5.1f}%  {o['super_class']}  (lipid: {o['is_lipid']})")
    n_lipid = sum(o["n_metabolites"] for o in out if o["level"] == "super_class" and o["is_lipid"] == "yes")
    print(f"lipids {n_lipid} / {total}; main classes: {sum(o['level'] == 'main_class' for o in out)}")


# run main() when the file is executed as a script
if __name__ == "__main__":
    main()
