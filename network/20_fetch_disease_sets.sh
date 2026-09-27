#!/usr/bin/env bash
# =====================================================================================================
# 20_fetch_disease_sets.sh — STEP 20f: DOWNLOAD THE PUBLISHED DISEASE TABLES USED BY STEP 20 (ONCE)
# =====================================================================================================
# PURPOSE: fetch the processed supplementary tables of five papers into $DISEASE_EXT so step 20 runs offline:
#   Kjærgaard et al. 2025, Cell 188:4106 (doi 10.1016/j.cell.2025.05.005)      -> kjaergaard_2025_cell/mmc1-3.xlsx
#   Needham et al. 2024, Cell Metab 36:2542 (doi 10.1016/j.cmet.2024.10.020)   -> needham_2024_cellmetab/mmc2-6.xlsx
#   Sun et al. 2023, Nature 622:329 (doi 10.1038/s41586-023-06592-6)          -> sun_2023_nature/MOESM3.xlsx (UK Biobank age)
#   Gadd et al. 2024, Nat Aging 4:1616 (doi 10.1038/s43587-024-00655-7)       -> gadd_2024_nataging/MOESM3.xlsx (UK Biobank incident T2D)
#   Larsen et al. 2023, Sci Adv 9:eadi7548 (doi 10.1126/sciadv.adi7548)        -> larsen_2023_sciadv/ (authors' GitHub
#     repository fpm-cbmr/HIIT_adipose_project at a fixed commit: cleaned protein matrix, sample names, limma script)
#   Supplementary files come from the publishers' public file servers (ars.els-cdn.com, static-content.springer.com); the GitHub files through `gh api`
#   (needs the GitHub CLI) at the pinned commit, so a rerun gets identical bytes. Ubaida-Mohien et al. 2019 is read
#   from the Amar et al. repository copy that step 16f fetches ($UBAIDA; see step 20).
# HOW TO RUN: bash network/20_fetch_disease_sets.sh      (run_all.sh step 20f; skipped when the files exist)
# OUTPUTS: $DISEASE_EXT/<paper>/...  (not committed; fingerprinted by step 0)
# =====================================================================================================
set -euo pipefail
X="${DISEASE_EXT:-$HOME/Desktop/output/hackathon-2026-track1/external/disease}"
LARSEN_COMMIT="c63dd02a3c840c7dba8352c1d3187bf08e786cdb"      # fpm-cbmr/HIIT_adipose_project, the commit used
UA="Mozilla/5.0"                                               # the file server refuses clients without a browser name
mkdir -p "$X/kjaergaard_2025_cell" "$X/needham_2024_cellmetab" "$X/larsen_2023_sciadv/data" "$X/larsen_2023_sciadv/R" "$X/sun_2023_nature" "$X/gadd_2024_nataging"
# Elsevier supplementary files: <pii>-mmc<i>.<ext>
for i in 1 2 3; do curl -sSfL -A "$UA" -o "$X/kjaergaard_2025_cell/mmc$i.xlsx" "https://ars.els-cdn.com/content/image/1-s2.0-S009286742500515X-mmc$i.xlsx"; done
for i in 2 3 4 5 6; do curl -sSfL -A "$UA" -o "$X/needham_2024_cellmetab/mmc$i.xlsx" "https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004169-mmc$i.xlsx"; done
# Springer Nature supplementary files (UK Biobank Olink): <article>_MOESM<i>_ESM.xlsx
curl -sSfL -A "$UA" -o "$X/sun_2023_nature/MOESM3.xlsx" "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx"
curl -sSfL -A "$UA" -o "$X/gadd_2024_nataging/MOESM3.xlsx" "https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-024-00655-7/MediaObjects/43587_2024_655_MOESM3_ESM.xlsx"
# Larsen 2023: authors' repository at the pinned commit
for f in data/Exprs_adipose_clean.txt data/File_names.txt R/Limma_and_GeneSet_analysis.R R/Preprocessing_and_Figure1.R README.md; do
  gh api "repos/fpm-cbmr/HIIT_adipose_project/contents/$f?ref=$LARSEN_COMMIT" -H "Accept: application/vnd.github.raw" > "$X/larsen_2023_sciadv/$f"
done
echo "$LARSEN_COMMIT" > "$X/larsen_2023_sciadv/COMMIT_SHA"
echo "disease tables in $X"
