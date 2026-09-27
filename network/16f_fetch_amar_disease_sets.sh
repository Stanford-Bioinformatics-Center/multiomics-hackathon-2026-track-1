#!/usr/bin/env bash
# =====================================================================================================
# 16f_fetch_amar_disease_sets.sh — STEP 16f: DOWNLOAD THE AMAR ET AL. 2024 DISEASE TABLES (ONCE)
# =====================================================================================================
# PURPOSE: fetch, into $AMAR_EXT, the published disease proteomics tables that step 16d turns into the disease scores
#   used by steps 17i, 18d, 18c, 18bc, 19 and 20. They come from the public repository of Amar et al. 2024 (MoTrPAC
#   rat training mitochondria paper), github.com/MoTrPAC/motrpac-rat-training-mitochondria, folder disease_datasets/, at
#   a PINNED commit (a rerun gets identical bytes), plus STRING v12's human protein list (used to pick the gene symbol
#   STRING knows when a table lists several). Before 2026-09-27 these scores were read from a file built outside this
#   repository (the Venus project, week 6); step 16d now builds the same table here.
#   The tables: Öhman 2021 (T2D muscle), Chae 2018 (T2D muscle, significant only), Coats 2018 (HCM heart, + its UniProt
#   map), Niu 2022 (NASH / cirrhosis liver), Stocks 2022 (ob/ob liver), Havlenova 2021 (rat HF heart), Park 2019 (mouse
#   MI heart), Yuan 2020 (NAFLD liver, significant only), and Ubaida-Mohien 2019 (muscle ageing; read by step 20).
# HOW TO RUN: bash network/16f_fetch_amar_disease_sets.sh      (run_all.sh step 16f; skipped when the files exist)
#   Needs the GitHub CLI (gh) and curl.
# OUTPUTS: $AMAR_EXT/disease_datasets/<file>, $AMAR_EXT/COMMIT_SHA, $AMAR_EXT/9606.protein.info.v12.0.txt.gz
#   (not committed; fingerprinted by step 0)
# =====================================================================================================
set -euo pipefail
X="${AMAR_EXT:-$HOME/Desktop/output/hackathon-2026-track1/external/amar_2024}"
AMAR_COMMIT="299e540ea2d2671df9cddbf883eaac90c1e7cb3a"          # MoTrPAC/motrpac-rat-training-mitochondria, 2023-09-27
FILES=(
  "ohman_2021_iscience.csv" "chae_2018_emm_stat.csv" "coats_2018_circgpm_stat.csv" "uniprot2gene_name_coats_bg.tsv"
  "niu_2022_molsystbiol_set.csv" "stocks_2022_molcellproteomics_stat.csv" "havlenova_2021_scirep.csv"
  "park_2019_celldeathdis.csv" "yuan_2020_jproteomics.csv" "ubaida_mohien_ 2019_elife_stat.csv"
)
mkdir -p "$X/disease_datasets"
for f in "${FILES[@]}"; do
  enc="${f// /%20}"                                             # one file name has a space in the repository
  gh api "repos/MoTrPAC/motrpac-rat-training-mitochondria/contents/disease_datasets/$enc?ref=$AMAR_COMMIT" \
    -H "Accept: application/vnd.github.raw" > "$X/disease_datasets/$f"
done
echo "$AMAR_COMMIT" > "$X/COMMIT_SHA"
curl -sSfL -o "$X/9606.protein.info.v12.0.txt.gz" "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
echo "Amar et al. disease tables + STRING protein info in $X"
