#!/usr/bin/env bash
# =====================================================================================================
# network/run_all.sh — RUN THE WHOLE NETWORK PIPELINE, IN ORDER, REPRODUCIBLY
# =====================================================================================================
# PURPOSE: one command from the MoTrPAC results (+ the team's mnet resource) to every table, figure and
#   interactive page, followed by the validation checks and the reproducibility manifest.
# HOW TO RUN (from the repository root):
#   bash network/run_all.sh                 # everything (about 15 minutes on a laptop)
#   bash network/run_all.sh 14 17i          # only steps 14 ... 17i (step labels as listed below)
#   REFRESH_ONLINE=1 bash network/run_all.sh   # also redo the online lookups (RefMet / ChEBI / PubChem, GlyGen)
# ENVIRONMENT (defaults in brackets; see network/README.md, section 4):
#   HACK_OUT [~/Desktop/output/hackathon-2026-track1/network]  tables       HACK_FIG [~/Desktop/output/hackathon]  figures
#   HACK_EXT [~/Desktop/output/hackathon-2026-track1/external/rhea]  Rhea cache
#   MNET_DIR [~/Desktop/output/hackathon/resources/mo_annotation]   the team's mnet resource (edges, PTM)
#   EDGE_SOURCE [mnet]  ("legacy" = the first curated STRING file + direct Rhea)
#   DISEASE_SCORES [~/Desktop/output/week_6/_shared/disease_scores.csv.gz]  Amar et al. 2024 disease sets (steps 17i, 18d, 19)
# ONLINE STEPS AND CACHES: step 1c queries RefMet / UniChem / PubChem (web services change over time), so it
#   is skipped when its output exists unless REFRESH_ONLINE=1; step 5 downloads Rhea once into HACK_EXT; the
#   inventory queries GlyGen once per protein and caches the answers (and GlyGen files) under
#   $HACK_OUT/inventory. With the caches in place a rerun needs no internet and reproduces every output.
# OUTPUTS: as each step documents; logs in $HACK_OUT/logs/<step>.log; network/ENVIRONMENT.md;
#   $HACK_OUT/99_manifest*.csv (the reproducibility check).
# =====================================================================================================
set -euo pipefail

# Folders and defaults (exported so every step sees the same values).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export HACK_OUT="${HACK_OUT:-$HOME/Desktop/output/hackathon-2026-track1/network}"
export HACK_FIG="${HACK_FIG:-$HOME/Desktop/output/hackathon}"
export HACK_EXT="${HACK_EXT:-$HOME/Desktop/output/hackathon-2026-track1/external/rhea}"
export MNET_DIR="${MNET_DIR:-$HOME/Desktop/output/hackathon/resources/mo_annotation}"
export EDGE_SOURCE="${EDGE_SOURCE:-mnet}"
# Fixed build date for files that embed one (the LaTeX PDF), so reruns are byte-identical.
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1790380800}"
REFRESH_ONLINE="${REFRESH_ONLINE:-0}"
mkdir -p "$HACK_OUT/logs" "$HACK_FIG"

# The steps, in dependency order: label, command.
STEPS=(
  "00|Rscript $HERE/00_environment.R"
  "01|Rscript $HERE/01_node_embeddings.R"
  "01b|Rscript $HERE/01b_metabolite_embeddings.R"
  "01c|python3 $HERE/01c_metabolite_ids.py"
  "01d|python3 $HERE/01d_metabolite_classes.py"
  "02|Rscript $HERE/02_string_edges.R"
  "03|Rscript $HERE/03_edge_weights.R"
  "05|Rscript $HERE/05_rhea_metabolite_protein.R"
  "06|Rscript $HERE/06_metabolite_network.R"
  "08|Rscript $HERE/08_metabolite_rule_experiments.R"
  "09|Rscript $HERE/resource/export_feature_lists.R"
  "10|Rscript $HERE/10_plot_arm_networks.R"
  "11|Rscript $HERE/11_plot_edge_difference.R"
  "12|Rscript $HERE/12_normalization_comparison.R"
  "13|Rscript $HERE/13_logfc_descriptive_stats.R"
  "14|Rscript $HERE/14_joint_network.R"
  "07|Rscript $HERE/07_hub_report.R"      # hub report (after 14: it includes the joint network)
  "15|Rscript $HERE/15_joint_network_classes.R"
  "inv1|python3 $HERE/inventory/glygen_protein_inventory.py"
  "inv2|Rscript $HERE/inventory/glygen_motrpac_inventory.R"
  "inv3|Rscript $HERE/inventory/export_phospho_features.R"
  "16|Rscript $HERE/16_annotated_networks.R"
  "17s|Rscript $HERE/17_filter_stats.R"
  "17i|Rscript $HERE/17_interactive_networks.R"
  "18|Rscript $HERE/18_t2d_lipid_classes.R"
  "18d|Rscript $HERE/18_disease_modules.R"
  "18c|Rscript $HERE/18_option_c_graphical_modules.R"
  "18bc|Rscript $HERE/18_option_bc_hybrid.R"
  "19|Rscript $HERE/19_t2d_stories.R"       # three T2D stories: endurance vs resistance (figures 19a-c)
  "99v|Rscript $HERE/99_validate_outputs.R"
  "99m|Rscript $HERE/99_manifest.R"
)
FROM="${1:-00}"; TO="${2:-99m}"

# Run the selected range; stop at the first failing step (its log is shown).
running=0
for s in "${STEPS[@]}"; do
  label="${s%%|*}"; cmd="${s#*|}"
  [ "$label" = "$FROM" ] && running=1
  if [ "$running" = 1 ]; then
    if [ "$label" = "01c" ] && [ "$REFRESH_ONLINE" != 1 ] && [ -f "$HACK_OUT/01c_metabolite_ids.csv" ]; then
      echo "[01c] skipped: cached web lookups in 01c_metabolite_ids.csv (REFRESH_ONLINE=1 to redo)"
    else
      start=$(date +%s)
      if ! $cmd > "$HACK_OUT/logs/$label.log" 2>&1; then
        echo "[$label] FAILED — last lines of $HACK_OUT/logs/$label.log:"; tail -20 "$HACK_OUT/logs/$label.log"; exit 1
      fi
      echo "[$label] ok ($(( $(date +%s) - start )) s)"
    fi
  fi
  [ "$label" = "$TO" ] && break
done
[ "$running" = 1 ] || { echo "unknown start step '$FROM'"; exit 1; }
echo "done. Validation: $HACK_OUT/logs/99v.log · reproducibility: $HACK_OUT/logs/99m.log"
