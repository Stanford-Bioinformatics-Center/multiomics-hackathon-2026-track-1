#!/usr/bin/env Rscript
# =====================================================================================================
# 01b_metabolite_embeddings.R — STEP 1b: ONE NODE PER METABOLITE, EACH WITH TWO 9-NUMBER RESPONSE VECTORS
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The metabolite counterpart of step 1 (01_node_embeddings.R, which does genes). For every metabolite
#   (a small molecule such as an amino acid, a lipid or a sugar) it writes down how that metabolite
#   changed after one bout of exercise, in three tissues, at three times after exercise. Metabolomics is
#   ONE molecular layer (genes had two: RNA and protein), so each vector has 3 tissues x 3 times = 9
#   numbers, half the length of a gene's 18. Each metabolite gets TWO vectors, one for endurance and one
#   for resistance, built completely separately so the two can later be compared honestly.
#
# THE DATA AND WHERE IT COMES FROM
#   Source: the R package MotrpacHumanPreSuspensionAnalysis (v0.2.4), MoTrPAC human pre-suspension
#   (pre-COVID) study. Tables ADIPOSE_METAB_DA, BLOOD_METAB_DA (blood = plasma), MUSCLE_METAB_DA: finished
#   differential-analysis results. We do NOT re-run statistics on raw data.
#
# QUALITY CONTROL AND FILTERING DONE UPSTREAM (by the consortium, before we see the data)
#   - Samples: first supervised exercise bout only (visit ADU_BAS); outlier samples from the
#     consortium's review (package object OUTLIERS) removed.
#   - Platforms: each tissue was measured on 10-12 assays ("platforms"): targeted panels (known compounds,
#     e.g. amino acids, TCA-cycle acids) and untargeted LC-MS (HILIC, ion-pairing, reversed-phase, lipid
#     reversed-phase, positive and negative mode).
#   - Per platform: features missing in more than 20% of samples removed; remaining gaps filled by
#     k-nearest-neighbour imputation; log2 transform; median/MAD normalisation (untargeted platforms).
#   - Names are standardised to RefMet (the Metabolomics Workbench reference nomenclature).
#   - One platform per metabolite per tissue: if a metabolite was measured on several platforms, the
#     one with the lowest technical coefficient of variation (CV, measured on repeated QC samples; the
#     most reproducible) is kept (package function .prioritize_metab_by_cv_da, table METABOLOMICS_CVS).
#     We checked: every metabolite appears on exactly one platform per tissue.
#   - Model (dream linear mixed model, one per metabolite):
#       targeted:   ~ 0 + group_timepoint + BMI + calculatedAge + codedsiteid + Sex + (1 | pid)
#       untargeted: the same + raw_intensity_post (each sample's overall signal level, a technical covariate)
#     "(1 | pid)" gives each participant their own baseline.
#
# WHAT WE FILTER OR CHOOSE HERE
#   1. Which comparison: the "delta-delta" contrast, (exercise group after - before) minus (control group
#      after - before) at the same time; EE-CON = endurance vs control, RE-CON = resistance vs control.
#   2. Which metabolites (450): measured in ALL three tissues. Tissues are matched by EXACT RefMet name.
#      Names are not lower-cased or stripped of punctuation, because in lipid names punctuation carries
#      meaning (e.g. "PC 16:0_18:1" vs "PC 16:0/18:1"); loosening the match would merge different molecules.
#   3. No duplicate handling is needed: the consortium already kept one platform per metabolite (above).
#   4. Scaling: each tissue's metabolomics block is divided by one number, its root-mean-square logFC,
#      pooled over the 450 metabolites, all three times and BOTH arms. One shared number per tissue is a
#      common unit: it puts tissues on one footing while keeping real differences between arms and between
#      times intact (a separate number per arm could hide a genuine arm difference).
#   5. Nothing is missing by design: all three tissues have all three times in both arms.
#
# EMBEDDING LAYOUT (9 numbers per metabolite per arm)
#   adipose 0.5/4/24 h, blood 0.5/4/24 h, muscle 0.5/4/24 h
#
# UNCERTAINTY, SAVED FOR LATER COMPARISON STEPS
#   Each number's standard error (same scale) and the correlation between a metabolite's endurance and
#   resistance estimates (correlated because both arms share one control group), exactly as for genes.
#
# TECH STACK
#   R 4.4; data.table (tables), MotrpacHumanPreSuspensionAnalysis (the data).
#
# OUTPUTS (folder $HACK_OUT, default ~/Desktop/output/hackathon-2026-track1/network; outside the repo)
#   01b_metab_nodes_EE.csv, 01b_metab_nodes_RE.csv      scaled 9-number vectors (the node tables)
#   01b_metab_nodes_EE_raw_logFC.csv, ..._RE_...        the same before scaling
#   01b_metab_nodes_EE_se.csv, 01b_metab_nodes_RE_se.csv standard errors, same scale
#   01b_metab_nodes_arm_corr.csv                         correlation between EE and RE estimates
#   01b_metab_scale_factors.csv                          the one divisor per tissue
#   01b_metab_nodes_provenance.csv                       which platform measured each metabolite in each tissue
# =====================================================================================================

# Load the two packages quietly: the MoTrPAC data package and data.table (fast tables).
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Where results are written (override with HACK_OUT); results never go inside the code repo.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Create the folder if needed; stay quiet if it already exists.
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# The three tissues, in the order they appear in the vector.
TISSUES <- c("adipose", "blood", "muscle")
# The three post-exercise times: short labels on the left, the package's own names on the right.
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
# The two comparisons we build vectors for: endurance vs control, resistance vs control.
ARMS <- c("EE-CON", "RE-CON")

# ---- load each tissue's metabolomics results at the 3 post-exercise times -----------------------------
# For each tissue, read its results table and stack the three into one long table called `da`.
da <- rbindlist(lapply(TISSUES, function(tis) {
  # read e.g. MUSCLE_METAB_DA and keep only the rows we need:
  #   the two arm-vs-control comparisons, plus endurance-vs-resistance (used only for its error bar),
  #   and only the three post-exercise times
  x <- as.data.table(get(paste0(toupper(tis), "_METAB_DA")))[
    contrast_category %in% c(ARMS, "EE-RE") & Timepoint %in% TPS,
    # the columns we keep: metabolite name, platform, comparison, time (short label), the change
    # (logFC), and the standard error. The package gives a 95% confidence interval, so the standard
    # error is its half-width divided by the matching t-distribution value.
    .(metabolite = as.character(feature_id), platform = as.character(platform),
      arm = as.character(contrast_category), tp = names(TPS)[match(as.character(Timepoint), TPS)],
      logFC, se = (CI.R - CI.L) / (2 * qt(0.975, degrees_of_freedom)))]
  # label every row with its tissue
  x[, tissue := tis]
}))
# Safety check: one platform per metabolite per tissue (the consortium's lowest-CV rule), so exactly one
# row per metabolite x tissue x comparison x time.
stopifnot(!anyDuplicated(da[, .(tissue, metabolite, arm, tp)]))
# Set aside the endurance-vs-resistance rows. EE-RE equals EE-CON minus RE-CON exactly, so it adds no new
# change values; we use only its standard error, below.
arm_diff <- da[arm == "EE-RE"]
# Keep only the two arm-vs-control comparisons in the main table.
da <- da[arm %in% ARMS]

# ---- universe: metabolites measured in all three tissues, both arms, all three times ----------------
# For each tissue, the metabolites present in all 6 arm x time results (2 arms x 3 times).
per_tissue <- da[, .(k = uniqueN(paste(arm, tp))), by = .(tissue, metabolite)][k == 6]
# Keep metabolites found in every tissue (exact RefMet-name match).
universe <- Reduce(intersect, split(per_tissue$metabolite, per_tissue$tissue))
# Safety check: this should be exactly 450; stop with an error if the data ever changes.
stopifnot(length(universe) == 450)
# Throw away rows for metabolites outside the universe.
da <- da[metabolite %in% universe]
arm_diff <- arm_diff[metabolite %in% universe]

# ---- correlation between the EE-CON and RE-CON estimates (they share the control group) -----------
# Var(EE - RE) = Var(EE) + Var(RE) - 2 Cov(EE, RE), so with all three standard errors known:
#   rho = (se_EE^2 + se_RE^2 - se_diff^2) / (2 se_EE se_RE)
# Put the EE and RE standard errors side by side, one row per metabolite x tissue x time.
se_w <- dcast(da, tissue + metabolite + tp ~ arm, value.var = "se")
# Add the standard error of the EE-RE difference next to them.
se_w <- arm_diff[, .(tissue, metabolite, tp, se_diff = se)][se_w, on = .(tissue, metabolite, tp)]
# Safety check: every row has its difference standard error.
stopifnot(!anyNA(se_w$se_diff))
# Compute the correlation; clamp into -1..1 (rounding safety).
se_w[, rho := pmin(1, pmax(-1, (`EE-CON`^2 + `RE-CON`^2 - se_diff^2) / (2 * `EE-CON` * `RE-CON`)))]
# Attach the correlation to the main table.
da <- se_w[, .(tissue, metabolite, tp, rho)][da, on = .(tissue, metabolite, tp)]

# ---- per-tissue scaling: one shared unit per tissue --------------------------------------------------
# The root-mean-square logFC of each tissue, over all 450 metabolites, all 3 times and BOTH arms.
scale_f <- da[, .(rms_logFC = sqrt(mean(logFC^2)), n_values = .N), by = tissue]
# Attach each tissue's divisor to its rows.
da <- scale_f[, .(tissue, rms_logFC)][da, on = "tissue"]
# Divide every change, and its standard error, by the tissue's divisor (signs and zero unchanged).
da[, `:=`(scaled = logFC / rms_logFC, scaled_se = se / rms_logFC)]
# Save the three divisors so anyone can undo or check the scaling, and show them.
fwrite(scale_f, file.path(OUT, "01b_metab_scale_factors.csv"))
print(scale_f)

# ---- 9-number vector, one file per arm --------------------------------------------------------------
# The 9 column names in their fixed order: tissue, then time (e.g. "blood_metab_4h").
dims <- CJ(tissue = TISSUES, tp = names(TPS), sorted = FALSE)[, dim := paste(tissue, "metab", tp, sep = "_")]

# Helper: turn the long table into a wide one (one row per metabolite, 9 columns) for one arm and value.
to_wide <- function(a, value) {
  # rows for this arm, labelled with their column name
  long <- da[arm == a][dims, on = .(tissue, tp)]
  # spread into one column per dimension
  wide <- dcast(long, metabolite ~ dim, value.var = value)
  # put the columns in the fixed order
  setcolorder(wide, c("metabolite", dims$dim))
  # sort rows alphabetically by metabolite name (same order in every file)
  setorder(wide, metabolite)
  # safety check: 450 metabolites, 1 + 9 columns, nothing missing
  stopifnot(nrow(wide) == 450, ncol(wide) == 1 + 9, !anyNA(wide))
  wide
}

# Write the files for each arm.
for (a in ARMS) {
  # short tag for file names: "EE" or "RE"
  tag <- sub("-CON", "", a)
  # the main node table: scaled values
  f <- file.path(OUT, sprintf("01b_metab_nodes_%s.csv", tag))
  fwrite(to_wide(a, "scaled"), f)
  # the same numbers before scaling, for reference
  fwrite(to_wide(a, "logFC"), file.path(OUT, sprintf("01b_metab_nodes_%s_raw_logFC.csv", tag)))
  # the standard error of every scaled number
  fwrite(to_wide(a, "scaled_se"), file.path(OUT, sprintf("01b_metab_nodes_%s_se.csv", tag)))
  # report progress
  message(sprintf("%s: 450 metabolites x 9 dims -> %s", a, f))
}

# The EE-RE correlation belongs to the metabolite x dimension, not to one arm, so it is written once.
fwrite(to_wide("EE-CON", "rho"), file.path(OUT, "01b_metab_nodes_arm_corr.csv"))

# Provenance: which platform measured each metabolite in each tissue.
prov <- dcast(unique(da[, .(metabolite, tissue, platform)]), metabolite ~ tissue, value.var = "platform")
setnames(prov, TISSUES, paste0("platform_", TISSUES))
fwrite(prov, file.path(OUT, "01b_metab_nodes_provenance.csv"))
message("provenance -> ", file.path(OUT, "01b_metab_nodes_provenance.csv"))
