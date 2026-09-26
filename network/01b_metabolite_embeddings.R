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
#   for resistance, each computed from that arm's own comparison with the control group. (Both arms share
#   the same control group, so their measurement errors are correlated; that correlation is saved.)
#
# THE DATA AND WHERE IT COMES FROM
#   Source: the R package MotrpacHumanPreSuspensionAnalysis (v0.2.4), MoTrPAC human pre-suspension
#   (pre-COVID) study. Tables ADIPOSE_METAB_DA, BLOOD_METAB_DA (blood = plasma), MUSCLE_METAB_DA: finished
#   differential-analysis results. We do NOT re-run statistics on raw data.
#
# QUALITY CONTROL AND FILTERING DONE UPSTREAM (by the consortium, before we see the data)
#   - Samples: first supervised exercise bout only (visit ADU_BAS); outlier samples from the
#     consortium's review (package object OUTLIERS) removed, per platform, before normalisation.
#   - Platforms: each tissue was measured on 10-12 assays ("platforms"): targeted panels (known compounds,
#     e.g. amino acids, TCA-cycle acids) and untargeted LC-MS (HILIC, ion-pairing, reversed-phase, lipid
#     reversed-phase, positive and negative mode).
#   - Per platform: zero or negative intensities set to missing; features missing in more than 20% of
#     samples removed; remaining gaps imputed (a minimum-value / k-nearest-neighbour hybrid; small targeted
#     panels of <= 12 features use half-minimum imputation); log2 transform; median/MAD normalisation for
#     untargeted platforms, applied only where overall sample intensity is not associated with sex or group.
#   - Names are standardised to RefMet (the Metabolomics Workbench reference nomenclature).
#   - One platform per metabolite per tissue: if a metabolite was measured on several platforms, the one
#     with the lowest technical coefficient of variation (CV; the most reproducible, per the consortium's
#     METABOLOMICS_CVS table) is kept. We checked it from the data: every metabolite appears on exactly one
#     platform per tissue, and for the 114 / 219 / 207 metabolites (adipose / blood / muscle) measured on
#     several platforms the kept platform is always the lowest-CV one.
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
#   4. NORMALISATION (critical QC step, team decision 2026-09-26): all metabolite log fold changes are
#      divided by the MAXIMUM absolute log fold change across all three tissues, all three times and both
#      arms (metabolomics is one ome, so one divisor). Values then lie between -1 and +1, the single
#      largest change is exactly +-1, and metabolites are on the same -1..+1 footing as the gene
#      embeddings. This matters because metabolite edge weights are dot products (step 6); a common,
#      shared divisor keeps every difference between arms, times and tissues intact. The divisor is set by
#      one measurement (recorded in 01b_metab_scale_factors.csv), and tissues keep their raw scale
#      differences relative to each other.
#   5. Nothing is missing by design: all three tissues have all three times in both arms.
#   6. Time labels: "0.5h" is the package's post_15_30_45_min bin and "4h" its post_3.5_4_hr bin; the exact
#      collection time within a bin differs between tissues.
#
# EMBEDDING LAYOUT (9 numbers per metabolite per arm)
#   adipose 0.5/4/24 h, blood 0.5/4/24 h, muscle 0.5/4/24 h
#
# UNCERTAINTY, SAVED FOR LATER COMPARISON STEPS
#   Each number's standard error (same scale) and the correlation between the ERRORS of a metabolite's
#   endurance and resistance estimates (correlated because both arms share one control group; a property
#   of the measurement, not a similarity of responses), exactly as for genes.
#
# TECH STACK
#   R 4.4; data.table (tables), MotrpacHumanPreSuspensionAnalysis (the data).
#
# OUTPUTS (folder $HACK_OUT, default ~/Desktop/output/hackathon-2026-track1/network; outside the repo)
#   01b_metab_nodes_EE.csv, 01b_metab_nodes_RE.csv      scaled 9-number vectors (the node tables)
#   01b_metab_nodes_EE_raw_logFC.csv, ..._RE_...        the same before scaling
#   01b_metab_nodes_EE_se.csv, 01b_metab_nodes_RE_se.csv standard errors, same scale
#   01b_metab_nodes_arm_corr.csv                         correlation between the errors of EE and RE estimates
#   01b_metab_scale_factors.csv                          the divisor (max |logFC|) and which value sets it
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
    # (logFC), and the standard error. The standard error is exactly logFC / t (the package's t statistic
    # is logFC divided by its standard error); the confidence-interval route is used only if t = 0.
    .(metabolite = as.character(feature_id), platform = as.character(platform),
      arm = as.character(contrast_category), tp = names(TPS)[match(as.character(Timepoint), TPS)],
      logFC, se = fifelse(!is.na(t) & t != 0, logFC / t, (CI.R - CI.L) / (2 * qt(0.975, degrees_of_freedom))))]
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
# Same for the set-aside endurance-vs-resistance rows.
arm_diff <- arm_diff[metabolite %in% universe]

# ---- correlation between the EE-CON and RE-CON estimates (they share the control group) -----------
# Var(EE - RE) = Var(EE) + Var(RE) - 2 Cov(EE, RE), so with all three standard errors known:
#   rho = (se_EE^2 + se_RE^2 - se_diff^2) / (2 se_EE se_RE)
# Put the EE and RE standard errors side by side, one row per metabolite x tissue x time.
se_w <- dcast(da, tissue + metabolite + tp ~ arm, value.var = "se")
# Add the standard error of the EE-RE difference next to them.
se_w <- arm_diff[, .(tissue, metabolite, tp, se_diff = se)][se_w, on = .(tissue, metabolite, tp)]
# Safety checks: every row has its difference standard error, and no row was duplicated by the join.
stopifnot(!anyNA(se_w$se_diff), !anyDuplicated(se_w[, .(tissue, metabolite, tp)]))
# Compute the correlation; clamp into -1..1 (rounding safety).
se_w[, rho := pmin(1, pmax(-1, (`EE-CON`^2 + `RE-CON`^2 - se_diff^2) / (2 * `EE-CON` * `RE-CON`)))]
# Attach the correlation to the main table.
da <- se_w[, .(tissue, metabolite, tp, rho)][da, on = .(tissue, metabolite, tp)]

# ---- normalisation (critical QC): divide by the maximum absolute log fold change ---------------------
# One divisor for metabolomics: the largest absolute logFC over all 450 metabolites, all tissues, all
# times and BOTH arms; also record which metabolite / tissue / time / arm sets it.
scale_f <- da[, .(ome = "metab", max_abs_logFC = max(abs(logFC)),
                  set_by = paste(metabolite, tissue, tp, arm)[which.max(abs(logFC))], n_values = .N)]
# Attach the divisor to every row.
da[, max_abs_logFC := scale_f$max_abs_logFC]
# Divide every change, and its standard error, by the divisor (values now lie in -1..+1; signs unchanged).
da[, `:=`(scaled = logFC / max_abs_logFC, scaled_se = se / max_abs_logFC)]
# Safety check: the largest absolute scaled value is exactly 1.
stopifnot(abs(max(abs(da$scaled)) - 1) < 1e-12)
# Save the divisor (and what sets it) so anyone can undo or check the normalisation.
fwrite(scale_f, file.path(OUT, "01b_metab_scale_factors.csv"))
# Show it on screen.
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
# Rename the columns to platform_adipose, platform_blood, platform_muscle.
setnames(prov, TISSUES, paste0("platform_", TISSUES))
# Save the provenance table.
fwrite(prov, file.path(OUT, "01b_metab_nodes_provenance.csv"))
# Report where it was written.
message("provenance -> ", file.path(OUT, "01b_metab_nodes_provenance.csv"))
