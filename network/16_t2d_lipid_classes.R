#!/usr/bin/env Rscript
# =====================================================================================================
# 16_t2d_lipid_classes.R — STEP 16: HOW DO T2D-RELEVANT LIPID CLASSES RESPOND TO EACH EXERCISE ARM?
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   In the metabolite networks, sphingolipids and fatty acyls behave differently between endurance and
#   resistance exercise. Both classes are elevated in type 2 diabetes (T2D) and insulin resistance. Does
#   either arm move these classes in the opposite (lower) direction, and where and when? This script
#   produces the tables behind the README section "Preliminary observations"; it adds no new statistics
#   beyond the consortium's estimates.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Per arm, tissue and time: the mean log fold change of each T2D-relevant metabolite class among our
#      450 metabolites (free fatty acids, acylcarnitines = "fatty esters", ceramides, sphingomyelins =
#      "phosphosphingolipids"), and how many individual metabolites change clearly down or up.
#   2. The same, metabolite by metabolite, for T2D-relevant species: ceramides C16:0 / C18:0 / C22:0,
#      palmitic, myristic and linoleic acid.
#   3. The consortium's clinical chemistry for these participants: total plasma free fatty acids (NEFA),
#      glycerol and lactate per arm and time (an independent assay from the metabolomics).
#   "Clearly" = |z| > 2, where z = log fold change / its standard error. This is a descriptive,
#   UNCORRECTED marker of values that stand out from their own measurement noise, not a test.
#
# HOW TO RUN
#   After steps 1b and 1c:   Rscript network/16_t2d_lipid_classes.R      (a few seconds)
#
# DATA AND PROVENANCE
#   MoTrPAC human pre-suspension study, R package MotrpacHumanPreSuspensionAnalysis v0.2.4.
#   Metabolites: unnormalised delta-delta log2 fold changes (arm change from pre-exercise minus control
#   change) and their standard errors from step 1b; RefMet classes from step 1c. Clinical chemistry: the
#   package's CLIN_CHEMISTRY_DA table (EE-CON and RE-CON contrasts; its adjusted p-values are the
#   consortium's). Upstream QC is described in the headers of steps 1 and 1b and in the README.
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis.
#
# INPUTS (files and columns)
#   $HACK_OUT/01b_metab_nodes_{EE,RE}_raw_logFC.csv   metabolite + 9 columns <tissue>_metab_<time> (log2 FC)
#   $HACK_OUT/01b_metab_nodes_{EE,RE}.csv, ..._se.csv normalised values and their standard errors (z = ratio)
#   $HACK_OUT/01c_metabolite_ids.csv                  metabolite, super_class, main_class, ...
#   package object CLIN_CHEMISTRY_DA                  feature_id, contrast_category, Timepoint, logFC, adj_p_value
#
# OUTPUTS (files, locations, columns)
#   $HACK_OUT/16_t2d_class_summary.csv     super_class, main_class, tissue, time, arm, n, mean_logFC,
#                                          n_clearly_down, n_clearly_up
#   $HACK_OUT/16_t2d_species.csv           metabolite, tissue, time, arm, logFC, z
#   $HACK_OUT/16_clinical_nefa_lactate.csv feature (NEFA, Glycerol, Lactate), arm, time, logFC, adj_p
#   Nothing is written inside the repository.
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   Plasma free fatty acids at 0.5 h: mean logFC endurance +0.02, resistance -0.24 (29 clearly down);
#   clinical NEFA at 0.5 h: endurance +0.02, resistance -0.58 (adj p ~5e-8). The script stops if any of
#   the 450 metabolites is missing, or if a requested class or species is absent.
#
# KNOWN LIMITS
#   Descriptive only: differences between the arms are not tested (the bootstrap test was removed).
#   Acylcarnitines rise during exercise because fat oxidation increases; that is not the same as the
#   incomplete oxidation associated with T2D. A single bout in sedentary adults says nothing directly about
#   the long-term effects of training, which is where exercise's protection against T2D comes from.
# =====================================================================================================

# Load the packages quietly.
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# The T2D-relevant RefMet main classes (all elevated in T2D / insulin resistance in the literature).
CLASSES <- c("Fatty acids", "Fatty esters", "Ceramides", "Phosphosphingolipids")
# The T2D-relevant individual species.
SPECIES <- c("Cer 18:1;O2/16:0", "Cer 18:1;O2/18:0", "Cer 18:1;O2/22:0", "Palmitic acid", "Myristic acid", "Linoleic acid")
# "Clearly" changed: |z| above this (descriptive, uncorrected).
Z_CUT <- 2

# ---- one long table: metabolite x tissue x time x arm, with logFC and z -----------------------------
# Helper: read one arm's raw logFC, normalised values and standard errors, and combine them.
arm_long <- function(arm) {
  # raw log2 fold changes
  r <- melt(fread(file.path(OUT, sprintf("01b_metab_nodes_%s_raw_logFC.csv", arm))), id.vars = "metabolite", value.name = "logFC")
  # normalised values
  n <- melt(fread(file.path(OUT, sprintf("01b_metab_nodes_%s.csv", arm))), id.vars = "metabolite", value.name = "norm")
  # their standard errors (same scale as the normalised values)
  s <- melt(fread(file.path(OUT, sprintf("01b_metab_nodes_%s_se.csv", arm))), id.vars = "metabolite", value.name = "se")
  # combine; z = value / standard error (the same whether computed on raw or normalised values)
  x <- merge(merge(r, n, by = c("metabolite", "variable")), s, by = c("metabolite", "variable"))
  # label arm, tissue and time from the column name (e.g. "blood_metab_0.5h")
  x[, `:=`(arm = arm, z = norm / se, tissue = sub("_.*", "", variable), time = sub(".*_", "", variable))]
  x[, .(metabolite, tissue, time, arm, logFC, z)]
}
# Both arms together.
d <- rbind(arm_long("EE"), arm_long("RE"))
# Attach RefMet classes from step 1c.
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))[, .(metabolite, super_class, main_class)]
# (join on the metabolite name)
d <- merge(d, ids, by = "metabolite")
# Safety check: all 450 metabolites are present, and every requested class and species exists.
stopifnot(uniqueN(d$metabolite) == 450, all(CLASSES %in% d$main_class), all(SPECIES %in% d$metabolite))

# ---- table 1: class summary per arm, tissue and time -------------------------------------------------
cls <- d[main_class %in% CLASSES,
         .(n = .N, mean_logFC = round(mean(logFC), 3), n_clearly_down = sum(z < -Z_CUT), n_clearly_up = sum(z > Z_CUT)),
         by = .(super_class, main_class, tissue, time, arm)]
# Order rows by class, tissue, time and arm.
cls <- cls[order(super_class, main_class, tissue, match(time, c("0.5h", "4h", "24h")), arm)]
# Save.
fwrite(cls, file.path(OUT, "16_t2d_class_summary.csv"))

# ---- table 2: T2D-relevant species -------------------------------------------------------------------
sp <- d[metabolite %in% SPECIES, .(metabolite, tissue, time, arm, logFC = round(logFC, 3), z = round(z, 2))]
# Order rows by species, tissue, time and arm.
sp <- sp[order(metabolite, tissue, match(time, c("0.5h", "4h", "24h")), arm)]
# Save.
fwrite(sp, file.path(OUT, "16_t2d_species.csv"))

# ---- table 3: clinical chemistry (independent assay) --------------------------------------------------
cc <- as.data.table(CLIN_CHEMISTRY_DA)[feature_id %in% c("NEFA", "Glycerol", "Lactate") & contrast_category %in% c("EE-CON", "RE-CON"),
        .(feature = feature_id, arm = sub("-CON", "", as.character(contrast_category)), time = as.character(Timepoint),
          logFC = round(logFC, 3), adj_p = signif(adj_p_value, 2))]
# Order rows by feature, time (as the study sequence) and arm.
TORDER <- c("pre_exercise", "during_20_min", "during_40_min", "post_10_min", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
# (apply that order)
cc <- cc[order(feature, match(time, TORDER), arm)]
# Save.
fwrite(cc, file.path(OUT, "16_clinical_nefa_lactate.csv"))

# ---- headline numbers on screen ----------------------------------------------------------------------
# Plasma free fatty acids at 0.5 h, per arm.
print(cls[main_class == "Fatty acids" & tissue == "blood" & time == "0.5h"])
# Clinical NEFA at 0.5 h, per arm.
print(cc[feature == "NEFA" & time == "post_15_30_45_min"])
