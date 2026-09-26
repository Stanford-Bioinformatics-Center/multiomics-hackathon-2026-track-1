#!/usr/bin/env Rscript
# =====================================================================================================
# 01_node_embeddings.R — STEP 1 OF THE NETWORK: ONE "NODE" PER GENE, EACH WITH TWO RESPONSE VECTORS
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   A network is made of nodes (here: genes) joined by edges (here: known protein interactions, step 2).
#   This script builds the nodes. For every gene it writes down how that gene changed after one bout of
#   exercise, in three tissues, at the RNA and the protein level, at three times after exercise. Those
#   18 numbers are the gene's "embedding vector" (its coordinates). Each gene gets TWO vectors: one for
#   the endurance arm and one for the resistance arm, each computed from that arm's own comparison with
#   the control group. (Both arms are compared with the SAME control group, so their measurement errors
#   are correlated; that correlation is saved for a future arm-comparison test.)
#
# THE DATA AND WHERE IT COMES FROM
#   Source: the R package MotrpacHumanPreSuspensionAnalysis (v0.2.4), from the MoTrPAC consortium's
#   human pre-suspension (pre-COVID) study of sedentary adults. It ships finished differential-analysis
#   tables ("*_DA" objects): for every measured feature, how much it changed and how confident we are.
#   We do NOT re-run any statistics on raw data; we read these published results.
#
# QUALITY CONTROL AND FILTERING DONE UPSTREAM (by the consortium, before we see the data)
#   - Samples: only the first supervised exercise bout (visit code ADU_BAS). Samples flagged by the
#     consortium's outlier review (package object OUTLIERS) are removed before normalisation, for every
#     assay. Reasons recorded there are mostly low RNA quality (RIN < 5), principal-component outliers and
#     blood contamination, plus a few genotype/sex-check failures.
#   - RNA-seq (all 3 tissues): genes kept only if expressed above 0.5 counts per million in at least 10%
#     of samples; counts normalised with TMM (edgeR); precision weights from voom; model fitted with dream
#     (variancePartition, a linear mixed model). Model: ~ 0 + group_timepoint + Batch + BMI +
#     calculatedAge + codedsiteid + pct_umi_dup + RIN + Sex + (1 | pid). "(1 | pid)" means each
#     participant gets their own baseline.
#   - MS proteomics (adipose, muscle; TMT labelling): values are log2 ratios to a pooled reference sample,
#     median-normalised, with technical batch effects removed (limma removeBatchEffect) and features with
#     too many missing values removed. A protein is then kept only if every group x timepoint that has
#     paired data has at least 3 participants measured both before and after exercise
#     (filter_paired_n). Model: ~ 0 + group_timepoint + BMI + calculatedAge + Sex + (1 | pid), dream.
#   - OLINK proteomics (blood): a targeted panel of ~1,400 proteins, QC-normalised; same model as MS.
#   - Every result table reports logFC (log2 fold change), its 95% confidence interval (CI.L, CI.R),
#     degrees of freedom, t, p-value and adjusted p-value.
#
# WHAT WE FILTER OR CHOOSE HERE
#   1. Which comparison: the "delta-delta" contrast, i.e. (exercise group after - before) minus
#      (control group after - before) at the same time. This removes changes that happen to anyone
#      over a day in the clinic (fasting, time of day) and keeps the part caused by exercise.
#      EE-CON = endurance vs control, RE-CON = resistance vs control.
#   2. Which genes (the "universe", 471 genes): a gene is kept only if it is measured in ALL six
#      tissue x layer combinations (RNA and protein in adipose, blood, muscle). Genes are matched across
#      layers by their Entrez gene ID. Blood protein is OLINK, a fixed ~1,400-protein panel, and that is
#      what limits the count to 471.
#   3. One measurement per gene per tissue x layer: if several transcripts or proteins map to the same
#      gene, we keep the one with the highest average value ("AveExpr"). For RNA that is the most highly
#      expressed transcript. For MS proteomics AveExpr is an average log-ratio to the pooled reference,
#      not abundance, so there it is simply a fixed tie-break (20 genes affected). Either way the rule
#      never looks at the exercise response, so it cannot cherry-pick responsive features.
#   4. NORMALISATION (critical QC step, team decision 2026-09-26): each ome's log fold changes are divided
#      by that ome's MAXIMUM absolute log fold change, taken across all tissues, timepoints and both arms
#      (RNA: one divisor for adipose + blood + muscle RNA; protein: one divisor for adipose MS + blood OLINK
#      + muscle MS). Every value then lies between -1 and +1, the single largest change in each ome is
#      exactly +-1, and RNA and protein are on the same -1..+1 footing.
#      Why this matters: the network's edge weights are DOT PRODUCTS of these vectors (step 3), i.e. sums
#      of products across the 16 dimensions. Without a common scale, the ome with the largest raw log fold
#      changes would dominate every sum simply because of its units, not its biology.
#      What it keeps: one divisor per ome is shared by all tissues, timepoints and both arms, so every
#      difference WITHIN an ome (between arms, times and tissues) is preserved exactly; signs and zero are
#      unchanged.
#      What to be aware of: (a) the divisor is set by a single measurement (the most extreme change), so a
#      new or corrected extreme value rescales that whole ome; the value and the feature that sets it are
#      written to 01_scale_factors.csv. (b) Because tissues share an ome's divisor, tissues keep their raw
#      scale differences within an ome (e.g. blood OLINK changes are typically larger than muscle MS
#      protein changes), so those dimensions weigh more in the dot products.
#   5. Missing by design: adipose tissue was biopsied at all three times, but adipose PROTEOMICS was only
#      profiled at 4 h, so its 0.5 h and 24 h columns are empty (NA) for every gene. They are left empty,
#      never filled in.
#   6. Time labels: "0.5h" is the package's post_15_30_45_min bin and "4h" its post_3.5_4_hr bin; the exact
#      collection time within a bin differs between tissues, so the same label is the same window, not
#      the same minute.
#
# EMBEDDING LAYOUT (18 numbers per gene per arm)
#   adipose rna 0.5/4/24 h, adipose prot 0.5/4/24 h, blood rna ..., blood prot ..., muscle rna ..., muscle prot ...
#
# UNCERTAINTY, SAVED FOR A FUTURE ARM-COMPARISON TEST (a bootstrap test was built and removed for now)
#   Every number is an estimate. We also save its standard error (how wide its error bar is) on the same
#   scale, and the correlation between the ERRORS of a gene's endurance and resistance estimates. That
#   correlation (median ~0.6) exists because both are compared against the SAME control group; it is a
#   property of the measurement, not a similarity between the two arms' biological responses.
#
# TECH STACK
#   R 4.4; packages data.table (tables), MotrpacHumanPreSuspensionAnalysis (the data).
#
# OUTPUTS (folder $HACK_OUT, default ~/Desktop/output/hackathon-2026-track1/network; kept outside the repo)
#   01_nodes_EE.csv, 01_nodes_RE.csv            the scaled 18-number vectors (the node tables)
#   01_nodes_EE_raw_logFC.csv, ..._RE_raw_...   the same before scaling
#   01_nodes_EE_se.csv, 01_nodes_RE_se.csv       standard error of every number, same scale
#   01_nodes_arm_corr.csv                        correlation between the errors of a gene's EE and RE estimates
#   01_scale_factors.csv                         the divisor for each ome (max |logFC|) and which value sets it
#   01_nodes_feature_provenance.csv              which transcript/protein was used for each gene
# =====================================================================================================

# Load the two packages quietly: the MoTrPAC data package and data.table (fast tables).
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Where results are written. A user can override it with the HACK_OUT environment variable;
# otherwise it defaults to a folder on this computer's Desktop (results never go inside the code repo).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Create that folder if it does not exist yet (and its parent folders); stay quiet if it already exists.
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# The three tissues, in the order they appear in the vector.
TISSUES <- c("adipose", "blood", "muscle")
# The two molecular layers. RNA tables share one naming pattern; protein tables differ by tissue (below).
OMES    <- c(rna = "TRNSCRPT", prot = NA)                       # prot object is tissue-specific
# The protein table for each tissue: MS proteomics for adipose and muscle, OLINK for blood.
PROT_OBJ <- c(adipose = "ADIPOSE_PROT_PR_DA", blood = "BLOOD_PROT_OL_DA", muscle = "MUSCLE_PROT_PR_DA")
# The three post-exercise times: short labels on the left, the package's own names on the right.
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
# The two comparisons we build vectors for: endurance vs control, resistance vs control.
ARMS <- c("EE-CON", "RE-CON")

# Build a lookup table that says which gene each measured feature (transcript or protein) belongs to.
map <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[
  # keep only the three kinds of feature we use: RNA transcripts, MS proteins, OLINK proteins
  assay %in% c("transcript-rna-seq", "prot-pr", "prot-ol"),
  # keep the columns we need and store them as plain text
  .(assay = as.character(assay), feature_id = as.character(feature_id),
    entrez_gene = as.character(entrez_gene), gene_symbol = as.character(gene_symbol))
# drop features that have no Entrez gene ID (they cannot be matched across layers)
])[!is.na(entrez_gene)]

# ---- load every tissue x layer, delta-delta contrasts at the 3 post-exercise times -----------------
# Make the list of six tissue x layer combinations (adipose-rna, adipose-prot, blood-rna, ...).
cells <- CJ(tissue = TISSUES, ome = names(OMES), sorted = FALSE)
# For each combination, read its results table and stack them all into one long table called `da`.
da <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  # the tissue and layer for this combination
  tis <- cells$tissue[i]; ome <- cells$ome[i]
  # the name of the package table to read (e.g. "MUSCLE_TRNSCRPT_DA" or "BLOOD_PROT_OL_DA")
  obj <- if (ome == "rna") paste0(toupper(tis), "_TRNSCRPT_DA") else PROT_OBJ[[tis]]
  # read that table and keep only the rows we need:
  #   the two arm-vs-control comparisons, plus endurance-vs-resistance (used only for its error bar),
  #   and only the three post-exercise times
  x <- as.data.table(get(obj))[contrast_category %in% c(ARMS, "EE-RE") & Timepoint %in% TPS,
         # the columns we keep: feature, which comparison, which time (short label), the change (logFC),
         # the feature's average value (AveExpr), and the standard error. The standard error is exactly
         # logFC / t, because the package's t statistic is logFC divided by its standard error. (Rebuilding
         # it from the confidence interval would need the model's moderated degrees of freedom, which the
         # table does not store, so that route is used only in the rare case t = 0.)
         .(assay, feature_id, arm = as.character(contrast_category),
           tp = names(TPS)[match(as.character(Timepoint), TPS)], logFC, AveExpr,
           se = fifelse(!is.na(t) & t != 0, logFC / t, (CI.R - CI.L) / (2 * qt(0.975, degrees_of_freedom))))]
  # label every row with its tissue and layer
  x[, `:=`(tissue = tis, ome = ome)]
  # attach the gene each feature belongs to (features without a gene are dropped here)
  merge(x, map, by = c("assay", "feature_id"))
}))
# Set aside the endurance-vs-resistance rows. The package's EE-RE estimate equals EE-CON minus RE-CON
# exactly (checked), so it adds no new change values; we use only its standard error, further below.
arm_diff <- da[arm == "EE-RE"]          # EE-RE = EE-CON - RE-CON exactly; used only for its SE
# Keep only the two arm-vs-control comparisons in the main table.
da <- da[arm %in% ARMS]

# ---- universe: genes present in all 6 tissue x layer combinations ----------------------------------
# For each tissue x layer, collect the set of genes measured there.
genes_by_cell <- da[, .(genes = list(unique(entrez_gene))), by = .(tissue, ome)]
# Keep only genes found in every one of the six sets (the overlap of all six).
universe <- Reduce(intersect, genes_by_cell$genes)
# Safety check: this should be exactly 471 genes; stop with an error if the data ever changes.
stopifnot(length(universe) == 471)
# Throw away rows for genes outside the universe.
da <- da[entrez_gene %in% universe]

# ---- one feature per gene per tissue x layer: the highest average value ("AveExpr") ------------------
# List every candidate feature for each gene in each tissue x layer, with its AveExpr (the same in every
# comparison for a feature; averaging just collapses the duplicates). For RNA, AveExpr is average
# expression, so this picks the most highly expressed transcript. For MS proteomics it is the average
# log-ratio to the pooled reference sample, not abundance; the rule is then an arbitrary but fixed
# tie-break. Either way it never looks at the exercise response.
feat_expr <- unique(da[, .(tissue, ome, entrez_gene, feature_id, AveExpr)])[
  , .(AveExpr = mean(AveExpr)), by = .(tissue, ome, entrez_gene, feature_id)]
# Sort so the highest-AveExpr feature comes first for each gene (ties broken by feature name, so the
# choice is always the same on every run).
setorder(feat_expr, tissue, ome, entrez_gene, -AveExpr, feature_id)
# Record how many candidates each gene had (for the provenance file).
feat_expr[, n_candidates := .N, by = .(tissue, ome, entrez_gene)]
# Keep the first (highest-AveExpr) feature for each gene in each tissue x layer.
chosen <- feat_expr[, .SD[1], by = .(tissue, ome, entrez_gene)]
# Keep only the rows of the chosen features in the main table.
da <- da[chosen[, .(tissue, ome, entrez_gene, feature_id)], on = .(tissue, ome, entrez_gene, feature_id)]
# Safety check: now there must be exactly one value per gene, tissue, layer, arm and time.
stopifnot(!anyDuplicated(da[, .(tissue, ome, entrez_gene, arm, tp)]))

# ---- correlation between the EE-CON and RE-CON estimates (they share the control group) -----------
# Standard statistics: the variance of a difference is Var(EE) + Var(RE) - 2 x Cov(EE, RE).
# We know all three standard errors, so we can solve for the covariance and turn it into a correlation:
#   rho = (se_EE^2 + se_RE^2 - se_diff^2) / (2 se_EE se_RE)
# Keep the EE-RE rows for exactly the features chosen above.
arm_diff <- arm_diff[chosen[, .(tissue, ome, entrez_gene, feature_id)], on = .(tissue, ome, entrez_gene, feature_id), nomatch = NULL]
# Put the EE and RE standard errors side by side, one row per gene x tissue x layer x time.
se_w <- dcast(da, tissue + ome + entrez_gene + tp ~ arm, value.var = "se")
# Add the standard error of the EE-RE difference next to them.
se_w <- arm_diff[, .(tissue, ome, entrez_gene, tp, se_diff = se)][se_w, on = .(tissue, ome, entrez_gene, tp)]
# Safety checks: every row has its difference standard error, and no row was duplicated by the join.
stopifnot(!anyNA(se_w$se_diff), !anyDuplicated(se_w[, .(tissue, ome, entrez_gene, tp)]))
# Compute the correlation with the formula above; clamp it into the valid range -1..1 (rounding safety).
se_w[, rho := pmin(1, pmax(-1, (`EE-CON`^2 + `RE-CON`^2 - se_diff^2) / (2 * `EE-CON` * `RE-CON`)))]
# Attach the correlation to the main table.
da <- se_w[, .(tissue, ome, entrez_gene, tp, rho)][da, on = .(tissue, ome, entrez_gene, tp)]

# ---- normalisation (critical QC): divide each ome by its maximum absolute log fold change -------------
# For each ome (RNA, protein), one divisor: the largest absolute logFC over all 471 genes, all tissues,
# all timepoints and BOTH arms together; also record which gene / tissue / time / arm sets it.
scale_f <- da[, .(max_abs_logFC = max(abs(logFC)),
                  set_by = paste(gene_symbol, tissue, tp, arm)[which.max(abs(logFC))],
                  n_values = .N), by = ome]
# Attach each ome's divisor to its rows.
da <- scale_f[, .(ome, max_abs_logFC)][da, on = "ome"]
# Divide every change, and its standard error, by the ome's divisor (values now lie in -1..+1; signs and
# zero are unchanged).
da[, `:=`(scaled = logFC / max_abs_logFC, scaled_se = se / max_abs_logFC)]
# Safety check: the largest absolute scaled value in each ome is exactly 1.
stopifnot(all(abs(da[, max(abs(scaled)), by = ome]$V1 - 1) < 1e-12))
# Save the two divisors (and what sets them) so anyone can undo or check the normalisation.
fwrite(scale_f, file.path(OUT, "01_scale_factors.csv"))
# Show them on screen.
print(scale_f)

# ---- 18-number vector, one file per arm ---------------------------------------------------------------
# The 18 column names in their fixed order: tissue, then layer, then time (e.g. "muscle_prot_24h").
dims <- CJ(tissue = TISSUES, ome = names(OMES), tp = names(TPS), sorted = FALSE)[
  , dim := paste(tissue, ome, tp, sep = "_")]
# One gene symbol (human-readable name, e.g. "HSPB1") per Entrez ID.
symbols <- unique(map[entrez_gene %in% universe, .(entrez_gene, gene_symbol)])[
  , .(gene_symbol = gene_symbol[1]), by = entrez_gene]

# Helper: turn the long table into a wide one (one row per gene, 18 columns) for one arm and one value.
to_wide <- function(a, value) {
  # rows for this arm, lined up against the 18 column names
  long <- da[arm == a][dims, on = .(tissue, ome, tp)]
  # spread into one column per dimension
  wide <- dcast(long[!is.na(entrez_gene)], entrez_gene ~ dim, value.var = value)
  # add back the columns that exist for no gene (adipose protein 0.5 h and 24 h) as empty (NA)
  for (d in setdiff(dims$dim, names(wide))) set(wide, j = d, value = NA_real_)
  # add the gene symbol
  wide <- symbols[wide, on = "entrez_gene"]
  # put the columns in the fixed order
  setcolorder(wide, c("entrez_gene", "gene_symbol", dims$dim))
  # sort rows alphabetically by gene symbol (same order in every file)
  setorder(wide, gene_symbol)
  # safety check: 471 genes and 2 + 18 columns
  stopifnot(nrow(wide) == 471, ncol(wide) == 2 + 18)
  wide
}

# Write the files for each arm.
for (a in ARMS) {
  # short tag for file names: "EE" or "RE"
  tag <- sub("-CON", "", a)
  # the main node table: scaled values
  f <- file.path(OUT, sprintf("01_nodes_%s.csv", tag))
  fwrite(to_wide(a, "scaled"), f)
  # the same numbers before scaling, for reference
  fwrite(to_wide(a, "logFC"), file.path(OUT, sprintf("01_nodes_%s_raw_logFC.csv", tag)))
  # the standard error of every scaled number
  fwrite(to_wide(a, "scaled_se"), file.path(OUT, sprintf("01_nodes_%s_se.csv", tag)))
  # report progress on screen
  message(sprintf("%s: 471 nodes x 18 dims -> %s", a, f))
}

# The EE-RE correlation belongs to the gene x dimension, not to one arm, so it is written once.
fwrite(to_wide("EE-CON", "rho"), file.path(OUT, "01_nodes_arm_corr.csv"))   # same for both arms

# Provenance: for each gene, which feature was used in each tissue x layer and how many it was chosen from.
prov <- dcast(chosen, entrez_gene ~ paste(tissue, ome, sep = "_"),
              value.var = c("feature_id", "n_candidates"))
# add gene symbols
prov <- symbols[prov, on = "entrez_gene"]
# save and report
fwrite(prov, file.path(OUT, "01_nodes_feature_provenance.csv"))
# Report where the provenance file was written.
message("provenance -> ", file.path(OUT, "01_nodes_feature_provenance.csv"))
