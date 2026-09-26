#!/usr/bin/env Rscript
# =====================================================================================================
# 13_logfc_descriptive_stats.R — STEP 13: DESCRIPTIVE STATISTICS OF THE LOG FOLD CHANGES, AS A LATEX PDF
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   How large are the exercise-induced log fold changes in each molecular layer ("ome"), and do they differ
#   between the two study arms? The answer justifies the normalisation used before the network dot
#   products (README, "Critical QC step"): the omes have clearly different ranges, so without a common
#   scale one ome would dominate the edge weights through its units alone.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Summarises the log fold changes that the networks are built from, for each ome (RNA, protein,
#   metabolites), in two tables:
#     Table 1: pooled across the two study arms (endurance vs control and resistance vs control together);
#     Table 2: by study arm (each arm separately).
#   For each: number of values, minimum, maximum, mean and standard deviation (SD) of the log fold change.
#   The tables are written as a LaTeX document and compiled to PDF, and also saved as a CSV.
#
# HOW TO RUN
#   After steps 1 and 1b have been run (they write the raw log fold change tables):
#     Rscript network/13_logfc_descriptive_stats.R
#   Runtime: a few seconds (the first run may take longer while TinyTeX installs LaTeX packages).
#
# DATA AND PROVENANCE
#   The UNNORMALISED log fold changes (log2 scale) of the network features, exactly as they enter the step 1
#   normalisation: MoTrPAC human pre-suspension study, R package MotrpacHumanPreSuspensionAnalysis v0.2.4,
#   delta-delta contrasts (arm change from pre-exercise minus control change) at 0.5 / 4 / 24 h, in
#   adipose, blood and muscle, for
#     - RNA: the 471 network genes (RNA-seq, 3 tissues);
#     - protein: the same 471 genes (adipose MS at 4 h only, blood OLINK, muscle MS);
#     - metabolites: the 450 network metabolites.
#   Upstream QC of these values is documented in the headers of steps 1 and 1b and in the README. The
#   maximum absolute values in Table 1 are the step 1 / 1b normalisation divisors.
#
# TECH STACK
#   R 4.4; data.table; tinytex (compiles the LaTeX file with pdflatex; TinyTeX must be installed:
#   install.packages("tinytex"); tinytex::install_tinytex()). LaTeX packages: booktabs, caption.
#
# INPUTS (files and columns)
#   $HACK_OUT/01_nodes_EE_raw_logFC.csv, 01_nodes_RE_raw_logFC.csv     (step 1)
#     entrez_gene, gene_symbol, then 18 columns <tissue>_<rna|prot>_<0.5h|4h|24h>; adipose_prot_0.5h and
#     adipose_prot_24h are empty by design
#   $HACK_OUT/01b_metab_nodes_EE_raw_logFC.csv, 01b_metab_nodes_RE_raw_logFC.csv   (step 1b)
#     metabolite, then 9 columns <tissue>_metab_<0.5h|4h|24h>
#
# OUTPUTS (files, locations, columns)
#   $HACK_FIG/13_logfc_descriptive_stats.pdf   the two tables as a PDF (default ~/Desktop/output/hackathon)
#   $HACK_FIG/13_logfc_descriptive_stats.tex   its LaTeX source
#   $HACK_OUT/13_logfc_descriptive_stats.csv   table (1 or 2), ome, arm, n, min, max, mean, sd
#   Nothing is written inside the repository.
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   n: RNA 8,478 (471 genes x 3 tissues x 3 times x 2 arms), protein 6,594 (471 x 7 observed dims x 2),
#   metabolites 8,100 (450 x 9 x 2). Pooled max |logFC|: RNA 2.746, protein 1.763, metabolites 4.204.
#   The script stops with an error if the counts differ from features x dimensions x arms, or if the
#   pooled maximum |logFC| of an ome differs from the divisor recorded by step 1 / 1b.
#
# KNOWN LIMITS
#   Values are pooled over tissues and timepoints within an ome; tissues can differ (e.g. blood OLINK vs
#   muscle MS protein). Mean and SD describe the spread of estimates, which includes measurement noise; they
#   are not tests of exercise effects.
# =====================================================================================================

# Load data.table quietly.
suppressMessages(library(data.table))

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where the PDF goes (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Readable names for the arms and omes.
ARM_NAME <- c(EE = "Endurance vs control", RE = "Resistance vs control")
# (and for the three omes)
OME_NAME <- c(rna = "RNA (RNA-seq)", prot = "Protein (MS; blood OLINK)", metab = "Metabolites")

# ---- gather every log fold change in one long table ---------------------------------------------
# Helper: read one raw-logFC table and reshape it to one row per value (feature, dimension, value).
long <- function(file, id_cols, arm) {
  x <- fread(file.path(OUT, file))
  # the dimension columns (everything after the identifier columns)
  dims <- names(x)[-(seq_len(id_cols))]
  # one row per feature x dimension; empty cells (adipose protein 0.5 / 24 h) are dropped
  v <- melt(x, id.vars = names(x)[seq_len(id_cols)], measure.vars = dims, variable.name = "dim",
            value.name = "logFC", na.rm = TRUE)
  # the ome of each dimension, from its name ("muscle_rna_4h" -> "rna", "blood_metab_4h" -> "metab")
  v[, ome := sub("^[a-z]+_([a-z]+)_.*$", "\\1", as.character(dim))]
  # label the arm
  v[, arm := arm]
  v[, .(ome, arm, logFC)]
}
# All values: genes (two identifier columns) and metabolites (one), both arms.
d <- rbind(long("01_nodes_EE_raw_logFC.csv", 2, "EE"), long("01_nodes_RE_raw_logFC.csv", 2, "RE"),
           long("01b_metab_nodes_EE_raw_logFC.csv", 1, "EE"), long("01b_metab_nodes_RE_raw_logFC.csv", 1, "RE"))
# Safety check: exactly the three omes.
stopifnot(setequal(unique(d$ome), names(OME_NAME)))

# ---- the statistics -------------------------------------------------------------------------------
# Helper: n, min, max, mean and SD of a set of values.
stats <- function(v) list(n = length(v), min = min(v), max = max(v), mean = mean(v), sd = sd(v))
# Table 1: pooled across the two arms, per ome.
t1 <- d[, stats(logFC), by = ome][, arm := "Both arms (pooled)"]
# Table 2: per ome and arm.
t2 <- d[, stats(logFC), by = .(ome, arm)][, arm := ARM_NAME[arm]]
# Validation 1: each ome has features x observed dimensions x 2 arms values (no value lost or duplicated).
expected_n <- c(rna = 471 * 9 * 2, prot = 471 * 7 * 2, metab = 450 * 9 * 2)
# (stop if any ome's count differs)
stopifnot(all(t1$n == expected_n[t1$ome]))
# Validation 2: the pooled maximum |logFC| of each ome equals the normalisation divisor of step 1 / 1b.
div <- rbind(fread(file.path(OUT, "01_scale_factors.csv"))[, .(ome, max_abs_logFC)],
             fread(file.path(OUT, "01b_metab_scale_factors.csv"))[, .(ome, max_abs_logFC)])
# (the pooled maximum |logFC| of each ome)
pooled_max <- d[, .(m = max(abs(logFC))), by = ome]
# (stop if it differs from the recorded divisor)
stopifnot(isTRUE(all.equal(pooled_max$m, div$max_abs_logFC[match(pooled_max$ome, div$ome)])))
# Fixed row order: RNA, protein, metabolites (and endurance before resistance).
ord <- function(x) x[order(match(ome, names(OME_NAME)), arm)]
# (apply that order to both tables)
t1 <- ord(t1); t2 <- ord(t2)
# Save both as one CSV (a column says which table each row belongs to).
fwrite(rbind(t1[, table := "1: across arms"], t2[, table := "2: by arm"])[, .(table, ome, arm, n, min, max, mean, sd)],
       file.path(OUT, "13_logfc_descriptive_stats.csv"))

# ---- the LaTeX document ---------------------------------------------------------------------------
# Helper: format a number with 3 decimals in LaTeX math mode (so the minus sign is typeset properly);
# values that round to zero are shown as 0.000 rather than -0.000.
f3 <- function(x) { x <- round(x, 3); x[x == 0] <- 0; sprintf("$%s$", formatC(x, format = "f", digits = 3)) }
# Helper: format a count with thousands separators.
fn <- function(x) formatC(x, format = "d", big.mark = ",")
# Rows of Table 1: ome, n, min, max, mean, SD.
rows1 <- t1[, sprintf("%s & %s & %s & %s & %s & %s \\\\", OME_NAME[ome], fn(n), f3(min), f3(max), f3(mean), f3(sd))]
# Rows of Table 2: ome (only on the first row of each ome), arm, n, min, max, mean, SD.
t2[, first := !duplicated(ome)]
# (one LaTeX row per ome x arm; a small gap after each ome's pair of rows)
rows2 <- t2[, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\%s", fifelse(first, OME_NAME[ome], ""), arm, fn(n),
                      f3(min), f3(max), f3(mean), f3(sd), fifelse(!first & ome != "metab", " \\addlinespace", ""))]
# The whole document, line by line.
tex <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=2cm]{geometry}",
  "\\usepackage{booktabs}",
  "\\usepackage[font=small,labelfont=bf]{caption}",
  "\\usepackage[T1]{fontenc}",
  "\\begin{document}",
  "\\section*{Descriptive statistics of log fold changes}",
  paste0("Stanford Multi-omics Hackathon 2026, Track 1 --- endurance vs resistance networks. Generated ",
         format(Sys.Date(), "%Y-%m-%d"), "."),
  "",
  paste("Values are the unnormalised log$_2$ fold changes of the network features: the delta-delta contrast",
        "(exercise arm's change from pre-exercise minus the control group's change) at 0.5, 4 and 24\\,h after",
        "exercise, in adipose, blood and muscle, for the 471 network genes (RNA and protein) and the 450 network",
        "metabolites (MoTrPAC human pre-suspension study). Adipose protein exists at 4\\,h only.",
        "These are the values that are normalised (each ome divided by its maximum absolute log fold change)",
        "before the network edge weights are computed."),
  "",
  "\\begin{table}[h]",
  "\\centering",
  "\\caption{Log fold changes per ome, pooled across the two study arms (endurance vs control and resistance vs control).}",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Ome & $n$ & Min & Max & Mean & SD \\\\",
  "\\midrule",
  rows1,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\begin{table}[h]",
  "\\centering",
  "\\caption{Log fold changes per ome, by study arm.}",
  "\\begin{tabular}{llrrrrr}",
  "\\toprule",
  "Ome & Study arm & $n$ & Min & Max & Mean & SD \\\\",
  "\\midrule",
  rows2,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  paste("\\noindent\\small $n$ = number of log fold changes (features $\\times$ tissues $\\times$ timepoints, and arms",
        "where pooled). SD = standard deviation. Source: \\texttt{network/13\\_logfc\\_descriptive\\_stats.R}."),
  "\\end{document}")
# Write the .tex file next to where the PDF will go.
tex_file <- file.path(FIG, "13_logfc_descriptive_stats.tex")
# (write it)
writeLines(tex, tex_file)
# Compile it to PDF with pdflatex (TinyTeX installs any missing LaTeX package automatically).
old <- setwd(FIG); on.exit(setwd(old), add = TRUE)
# (compile; returns the PDF file name)
pdf <- tinytex::pdflatex(basename(tex_file))
# Report where it was written, and show the numbers.
message("-> ", file.path(FIG, pdf))
# (show both tables on screen)
print(rbind(t1, t2[, !"first"], fill = TRUE))
