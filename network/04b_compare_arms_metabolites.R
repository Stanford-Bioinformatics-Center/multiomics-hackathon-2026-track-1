#!/usr/bin/env Rscript
# =====================================================================================================
# 04b_compare_arms_metabolites.R — STEP 4b: WHAT DIFFERS BETWEEN THE ENDURANCE AND RESISTANCE
#                                  METABOLITE NETWORKS, BEYOND MEASUREMENT NOISE?
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The metabolite counterpart of step 4, using exactly the same method. The two metabolite networks
#   (step 6) have the same edges (shared Rhea enzyme among our 471 genes + same RefMet super class) and
#   different weights (dot products of each metabolite pair's 9-number response vectors, per arm). This
#   script compares them:
#   - per EDGE: w_diff = w_EE - w_RE;
#   - per METABOLITE: its "strength" = sum of the weights of its edges, per arm, and the difference;
#   - connected groups of edges that differ ("differential subnetworks").
#   For each difference it asks whether it is larger than the measurement noise.
#
# WHERE THE ERROR BARS COME FROM: Monte Carlo error propagation (a parametric bootstrap), as in step 4
#   1. In each of B = 10,000 rounds, every metabolite's endurance and resistance value in every dimension is
#      redrawn from a normal distribution centred on its estimate, with its standard error (step 1b).
#   2. A metabolite's endurance and resistance values are drawn TOGETHER with the correlation of their
#      errors (median ~0.65), because both arms are compared with the same control group. Ignoring it would
#      overstate the noise in each difference and hide real differences.
#   3. Every edge weight, difference and metabolite strength is recomputed in that round.
#   p = 2 x min(share of rounds <= 0, share of rounds >= 0) (two-sided; smallest possible 2/10,001);
#   Benjamini-Hochberg false discovery rate (FDR) within edges and, separately, within metabolites; the 95%
#   interval is the middle 95% of the rounds (Efron & Tibshirani 1993).
#
# ASSUMPTIONS AND LIMITS (same as step 4)
#   - Errors are treated as independent between metabolites and between the 9 dimensions of one
#     metabolite (in reality they share people and control groups; not modelled, so p-values may be
#     somewhat too small).
#   - The redraws centre on the observed estimates and add noise to them, so the between-arm correlation
#     across redraws is pulled down; it is reported as "correlation under re-measurement noise", not as a
#     confidence interval.
#   - With 122 edges, about 6 would pass an uncorrected p < 0.05 by chance alone.
#   - "Not significant" is not evidence that the arms are the same.
#
# TECH STACK
#   R 4.4; data.table, Matrix (sparse metabolite x edge matrix), igraph (connected components).
#
# INPUTS:  $HACK_OUT/01b_metab_nodes_{EE,RE}.csv, 01b_metab_nodes_{EE,RE}_se.csv, 01b_metab_nodes_arm_corr.csv,
#          $HACK_OUT/06_metabolite_edges.csv
# OUTPUTS: $HACK_OUT/04b_metab_edge_diff.csv         per edge: difference, 95% interval, p, FDR
#          $HACK_OUT/04b_metab_node_diff.csv         per metabolite with >= 1 edge: strength per arm, difference, p, FDR
#          $HACK_OUT/04b_metab_diff_subnetworks.csv  connected groups of differing edges (FDR < 0.1; exploratory p < 0.05)
#          $HACK_OUT/04b_metab_summary.csv           headline numbers
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph) })

# Results folder (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Number of bootstrap rounds (override with N_BOOT).
B <- as.integer(Sys.getenv("N_BOOT", unset = "10000"))
# The false-discovery-rate threshold for calling a difference significant.
FDR <- 0.1
# Fix the random-number seed so every run gives exactly the same numbers.
set.seed(20260926)

# Read the metabolite edges and their weights from step 6.
edges <- fread(file.path(OUT, "06_metabolite_edges.csv"))
# Helper: read a step-1b table as a matrix (one row per metabolite, 9 columns).
read_mat <- function(file) {
  x <- fread(file.path(OUT, file))
  m <- as.matrix(x[, -1]); rownames(m) <- x$metabolite
  m
}
# The vectors (ZE endurance, ZR resistance).
ZE <- read_mat("01b_metab_nodes_EE.csv"); ZR <- read_mat("01b_metab_nodes_RE.csv")
# Their standard errors (SE endurance, SR resistance) and the correlation of their errors (RHO).
SE <- read_mat("01b_metab_nodes_EE_se.csv"); SR <- read_mat("01b_metab_nodes_RE_se.csv"); RHO <- read_mat("01b_metab_nodes_arm_corr.csv")
# Safety check: all five tables line up (same metabolites, same dimensions) and have no gaps.
stopifnot(identical(dimnames(ZE), dimnames(ZR)), identical(dimnames(ZE), dimnames(SE)), identical(dimnames(ZE), dimnames(SR)),
          identical(dimnames(ZE), dimnames(RHO)), !anyNA(ZE), !anyNA(ZR), !anyNA(SE), !anyNA(SR), !anyNA(RHO))
# For every edge, the row numbers of its two metabolites (ia = first, ib = second).
mets <- rownames(ZE); ia <- match(edges$metabolite_a, mets); ib <- match(edges$metabolite_b, mets)
# Safety check: every edge's metabolites were found.
stopifnot(!anyNA(ia), !anyNA(ib))

# Helper: every edge's weight (dot product) in each arm, from a given pair of vector tables.
weights <- function(E, R) list(EE = unname(rowSums(E[ia, , drop = FALSE] * E[ib, , drop = FALSE])),
                               RE = unname(rowSums(R[ia, , drop = FALSE] * R[ib, , drop = FALSE])))
# The observed weights (from the measured values).
obs <- weights(ZE, ZR)
# Safety check: they match what step 6 wrote.
stopifnot(all.equal(obs$EE, edges$w_EE), all.equal(obs$RE, edges$w_RE))
# The same fixed 0-1 transform as step 6 (s computed once from the observed weights, kept fixed).
SIG_SCALE <- median(abs(c(obs$EE, obs$RE)))
# The 0-1 transform itself: sigmoid of (weight / s).
sig <- function(w) plogis(w / SIG_SCALE)
# Safety check: matches step 6's 0-1 values in both arms.
stopifnot(all.equal(sig(obs$EE), edges$sig_EE), all.equal(sig(obs$RE), edges$sig_RE))

# A metabolite-by-edge table of 1s: row = metabolite, column = edge, 1 if the metabolite is on the edge.
inc <- Matrix::sparseMatrix(i = c(ia, ib), j = rep(seq_along(ia), 2), x = 1, dims = c(length(mets), length(ia)))
# Helper: each metabolite's strength = the sum of the weights of its edges.
strength <- function(w) as.numeric(inc %*% w)

# Observed per-metabolite differences in strength (EE - RE), on raw and on 0-1 weights.
obs_ndiff <- strength(obs$EE) - strength(obs$RE)
# The same on the 0-1 weights.
obs_ndiff_sig <- strength(sig(obs$EE)) - strength(sig(obs$RE))
# Observed correlation between the two arms' edge weights.
obs_cor <- cor(obs$EE, obs$RE)

# ---- the bootstrap --------------------------------------------------------------------------------
# Sizes: number of metabolites and dimensions.
n_m <- nrow(ZE); n_d <- ncol(ZE)
# Storage for every round: edge differences, metabolite differences (raw and 0-1), between-arm correlation.
be <- matrix(0, B, length(ia)); bn <- matrix(0, B, n_m); bns <- matrix(0, B, n_m); bcor <- numeric(B)
# Repeat B times (one round per pass).
for (b in seq_len(B)) {
  # two independent sets of standard-normal random numbers, one per metabolite x dimension
  u1 <- matrix(rnorm(n_m * n_d), n_m); u2 <- matrix(rnorm(n_m * n_d), n_m)
  # redrawn endurance values: estimate + standard error x noise
  E <- ZE + SE * u1
  # redrawn resistance values, with noise correlated to the endurance noise by RHO
  R <- ZR + SR * (RHO * u1 + sqrt(1 - RHO^2) * u2)
  # recompute every edge weight in this round
  w <- weights(E, R)
  # store this round's edge differences, metabolite differences and correlation
  be[b, ]  <- w$EE - w$RE
  bn[b, ]  <- strength(w$EE) - strength(w$RE)
  bns[b, ] <- strength(sig(w$EE)) - strength(sig(w$RE))
  bcor[b]  <- cor(w$EE, w$RE)
}
# Two-sided bootstrap p-value per column (never below 2/(B+1)).
boot_p <- function(M) pmax(2 / (B + 1), pmin(1, 2 * pmin(colMeans(M <= 0), colMeans(M >= 0))))
# The q-th percentile of each column (for the 95% interval).
ci <- function(M, q) apply(M, 2, quantile, q, names = FALSE)

# ---- edges -----------------------------------------------------------------------------------------
# 95% interval, p-value and FDR-adjusted p for each edge difference.
edges[, `:=`(w_diff_lo = ci(be, 0.025), w_diff_hi = ci(be, 0.975), p_boot = boot_p(be))][, fdr := p.adjust(p_boot, "BH")]
# Which arm has the HIGHER weight (signed: for negative weights, higher = closer to zero)...
edges[, higher_in := fifelse(w_diff > 0, "EE", "RE")]
# ...and which arm has the LARGER weight in size, ignoring sign (the stronger co-response).
edges[, larger_abs_in := fifelse(abs(w_EE) > abs(w_RE), "EE", "RE")]
# Whether the edge is positive in one arm and negative in the other.
edges[, sign_change := sign(w_EE) != sign(w_RE)]
# Most significant first (ties: largest difference first).
edges <- edges[order(p_boot, -abs(w_diff))]
# Save the per-edge table.
fwrite(edges, file.path(OUT, "04b_metab_edge_diff.csv"))

# ---- metabolites -----------------------------------------------------------------------------------
# Per metabolite: number of edges, strength per arm (raw and 0-1), differences, 95% interval, p-values.
nodes <- data.table(metabolite = mets, degree = as.integer(Matrix::rowSums(inc)),
                    strength_EE = strength(obs$EE), strength_RE = strength(obs$RE), delta_strength = obs_ndiff,
                    sig_strength_EE = strength(sig(obs$EE)), sig_strength_RE = strength(sig(obs$RE)),
                    delta_sig_strength = obs_ndiff_sig,
                    delta_lo = ci(bn, 0.025), delta_hi = ci(bn, 0.975),
                    p_boot = boot_p(bn), p_boot_sig = boot_p(bns))
# Only metabolites with at least one edge can differ in strength.
nodes <- nodes[degree > 0]
# FDR adjustment across metabolites.
nodes[, `:=`(fdr = p.adjust(p_boot, "BH"), fdr_sig = p.adjust(p_boot_sig, "BH"))]
# Most significant first (ties: largest difference first).
nodes <- nodes[order(p_boot, -abs(delta_strength))]
# Save the per-metabolite table.
fwrite(nodes, file.path(OUT, "04b_metab_node_diff.csv"))

# ---- differential subnetworks ---------------------------------------------------------------------
# Helper: connected groups among a set of edges, described per group.
modules <- function(sub, tier) {
  # nothing to do if the set is empty
  if (!nrow(sub)) return(NULL)
  # a small network from these edges; label each edge with its connected group
  g <- graph_from_data_frame(sub[, .(metabolite_a, metabolite_b)], directed = FALSE)
  sub <- copy(sub)[, module := components(g)$membership[metabolite_a]]
  # per group: size, edges higher in each arm, total difference, metabolites and edges
  m <- sub[, .(tier = tier, n_metabolites = uniqueN(c(metabolite_a, metabolite_b)), n_edges = .N,
               n_EE_higher = sum(higher_in == "EE"), n_RE_higher = sum(higher_in == "RE"),
               sum_w_diff = sum(w_diff),
               metabolites = paste(sort(unique(c(metabolite_a, metabolite_b))), collapse = ";"),
               edges = paste(sprintf("%s-%s(%s higher)", metabolite_a, metabolite_b, higher_in), collapse = ";")), by = module]
  # largest groups first, renumbered
  m[order(-n_edges, -abs(sum_w_diff))][, module := seq_len(.N)]
}
# Edges that pass the FDR threshold.
sig_edges <- edges[fdr < FDR]
# Groups among FDR-significant edges (primary) and among uncorrected p < 0.05 edges (exploratory).
mods <- rbind(modules(sig_edges, "fdr_lt_0.1"), modules(edges[p_boot < 0.05], "exploratory_p_lt_0.05"))
# An empty table with the right columns if there are no groups.
if (is.null(mods)) mods <- data.table(tier = character(), module = integer(), n_metabolites = integer(), n_edges = integer(),
                                      n_EE_higher = integer(), n_RE_higher = integer(), sum_w_diff = numeric(),
                                      metabolites = character(), edges = character())
# Save the subnetwork table.
fwrite(mods, file.path(OUT, "04b_metab_diff_subnetworks.csv"))

# ---- headline numbers ------------------------------------------------------------------------------
summ <- data.table(
  metric = c("edges", "bootstrap_draws", "cor_w_EE_w_RE", "cor_remeasured_2.5pct", "cor_remeasured_median", "cor_remeasured_97.5pct",
             "edges_sign_change", "edges_fdr_lt_0.1", "edges_fdr_EE_higher", "edges_fdr_RE_higher", "edges_p_lt_0.05",
             "edges_p_lt_0.05_expected_by_chance", "metabolites_with_edges", "metabolites_fdr_lt_0.1_signed",
             "metabolites_fdr_lt_0.1_sigmoid", "diff_subnetworks_fdr", "diff_subnetworks_exploratory"),
  value = c(nrow(edges), B, obs_cor, quantile(bcor, 0.025), median(bcor), quantile(bcor, 0.975),
            sum(edges$sign_change), nrow(sig_edges), sum(sig_edges$higher_in == "EE"), sum(sig_edges$higher_in == "RE"),
            sum(edges$p_boot < 0.05), 0.05 * nrow(edges), nrow(nodes), sum(nodes$fdr < FDR), sum(nodes$fdr_sig < FDR),
            sum(mods$tier == "fdr_lt_0.1"), sum(mods$tier != "fdr_lt_0.1")))
# Save and show.
fwrite(summ, file.path(OUT, "04b_metab_summary.csv"))
# Show them on screen.
print(summ)
