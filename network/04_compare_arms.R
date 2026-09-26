#!/usr/bin/env Rscript
# =====================================================================================================
# 04_compare_arms.R — STEP 4: ARE THE ENDURANCE AND RESISTANCE NETWORKS MORE ALIKE OR MORE DIFFERENT?
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The two networks have the same edges (from STRING) but different weights (from each arm's exercise
#   response, step 3). This script compares them honestly, without assuming they are alike or different:
#   - per EDGE: w_diff = w_EE - w_RE, the difference in strength between the arms;
#   - per GENE: its "strength" = the sum of the weights of all its edges, in each arm, and the difference;
#   - small groups of connected edges that all differ ("differential subnetworks").
#   For each difference it asks: is it bigger than the measurement noise? That needs error bars.
#
# WHERE THE ERROR BARS COME FROM: a parametric bootstrap
#   Every number in a gene's vector is an estimate with a standard error (from the consortium's models,
#   saved in step 1). The bootstrap asks "if the measurements had come out slightly differently, as their
#   error bars allow, how much would each difference move?":
#     1. In each of B = 10,000 rounds, re-draw every gene's endurance and resistance value in every
#        dimension from a normal distribution centred on the estimate, with its standard error.
#     2. The endurance and resistance values of a gene are drawn TOGETHER with their correlation
#        (median ~0.5-0.65), because both arms are compared against the same control group; ignoring
#        this would make the arms look more different than they are.
#     3. Recompute every edge weight, every difference and every gene strength in that round.
#   A difference is called real if, across the 10,000 rounds, it rarely crosses zero:
#     p = 2 x min(share of rounds <= 0, share of rounds >= 0)       (two-sided; smallest possible 2/10,001)
#   p-values are corrected for testing many edges (and, separately, many genes) with the
#   Benjamini-Hochberg false discovery rate (FDR). The 95% interval is the middle 95% of the rounds.
#   This is a standard parametric bootstrap (Efron & Tibshirani 1993, An Introduction to the Bootstrap).
#
# ASSUMPTIONS AND LIMITS
#   - Errors are treated as independent between genes (the consortium fitted each gene separately). In
#     reality genes measured in the same people are somewhat correlated; this is not modelled.
#   - The bootstrap centres on the observed estimates, so it gives error bars around what was measured.
#   - Most responses are small relative to their error (median |value| / standard error = 0.78), so the
#     test has limited power: "no significant difference" is NOT evidence that the arms are the same.
#
# A TEST WE TRIED AND REJECTED (kept here so nobody repeats it)
#   Shuffling which arm each gene's vector belongs to ("arm-label swap") cannot work per edge: an edge has
#   only two genes, and swapping both just flips the sign of w_diff, so half the shuffles reproduce the
#   observed difference exactly and no edge can ever reach p < ~0.5.
#
# TECH STACK
#   R 4.4; data.table (tables), Matrix (sparse gene-by-edge matrix for fast sums), igraph (components).
#
# INPUTS:  $HACK_OUT/01_nodes_{EE,RE}.csv, 01_nodes_{EE,RE}_se.csv, 01_nodes_arm_corr.csv (step 1),
#          $HACK_OUT/03_weighted_edges.csv (step 3)
# OUTPUTS: $HACK_OUT/04_edge_diff.csv         per edge: difference, 95% interval, p, FDR
#          $HACK_OUT/04_node_diff.csv         per gene with >= 1 edge: strength per arm, difference, p, FDR
#          $HACK_OUT/04_diff_subnetworks.csv  connected groups of differing edges (FDR < 0.1, and
#                                             an exploratory tier at uncorrected p < 0.05)
#          $HACK_OUT/04_summary.csv           headline numbers
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

# Read the weighted edges from step 3.
edges <- fread(file.path(OUT, "03_weighted_edges.csv"), colClasses = list(character = c("entrez_a", "entrez_b")))
# Helper: read a step-1 table as a matrix (one row per gene), dropping the columns empty for every gene.
read_mat <- function(file) {
  x <- fread(file.path(OUT, file), colClasses = list(character = "entrez_gene"))
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$entrez_gene
  m[, colSums(is.na(m)) < nrow(m), drop = FALSE]          # same 16 dims as 03
}
# The vectors (ZE endurance, ZR resistance), their standard errors (SE, SR) and the EE-RE correlation (RHO).
ZE <- read_mat("01_nodes_EE.csv"); ZR <- read_mat("01_nodes_RE.csv")
SE <- read_mat("01_nodes_EE_se.csv"); SR <- read_mat("01_nodes_RE_se.csv"); RHO <- read_mat("01_nodes_arm_corr.csv")
# Safety check: all five tables line up (same genes, same dimensions) and have no gaps.
stopifnot(identical(dimnames(ZE), dimnames(SE)), identical(dimnames(ZE), SR |> dimnames()),
          identical(dimnames(ZE), dimnames(RHO)), !anyNA(SE), !anyNA(SR), !anyNA(RHO))
# For every edge, the row numbers of its two genes (ia = first gene, ib = second gene).
genes <- rownames(ZE); ia <- match(edges$entrez_a, genes); ib <- match(edges$entrez_b, genes)
# Safety check: every edge's genes were found.
stopifnot(!anyNA(ia), !anyNA(ib))

# Helper: every edge's weight (dot product) in each arm, from a given pair of vector tables.
weights <- function(E, R) list(EE = unname(rowSums(E[ia, ] * E[ib, ])), RE = unname(rowSums(R[ia, ] * R[ib, ])))
# The observed weights (from the measured values).
obs <- weights(ZE, ZR)
# Safety check: they match what step 3 wrote.
stopifnot(all.equal(obs$EE, edges$w_EE), all.equal(obs$RE, edges$w_RE))   # matches 03
# The same fixed 0-1 transform as step 3 (the scale s is computed once, from the observed weights,
# and kept fixed in every bootstrap round, because it is a unit, not something being estimated).
SIG_SCALE <- median(abs(c(obs$EE, obs$RE)))          # fixed transform, same constant as 03
sig <- function(w) plogis(w / SIG_SCALE)
# Safety check: matches step 3's 0-1 values.
stopifnot(all.equal(sig(obs$EE), edges$sig_EE))

# A gene-by-edge table of 1s ("incidence matrix"): row = gene, column = edge, 1 if the gene is on the edge.
# Multiplying it by the edge weights adds up each gene's edge weights in one step (its "strength").
inc <- Matrix::sparseMatrix(i = c(ia, ib), j = rep(seq_along(ia), 2), x = 1, dims = c(length(genes), length(ia)))
strength <- function(w) as.numeric(inc %*% w)

# Observed per-gene differences in strength (EE - RE), on raw weights and on 0-1 weights.
obs_ndiff <- strength(obs$EE) - strength(obs$RE)
obs_ndiff_sig <- strength(sig(obs$EE)) - strength(sig(obs$RE))
# Observed correlation between the two arms' edge weights (1 = identical pattern, 0 = unrelated).
obs_cor <- cor(obs$EE, obs$RE)

# ---- the bootstrap --------------------------------------------------------------------------------
# Sizes: number of genes and dimensions.
n_g <- nrow(ZE); n_d <- ncol(ZE)
# Storage for every round: edge differences (be), gene differences on raw (bn) and 0-1 (bns) weights,
# and the between-arm correlation (bcor).
be <- matrix(0, B, length(ia)); bn <- matrix(0, B, n_g); bns <- matrix(0, B, n_g); bcor <- numeric(B)
for (b in seq_len(B)) {
  # two independent sets of standard-normal random numbers, one per gene x dimension
  u1 <- matrix(rnorm(n_g * n_d), n_g); u2 <- matrix(rnorm(n_g * n_d), n_g)
  # re-drawn endurance values: estimate + standard error x noise
  E <- ZE + SE * u1
  # re-drawn resistance values, with noise correlated to the endurance noise by RHO
  R <- ZR + SR * (RHO * u1 + sqrt(1 - RHO^2) * u2)
  # recompute every edge weight in this round
  w <- weights(E, R)
  # store this round's edge differences, gene differences and correlation
  be[b, ]  <- w$EE - w$RE
  bn[b, ]  <- strength(w$EE) - strength(w$RE)
  bns[b, ] <- strength(sig(w$EE)) - strength(sig(w$RE))
  bcor[b]  <- cor(w$EE, w$RE)
}
# Two-sided bootstrap p-value per column: how often the difference crosses zero (never below 2/(B+1)).
boot_p <- function(M) pmax(2 / (B + 1), pmin(1, 2 * pmin(colMeans(M <= 0), colMeans(M >= 0))))
# The q-th percentile of each column (used for the 95% interval).
ci <- function(M, q) apply(M, 2, quantile, q, names = FALSE)

# ---- edges -----------------------------------------------------------------------------------------
# 95% interval, p-value and false-discovery-rate-adjusted p (FDR) for each edge difference.
edges[, `:=`(w_diff_lo = ci(be, 0.025), w_diff_hi = ci(be, 0.975), p_boot = boot_p(be))][, fdr := p.adjust(p_boot, "BH")]
# Which arm the edge is stronger in.
edges[, stronger_in := fifelse(w_diff > 0, "EE", "RE")]
# Whether the edge is positive in one arm and negative in the other.
edges[, sign_change := sign(w_EE) != sign(w_RE)]
# Most significant first, then save.
edges <- edges[order(p_boot, -abs(w_diff))]
fwrite(edges, file.path(OUT, "04_edge_diff.csv"))

# ---- genes -----------------------------------------------------------------------------------------
# Gene ID -> symbol lookup, from the edge list.
sym <- unique(rbind(edges[, .(entrez = entrez_a, symbol = symbol_a)], edges[, .(entrez = entrez_b, symbol = symbol_b)]))
# Per gene: number of edges, strength in each arm (raw and 0-1), differences, 95% interval and p-values.
nodes <- data.table(entrez_gene = genes, degree = as.integer(Matrix::rowSums(inc)),
                    strength_EE = strength(obs$EE), strength_RE = strength(obs$RE), delta_strength = obs_ndiff,
                    sig_strength_EE = strength(sig(obs$EE)), sig_strength_RE = strength(sig(obs$RE)),
                    delta_sig_strength = obs_ndiff_sig,
                    delta_lo = ci(bn, 0.025), delta_hi = ci(bn, 0.975),
                    p_boot = boot_p(bn), p_boot_sig = boot_p(bns))
# Only genes with at least one edge can differ in strength.
nodes <- nodes[degree > 0]
# False-discovery-rate adjustment across genes.
nodes[, `:=`(fdr = p.adjust(p_boot, "BH"), fdr_sig = p.adjust(p_boot_sig, "BH"))]
# Add symbols, put ID and symbol first, most significant first, then save.
nodes[, gene_symbol := sym$symbol[match(entrez_gene, sym$entrez)]]
setcolorder(nodes, c("entrez_gene", "gene_symbol"))
nodes <- nodes[order(p_boot, -abs(delta_strength))]
fwrite(nodes, file.path(OUT, "04_node_diff.csv"))

# ---- differential subnetworks: connected groups of differing edges --------------------------------
# Helper: given a set of edges, find the groups that are connected to each other and describe each group.
modules <- function(sub, tier) {
  # nothing to do if the set is empty
  if (!nrow(sub)) return(NULL)
  # build a small network from these edges and label each edge with its connected group
  g <- igraph::graph_from_data_frame(sub[, .(symbol_a, symbol_b)], directed = FALSE)
  sub <- copy(sub)[, module := igraph::components(g)$membership[symbol_a]]
  # per group: size, how many edges are stronger in each arm, total difference, and the genes/edges
  m <- sub[, .(tier = tier, n_genes = uniqueN(c(symbol_a, symbol_b)), n_edges = .N,
               n_EE_stronger = sum(stronger_in == "EE"), n_RE_stronger = sum(stronger_in == "RE"),
               sum_w_diff = sum(w_diff),
               genes = paste(sort(unique(c(symbol_a, symbol_b))), collapse = ";"),
               edges = paste(sprintf("%s-%s(%s)", symbol_a, symbol_b, stronger_in), collapse = ";")), by = module]
  # largest groups first, renumbered 1, 2, 3, ...
  m[order(-n_edges, -abs(sum_w_diff))][, module := seq_len(.N)]
}
# Edges that pass the FDR threshold.
sig_edges <- edges[fdr < FDR]
# Groups among FDR-significant edges (primary) and among uncorrected p < 0.05 edges (exploratory only).
mods <- rbind(modules(sig_edges, "fdr_lt_0.1"), modules(edges[p_boot < 0.05], "exploratory_p_lt_0.05"))
# If there are no groups at all, write an empty table with the right column names.
if (is.null(mods)) mods <- data.table(tier = character(), module = integer(), n_genes = integer(), n_edges = integer(),
                                      n_EE_stronger = integer(), n_RE_stronger = integer(), sum_w_diff = numeric(),
                                      genes = character(), edges = character())
fwrite(mods, file.path(OUT, "04_diff_subnetworks.csv"))

# ---- headline numbers ------------------------------------------------------------------------------
summ <- data.table(
  metric = c("edges", "bootstrap_draws", "cor_w_EE_w_RE", "cor_boot_lo", "cor_boot_hi",
             "edges_sign_change", "edges_fdr_lt_0.1", "edges_fdr_EE_stronger", "edges_fdr_RE_stronger", "edges_p_lt_0.05",
             "nodes_with_edges", "nodes_fdr_lt_0.1_signed", "nodes_fdr_lt_0.1_sigmoid", "diff_subnetworks_fdr", "diff_subnetworks_exploratory"),
  value = c(nrow(edges), B, obs_cor, quantile(bcor, 0.025), quantile(bcor, 0.975),
            sum(edges$sign_change), nrow(sig_edges), sum(sig_edges$stronger_in == "EE"), sum(sig_edges$stronger_in == "RE"), sum(edges$p_boot < 0.05),
            nrow(nodes), sum(nodes$fdr < FDR), sum(nodes$fdr_sig < FDR), sum(mods$tier == "fdr_lt_0.1"), sum(mods$tier != "fdr_lt_0.1")))
# Save and show.
fwrite(summ, file.path(OUT, "04_summary.csv"))
print(summ)
