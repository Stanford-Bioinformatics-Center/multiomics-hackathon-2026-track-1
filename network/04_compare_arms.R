#!/usr/bin/env Rscript
# 04_compare_arms.R — compare the EE and RE weighted networks (same STRING topology, per-arm weights).
#
# Edge level:  w_diff = w_EE - w_RE (signed dot products from 03).
# Node level:  strength = sum of incident edge weights, per arm, both signed (sum w) and positive
#              (sum sigmoid(w), the weight used for positive-only methods); delta = EE - RE.
# Differential subnetworks: connected components of the edges with FDR < 0.1.
#
# Uncertainty (parametric bootstrap): every embedding value is an estimate with a standard error.
# In each of B draws, each gene's (EE, RE) value in each dimension is resampled from a bivariate
# normal centred on the estimate, with the two SEs and the EE-RE correlation from step 1 (the arms
# share one control group, rho ~0.5-0.65). All weights and differences are recomputed per draw.
# Two-sided bootstrap p = 2 x min(share of draws <= 0, share >= 0), floored at 2/(B+1); BH within edges
# and within nodes; 95% percentile intervals are reported. Errors are treated as independent across
# genes (the models are fitted gene by gene), which ignores between-gene correlation.
#
# Why not an arm-label swap null: swapping a gene's EE and RE vectors only flips the sign of an
# edge's w_diff when both ends swap, so half of all draws reproduce |w_diff| exactly and no edge can
# reach p < ~0.5. It cannot detect anything at the edge level.
#
# Inputs:  $HACK_OUT/01_nodes_{EE,RE}.csv, 01_nodes_{EE,RE}_se.csv, 01_nodes_arm_corr.csv,
#          $HACK_OUT/03_weighted_edges.csv
# Outputs: $HACK_OUT/04_edge_diff.csv, 04_node_diff.csv, 04_diff_subnetworks.csv, 04_summary.csv

suppressMessages({ library(data.table); library(igraph) })

OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
B <- as.integer(Sys.getenv("N_BOOT", unset = "10000"))
FDR <- 0.1
set.seed(20260926)

edges <- fread(file.path(OUT, "03_weighted_edges.csv"), colClasses = list(character = c("entrez_a", "entrez_b")))
read_mat <- function(file) {
  x <- fread(file.path(OUT, file), colClasses = list(character = "entrez_gene"))
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$entrez_gene
  m[, colSums(is.na(m)) < nrow(m), drop = FALSE]          # same 16 dims as 03
}
ZE <- read_mat("01_nodes_EE.csv"); ZR <- read_mat("01_nodes_RE.csv")
SE <- read_mat("01_nodes_EE_se.csv"); SR <- read_mat("01_nodes_RE_se.csv"); RHO <- read_mat("01_nodes_arm_corr.csv")
stopifnot(identical(dimnames(ZE), dimnames(SE)), identical(dimnames(ZE), SR |> dimnames()),
          identical(dimnames(ZE), dimnames(RHO)), !anyNA(SE), !anyNA(SR), !anyNA(RHO))
genes <- rownames(ZE); ia <- match(edges$entrez_a, genes); ib <- match(edges$entrez_b, genes)
stopifnot(!anyNA(ia), !anyNA(ib))

weights <- function(E, R) list(EE = unname(rowSums(E[ia, ] * E[ib, ])), RE = unname(rowSums(R[ia, ] * R[ib, ])))
obs <- weights(ZE, ZR)
stopifnot(all.equal(obs$EE, edges$w_EE), all.equal(obs$RE, edges$w_RE))   # matches 03

# incidence: node x edge, so strength = inc %*% w
inc <- Matrix::sparseMatrix(i = c(ia, ib), j = rep(seq_along(ia), 2), x = 1, dims = c(length(genes), length(ia)))
strength <- function(w) as.numeric(inc %*% w)

obs_ndiff <- strength(obs$EE) - strength(obs$RE)
obs_ndiff_sig <- strength(plogis(obs$EE)) - strength(plogis(obs$RE))
obs_cor <- cor(obs$EE, obs$RE)

# bootstrap draws: (eE, eR) with sd (SE, SR) and correlation RHO, per gene x dimension
n_g <- nrow(ZE); n_d <- ncol(ZE)
be <- matrix(0, B, length(ia)); bn <- matrix(0, B, n_g); bns <- matrix(0, B, n_g); bcor <- numeric(B)
for (b in seq_len(B)) {
  u1 <- matrix(rnorm(n_g * n_d), n_g); u2 <- matrix(rnorm(n_g * n_d), n_g)
  E <- ZE + SE * u1
  R <- ZR + SR * (RHO * u1 + sqrt(1 - RHO^2) * u2)
  w <- weights(E, R)
  be[b, ]  <- w$EE - w$RE
  bn[b, ]  <- strength(w$EE) - strength(w$RE)
  bns[b, ] <- strength(plogis(w$EE)) - strength(plogis(w$RE))
  bcor[b]  <- cor(w$EE, w$RE)
}
boot_p <- function(M) pmax(2 / (B + 1), pmin(1, 2 * pmin(colMeans(M <= 0), colMeans(M >= 0))))
ci <- function(M, q) apply(M, 2, quantile, q, names = FALSE)

# ---- edges ---------------------------------------------------------------------------
edges[, `:=`(w_diff_lo = ci(be, 0.025), w_diff_hi = ci(be, 0.975), p_boot = boot_p(be))][, fdr := p.adjust(p_boot, "BH")]
edges[, stronger_in := fifelse(w_diff > 0, "EE", "RE")]
edges[, sign_change := sign(w_EE) != sign(w_RE)]
edges <- edges[order(p_boot, -abs(w_diff))]
fwrite(edges, file.path(OUT, "04_edge_diff.csv"))

# ---- nodes ---------------------------------------------------------------------------
sym <- unique(rbind(edges[, .(entrez = entrez_a, symbol = symbol_a)], edges[, .(entrez = entrez_b, symbol = symbol_b)]))
nodes <- data.table(entrez_gene = genes, degree = as.integer(Matrix::rowSums(inc)),
                    strength_EE = strength(obs$EE), strength_RE = strength(obs$RE), delta_strength = obs_ndiff,
                    sig_strength_EE = strength(plogis(obs$EE)), sig_strength_RE = strength(plogis(obs$RE)),
                    delta_sig_strength = obs_ndiff_sig,
                    delta_lo = ci(bn, 0.025), delta_hi = ci(bn, 0.975),
                    p_boot = boot_p(bn), p_boot_sig = boot_p(bns))
nodes <- nodes[degree > 0]
nodes[, `:=`(fdr = p.adjust(p_boot, "BH"), fdr_sig = p.adjust(p_boot_sig, "BH"))]
nodes[, gene_symbol := sym$symbol[match(entrez_gene, sym$entrez)]]
setcolorder(nodes, c("entrez_gene", "gene_symbol"))
nodes <- nodes[order(p_boot, -abs(delta_strength))]
fwrite(nodes, file.path(OUT, "04_node_diff.csv"))

# ---- differential subnetworks: components of FDR < 0.1 edges (primary) and p < 0.05 (exploratory)
modules <- function(sub, tier) {
  if (!nrow(sub)) return(NULL)
  g <- igraph::graph_from_data_frame(sub[, .(symbol_a, symbol_b)], directed = FALSE)
  sub <- copy(sub)[, module := igraph::components(g)$membership[symbol_a]]
  m <- sub[, .(tier = tier, n_genes = uniqueN(c(symbol_a, symbol_b)), n_edges = .N,
               n_EE_stronger = sum(stronger_in == "EE"), n_RE_stronger = sum(stronger_in == "RE"),
               sum_w_diff = sum(w_diff),
               genes = paste(sort(unique(c(symbol_a, symbol_b))), collapse = ";"),
               edges = paste(sprintf("%s-%s(%s)", symbol_a, symbol_b, stronger_in), collapse = ";")), by = module]
  m[order(-n_edges, -abs(sum_w_diff))][, module := seq_len(.N)]
}
sig <- edges[fdr < FDR]
mods <- rbind(modules(sig, "fdr_lt_0.1"), modules(edges[p_boot < 0.05], "exploratory_p_lt_0.05"))
if (is.null(mods)) mods <- data.table(tier = character(), module = integer(), n_genes = integer(), n_edges = integer(),
                                      n_EE_stronger = integer(), n_RE_stronger = integer(), sum_w_diff = numeric(),
                                      genes = character(), edges = character())
fwrite(mods, file.path(OUT, "04_diff_subnetworks.csv"))

# ---- summary -------------------------------------------------------------------------
summ <- data.table(
  metric = c("edges", "bootstrap_draws", "cor_w_EE_w_RE", "cor_boot_lo", "cor_boot_hi",
             "edges_sign_change", "edges_fdr_lt_0.1", "edges_fdr_EE_stronger", "edges_fdr_RE_stronger", "edges_p_lt_0.05",
             "nodes_with_edges", "nodes_fdr_lt_0.1_signed", "nodes_fdr_lt_0.1_sigmoid", "diff_subnetworks_fdr", "diff_subnetworks_exploratory"),
  value = c(nrow(edges), B, obs_cor, quantile(bcor, 0.025), quantile(bcor, 0.975),
            sum(edges$sign_change), nrow(sig), sum(sig$stronger_in == "EE"), sum(sig$stronger_in == "RE"), sum(edges$p_boot < 0.05),
            nrow(nodes), sum(nodes$fdr < FDR), sum(nodes$fdr_sig < FDR), sum(mods$tier == "fdr_lt_0.1"), sum(mods$tier != "fdr_lt_0.1")))
fwrite(summ, file.path(OUT, "04_summary.csv"))
print(summ)
