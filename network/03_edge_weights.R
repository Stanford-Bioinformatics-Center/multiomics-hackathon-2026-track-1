#!/usr/bin/env Rscript
# =====================================================================================================
# 03_edge_weights.R — STEP 3: HOW STRONG IS EACH EDGE, SEPARATELY IN ENDURANCE AND IN RESISTANCE
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Step 2 decided WHICH genes are connected (the same edges for both arms). This step gives every edge
#   a strength ("weight") in each arm, from the exercise data. The weight is the dot product of the two
#   genes' response vectors from step 1: multiply the two genes' values dimension by dimension and add
#   up the 16 products. That gives two weighted networks with the same edges but different weights:
#       w_EE(u,v) = z_u(EE) . z_v(EE)        w_RE(u,v) = z_u(RE) . z_v(RE)        w_diff = w_EE - w_RE
#
# HOW TO READ A WEIGHT
#   - Large positive: on balance across the 16 dimensions, the two genes respond in the same direction
#     (both up or both down), and strongly.
#   - Large negative: on balance, they respond in opposite directions, and strongly.
#   - Near zero: either at least one gene barely responds, OR the two respond in the same direction in
#     some dimensions and opposite in others, so the products cancel out.
#   The dot product runs over the WHOLE vector at once (all tissues x RNA and protein x every time), not
#   tissue by tissue. The adipose protein 0.5 h and 24 h columns exist for no gene, so they are left out;
#   every edge uses the same 16 dimensions.
#
# METHOD AND PRECEDENT FOR THE DOT PRODUCT
#   The node-embedding "encoder-decoder" framework (Hamilton, Ying & Leskovec 2017, IEEE Data Eng.
#   Bulletin; Stanford CS224W lecture on node embeddings) scores the similarity of two nodes u, v as the
#   dot product of their vectors, z_u . z_v. Our vectors are measured exercise responses rather than
#   learned coordinates, so we use the dot product directly as the edge weight. The raw weight w is the
#   primary quantity: it is signed, untransformed, and on the same scale in both arms (step 1 divides each
#   ome by one shared maximum |logFC|, the critical normalisation that makes these dot products meaningful
#   across RNA and protein), so w_diff can be read directly.
#
# THE 0-1 VERSION (sig_EE, sig_RE) AND WHY IT IS DIVIDED BY THE MEDIAN
#   Some network methods (random walks, community detection, shortest paths) need weights that are all
#   positive. For those we provide sig = sigmoid(w / s) = 1 / (1 + exp(-w / s)), which maps any number
#   to 0-1: negative weights fall below 0.5, zero lands exactly on 0.5, positive weights go above 0.5.
#   The sigmoid of a dot product is the standard "edge probability" decoder of embedding methods such as
#   LINE (Tang et al. 2015) and graph autoencoders (Kipf & Welling 2016), and it appears in the
#   negative-sampling training of DeepWalk / node2vec.
#   Why divide by s first: those methods LEARN vectors whose dot products naturally sit in the few-units
#   range where the sigmoid is informative. Our vectors are not learned: after the step 1 normalisation
#   (values in -1..+1) their dot products are small (about -0.3 to +0.5), so a plain sigmoid(w) would
#   squeeze every weight into roughly 0.43-0.62 and hide the differences between edges. Dividing by s
#   spreads them over the 0-1 range (about 7% then end up above 0.99 or below 0.01). Dividing by s is
#   called "temperature scaling" (Hinton, Vinyals & Dean 2015; Guo et al. 2017 for calibration): it
#   changes how steep the sigmoid is, not the order of the edges or their signs.
#   Why the median: by analogy with the "median heuristic", the standard default for setting the width of
#   a similarity kernel from the median distance between data points (Gretton et al. 2012, JMLR, kernel
#   two-sample test), we set s from the median size of the weights. It is robust to the few very large
#   weights. This is a HEURISTIC, not a derived optimum: it makes the typical edge land at
#   sigmoid(+-1) = 0.27 / 0.73.
#   Why ONE s for both arms: s is a unit of measurement, like choosing centimetres over inches. Using the
#   same s for both arms means a given weight maps to the same 0-1 value in either arm, so any real
#   difference between the arms is carried through. Giving each arm its own s would do the opposite: it
#   would rescale each network to its own typical edge and could hide a genuine overall difference.
#   Caveat: the sigmoid is non-linear, so it compresses big differences. Conclusions about whether the
#   arms are more similar or different should be drawn from the raw w, with the 0-1
#   version used only where a method requires positive weights.
#
# TECH STACK
#   R 4.4; data.table.
#
# INPUTS:  $HACK_OUT/01_nodes_{EE,RE}.csv (step 1), $HACK_OUT/02_edges.csv (step 2)
# OUTPUT:  $HACK_OUT/03_weighted_edges.csv   one row per edge: w_EE, w_RE, w_diff, sig_EE, sig_RE
# =====================================================================================================

# Load data.table quietly.
suppressMessages(library(data.table))

# Results folder (override with HACK_OUT); same folder steps 1 and 2 wrote to.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))

# Read the edge list from step 2 (gene IDs kept as text).
edges <- fread(file.path(OUT, "02_edges.csv"), colClasses = list(character = c("entrez_a", "entrez_b")))

# Helper: read one arm's vectors from step 1 as a matrix (one row per gene, one column per dimension).
embedding <- function(arm) {
  # read the node table
  x <- fread(file.path(OUT, sprintf("01_nodes_%s.csv", arm)), colClasses = list(character = "entrez_gene"))
  # keep only the numeric dimensions and label the rows with gene IDs
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$entrez_gene
  # columns that are empty for every gene (adipose protein 0.5 h and 24 h)
  absent <- colSums(is.na(m)) == nrow(m)                 # dimensions no gene has
  # safety check: nothing else may be missing; if it were, that would be a bug upstream
  stopifnot(!anyNA(m[, !absent]))                        # anything else missing would be a bug
  # return the matrix without the empty columns (16 dimensions)
  m[, !absent, drop = FALSE]
}
# Read both arms.
Z <- list(EE = embedding("EE"), RE = embedding("RE"))
# Safety check: same dimensions and same genes, in the same order, in both arms.
stopifnot(identical(colnames(Z$EE), colnames(Z$RE)), identical(rownames(Z$EE), rownames(Z$RE)))
# Report which dimensions go into the dot product.
message(sprintf("dot product over %d dims: %s", ncol(Z$EE), paste(colnames(Z$EE), collapse = ", ")))

# For each arm: the weight of every edge = sum over dimensions of (gene A value x gene B value).
for (arm in names(Z)) {
  M <- Z[[arm]]
  set(edges, j = paste0("w_", arm), value = rowSums(M[edges$entrez_a, ] * M[edges$entrez_b, ]))
}
# The scale s for the 0-1 version: the median size of a weight, pooled over BOTH arms (one shared unit).
SIG_SCALE <- median(abs(c(edges$w_EE, edges$w_RE)))
# The 0-1 version of each weight: sigmoid(w / s).
edges[, `:=`(sig_EE = plogis(w_EE / SIG_SCALE), sig_RE = plogis(w_RE / SIG_SCALE))]
# Report s.
message(sprintf("sigmoid scale s = median |w| over both arms = %.3f", SIG_SCALE))
# The per-edge difference between arms (positive = the endurance weight is higher; for negative weights
# "higher" means closer to zero, so this is not the same as "stronger").
edges[, w_diff := w_EE - w_RE]

# Keep the columns we report, sort by STRING score, and save.
out <- edges[, .(entrez_a, symbol_a, entrez_b, symbol_b, combined_score, w_EE, w_RE, w_diff, sig_EE, sig_RE)]
# Sort rows: highest STRING score first, then alphabetically.
setorder(out, -combined_score, symbol_a, symbol_b)
# Save the weighted edge list.
fwrite(out, file.path(OUT, "03_weighted_edges.csv"))

# Print a short summary per arm: typical weight, range, how many are negative, and how many 0-1 values
# are pinned at the extremes (a check that the sigmoid is not saturated).
for (arm in names(Z)) {
  w <- out[[paste0("w_", arm)]]
  s <- out[[paste0("sig_", arm)]]
  message(sprintf("%s: %d edges, median %.2f, range %.1f..%.1f, %d negative; sigmoid > 0.99 for %d, < 0.01 for %d",
                  arm, length(w), median(w), min(w), max(w), sum(w < 0), sum(s > 0.99), sum(s < 0.01)))
}
# How alike the two arms' weights are overall, and how many edges flip sign between arms.
message(sprintf("cor(w_EE, w_RE) = %.2f; %d edges change sign between arms",
                cor(out$w_EE, out$w_RE), sum(sign(out$w_EE) != sign(out$w_RE))))
# Report where it was written.
message("-> ", file.path(OUT, "03_weighted_edges.csv"))
