#!/usr/bin/env Rscript
# 03_edge_weights.R — weight each STRING-gated edge by the dot product of its endpoints' embeddings,
# separately per arm, giving one weighted network for EE and one for RE.
#
#   w_EE(u,v) = z_u^EE . z_v^EE      w_RE(u,v) = z_u^RE . z_v^RE
#
# The dot product runs over the WHOLE embedding: all tissues x RNA + protein x every timepoint that
# exists. The adipose protein 0.5 h and 24 h dimensions do not exist for any gene (NA in every node),
# so they are left out; that leaves 16 dimensions, identical for every edge and both arms.
# Weights are signed (positive = endpoints respond in the same direction) and untransformed.
# For methods that need positive weights, sig_<arm> = sigmoid(w) = 1 / (1 + exp(-w)) maps them to 0-1
# (negative -> below 0.5), the decoder used by DeepWalk / node2vec.
# Topology is unchanged: the edge set is exactly 02_edges.csv.
#
# Inputs:  $HACK_OUT/01_nodes_{EE,RE}.csv, $HACK_OUT/02_edges.csv
# Output:  $HACK_OUT/03_weighted_edges.csv

suppressMessages(library(data.table))

OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))

edges <- fread(file.path(OUT, "02_edges.csv"), colClasses = list(character = c("entrez_a", "entrez_b")))

embedding <- function(arm) {
  x <- fread(file.path(OUT, sprintf("01_nodes_%s.csv", arm)), colClasses = list(character = "entrez_gene"))
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$entrez_gene
  absent <- colSums(is.na(m)) == nrow(m)                 # dimensions no gene has
  stopifnot(!anyNA(m[, !absent]))                        # anything else missing would be a bug
  m[, !absent, drop = FALSE]
}
Z <- list(EE = embedding("EE"), RE = embedding("RE"))
stopifnot(identical(colnames(Z$EE), colnames(Z$RE)), identical(rownames(Z$EE), rownames(Z$RE)))
message(sprintf("dot product over %d dims: %s", ncol(Z$EE), paste(colnames(Z$EE), collapse = ", ")))

for (arm in names(Z)) {
  M <- Z[[arm]]
  w <- rowSums(M[edges$entrez_a, ] * M[edges$entrez_b, ])
  set(edges, j = paste0("w_", arm), value = w)
  set(edges, j = paste0("sig_", arm), value = plogis(w))
}
edges[, w_diff := w_EE - w_RE]

out <- edges[, .(entrez_a, symbol_a, entrez_b, symbol_b, combined_score, w_EE, w_RE, w_diff, sig_EE, sig_RE)]
setorder(out, -combined_score, symbol_a, symbol_b)
fwrite(out, file.path(OUT, "03_weighted_edges.csv"))

for (arm in names(Z)) {
  w <- out[[paste0("w_", arm)]]
  s <- out[[paste0("sig_", arm)]]
  message(sprintf("%s: %d edges, median %.2f, range %.1f..%.1f, %d negative; sigmoid > 0.99 for %d, < 0.01 for %d",
                  arm, length(w), median(w), min(w), max(w), sum(w < 0), sum(s > 0.99), sum(s < 0.01)))
}
message(sprintf("cor(w_EE, w_RE) = %.2f; %d edges change sign between arms",
                cor(out$w_EE, out$w_RE), sum(sign(out$w_EE) != sign(out$w_RE))))
message("-> ", file.path(OUT, "03_weighted_edges.csv"))
