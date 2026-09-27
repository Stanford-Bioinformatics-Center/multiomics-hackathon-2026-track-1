# =====================================================================================================
# methods.R — reading a WeightedNetwork: edge table, node strength, arm-specific edges, printing
# =====================================================================================================

#' The edges of a network as a data.table.
#' @param net `WeightedNetwork`.
#' @return data.table with a, b, edge_type, source, evidence, w_EE, w_RE, w_diff, n_dims.
edge_table <- function(net) {
  if (!is(net, "WeightedNetwork")) exnet_abort("'net' must be a WeightedNetwork")
  data.table::as.data.table(net@edges)
}

#' Node strength per arm: the sum of |weight| over a node's edges (the weighted-degree definition of a hub).
#' @param net `WeightedNetwork`.
#' @return data.table node, degree, strength_EE, strength_RE, sorted by the larger strength.
node_strength <- function(net) {
  e <- edge_table(net)
  # Each edge counts once for each of its two ends.
  long <- data.table::rbindlist(list(e[, list(node = a, w_EE, w_RE)], e[, list(node = b, w_EE, w_RE)]))
  s <- long[, list(degree = .N, strength_EE = sum(abs(w_EE)), strength_RE = sum(abs(w_RE))), by = "node"]
  s[order(-pmax(s$strength_EE, s$strength_RE))]
}

#' Arm-specific edges: strongly co-regulated in ONE arm but not the other.
#' An edge is "endurance-specific" if w_EE >= tau and w_RE < tau, "resistance-specific" if w_RE >= tau and
#' w_EE < tau, "both" if both reach tau, otherwise "neither". Only positive weights (co-regulation in the same
#' direction) can reach tau. By default tau is the q-quantile of |w| pooled over both arms, one shared bar.
#' @param net `WeightedNetwork`.
#' @param q `numeric(1)` in (0, 1), quantile used when `tau` is not given (default 0.75).
#' @param tau `numeric(1)` or NULL, an explicit threshold (overrides `q`).
#' @return data.table of the edges with columns tau and specificity.
arm_specific_edges <- function(net, q = 0.75, tau = NULL) {
  e <- edge_table(net)
  exnet_check(checkmate::check_number(q, lower = 0, upper = 1), "q")
  if (!is.null(tau)) exnet_check(checkmate::check_number(tau, lower = 0, finite = TRUE), "tau")
  # One threshold for both arms, so neither arm is judged by an easier bar.
  t <- if (is.null(tau)) unname(stats::quantile(abs(c(e$w_EE, e$w_RE)), q)) else tau
  e[, tau := t]
  e[, specificity := ifelse(w_EE >= t & w_RE >= t, "both", ifelse(w_EE >= t, "endurance-specific", ifelse(w_RE >= t, "resistance-specific", "neither")))]
  e[]
}

# Printing: a short, readable summary instead of dumping the slots.
setMethod("show", "Embedding", function(object) {
  cat(sprintf("<Embedding> arm %s, %d %s x %d dimensions (%d empty)\n", object@arm, nrow(object@values), object@kind, ncol(object@values),
              sum(colSums(!is.na(object@values)) == 0)))
})
setMethod("show", "PhysicalEdges", function(object) {
  cat(sprintf("<PhysicalEdges> %d undirected edges; sources: %s\n", nrow(object@edges), paste(names(table(object@edges$source)), collapse = ", ")))
})
setMethod("show", "WeightedNetwork", function(object) {
  e <- object@edges
  cat(sprintf("<WeightedNetwork> %d edges (%s)\n", nrow(e), paste(sprintf("%s: %d", names(table(e$edge_type)), as.integer(table(e$edge_type))), collapse = "; ")))
  cat(sprintf("  cor(w_EE, w_RE) = %.3f; mean w_EE = %.4f, mean w_RE = %.4f\n", if (nrow(e) > 2) stats::cor(e$w_EE, e$w_RE) else NA, mean(e$w_EE), mean(e$w_RE)))
})
