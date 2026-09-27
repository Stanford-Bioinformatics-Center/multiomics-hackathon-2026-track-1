# =====================================================================================================
# build.R — constructors: typed inputs in, a validated weighted network out
# =====================================================================================================
# THE RULE THIS ENGINE IMPLEMENTS (hard gate x soft weight)
#   For molecules u and v and exercise arm a (EE or RE):
#       edge(u, v) exists   <=>  a physical database links u and v            (HARD: STRING >= 700, Rhea)
#       w_a(u, v)           =    sum over dimensions d of z_u,d(a) * z_v,d(a)  (SOFT: dot product)
#   where z_u(a) is u's normalised exercise-response vector in arm a (one value per tissue x ome x time).
#   Similar responses without a physical link never create an edge; a physical link between molecules
#   that respond differently gets a weight near zero or negative. w_diff = w_EE - w_RE compares the arms.
#   Metabolite - protein edges: the metabolite's vector (one value per tissue x time) is DOUBLED so each
#   value sits under both the RNA and the protein column of the same tissue and time, then one dot product
#   is taken with the gene's vector (the metabolite counts equally against RNA and protein).
# =====================================================================================================

#' Create a validated Embedding (one arm's response vectors).
#' @param values numeric matrix, rows named by molecule, columns named by response dimension.
#' @param arm `character(1)` "EE" or "RE".
#' @param kind `character(1)` label for the rows ("gene", "metabolite", ...).
#' @return an `Embedding`.
embedding <- function(values, arm, kind = "gene") {
  # Accept a data.frame by converting it, but insist on a numeric matrix underneath.
  if (is.data.frame(values)) values <- as.matrix(values)
  exnet_check(checkmate::check_matrix(values, mode = "numeric", min.rows = 1, min.cols = 1), "values")
  exnet_check(checkmate::check_choice(arm, c("EE", "RE")), "arm")
  exnet_check(checkmate::check_string(kind, min.chars = 1), "kind")
  # new() runs the class validity rules (names, finiteness, where NAs may be).
  new("Embedding", values = values, arm = arm, kind = kind)
}

#' Create validated PhysicalEdges (the hard layer).
#' @param a,b `character` molecule IDs at each end (same length).
#' @param source `character` database supporting each edge (recycled if length 1).
#' @param evidence `character` free-text evidence (e.g. STRING score, Rhea reaction IDs); optional.
#' @return a `PhysicalEdges` object.
physical_edges <- function(a, b, source, evidence = "") {
  exnet_check(checkmate::check_character(a, any.missing = FALSE, min.len = 1), "a")
  exnet_check(checkmate::check_character(b, any.missing = FALSE, len = length(a)), "b")
  exnet_check(checkmate::check_character(source, any.missing = FALSE), "source")
  # Build the edge table (recycling a single source / evidence value over all edges).
  e <- data.frame(a = a, b = b, source = rep_len(source, length(a)), evidence = rep_len(as.character(evidence), length(a)), stringsAsFactors = FALSE)
  new("PhysicalEdges", edges = e)
}

#' Double a metabolite embedding into the gene dimension space (for metabolite - protein edges).
#' @param metab `Embedding` of metabolites, columns named "<tissue>_metab_<time>".
#' @param genes `Embedding` of genes (same arm), columns named "<tissue>_<rna|prot>_<time>".
#' @return an `Embedding` of the metabolites with the genes' columns: each metabolite value sits under both
#'   the RNA and the protein column of its tissue and time. A gene column that is empty for every gene
#'   (a dimension nobody measured, e.g. adipose protein at 0.5 h) is left empty here too, so it can never
#'   contribute to a metabolite - protein dot product.
double_embedding <- function(metab, genes) {
  if (!is(metab, "Embedding") || !is(genes, "Embedding")) exnet_abort("'metab' and 'genes' must be Embedding objects")
  if (metab@arm != genes@arm) exnet_abort("'metab' and 'genes' must be the same arm", "exnet_arm_error")
  gene_dims <- colnames(genes@values)
  # For each gene column, the metabolite column with the same tissue and time.
  src <- sub("_(rna|prot)_", "_metab_", gene_dims)
  miss <- setdiff(src, colnames(metab@values))
  if (length(miss)) exnet_abort(sprintf("no metabolite column for gene dimension(s): %s", paste(unique(miss), collapse = ", ")), "exnet_dimension_error")
  # Copy the columns (each metabolite column used once per ome) and rename them to the gene columns.
  v <- metab@values[, src, drop = FALSE]; colnames(v) <- gene_dims
  # Blank the dimensions that no gene has.
  empty <- colSums(!is.na(genes@values)) == 0; v[, empty] <- NA_real_
  embedding(v, metab@arm, "metabolite (doubled)")
}

#' Row-bind embeddings of the same arm and dimensions (e.g. doubled metabolites + genes).
#' @param ... `Embedding` objects.
#' @return one `Embedding`.
bind_embeddings <- function(...) {
  parts <- list(...)
  if (!all(vapply(parts, is, logical(1), "Embedding"))) exnet_abort("bind_embeddings() takes Embedding objects")
  if (length(unique(vapply(parts, function(p) p@arm, ""))) != 1) exnet_abort("all embeddings must be the same arm", "exnet_arm_error")
  if (length(unique(lapply(parts, function(p) colnames(p@values)))) != 1) exnet_abort("all embeddings must have the same dimensions", "exnet_dimension_error")
  embedding(do.call(rbind, lapply(parts, function(p) p@values)), parts[[1]]@arm, paste(unique(vapply(parts, function(p) p@kind, "")), collapse = " + "))
}

# Internal: the dot product of the two ends of every edge in one arm; empty (all-missing) dimensions add nothing.
.edge_dots <- function(E, a, b) {
  # Every edge end must have a vector; a molecule without one would silently get weight 0.
  miss <- setdiff(unique(c(a, b)), rownames(E@values))
  if (length(miss)) exnet_abort(sprintf("%d edge molecule(s) have no response vector in arm %s, e.g. %s", length(miss), E@arm, paste(head(miss, 5), collapse = ", ")), "exnet_missing_node_error")
  # Element-wise product of the two vectors, summed over dimensions (na.rm skips only whole empty columns,
  # which the Embedding validity rule guarantees are the only missing values).
  unname(rowSums(E@values[a, , drop = FALSE] * E@values[b, , drop = FALSE], na.rm = TRUE))
}

#' Build a weighted network: hard (physical) edges, soft (dot-product) weights, both arms.
#' @param hard `PhysicalEdges` the edges that are allowed to exist.
#' @param EE,RE `Embedding` response vectors of every molecule on the edges, endurance and resistance arm,
#'   with identical rows and columns (same molecules, same dimensions, same order).
#' @param edge_type `character(1)` label for these edges, e.g. "protein - protein".
#' @return a validated `WeightedNetwork`.
build_network <- function(hard, EE, RE, edge_type = "edge") {
  if (!is(hard, "PhysicalEdges")) exnet_abort("'hard' must be a PhysicalEdges object (use physical_edges())")
  if (!is(EE, "Embedding") || !is(RE, "Embedding")) exnet_abort("'EE' and 'RE' must be Embedding objects (use embedding())")
  if (EE@arm != "EE" || RE@arm != "RE") exnet_abort("'EE' must be the endurance arm and 'RE' the resistance arm", "exnet_arm_error")
  # The two arms must describe the same molecules in the same dimension space, or w_diff would compare unlike things.
  if (!identical(colnames(EE@values), colnames(RE@values))) exnet_abort("EE and RE must have identical dimensions in the same order", "exnet_dimension_error")
  if (!identical(sort(rownames(EE@values)), sort(rownames(RE@values)))) exnet_abort("EE and RE must cover the same molecules", "exnet_dimension_error")
  exnet_check(checkmate::check_string(edge_type, min.chars = 1), "edge_type")
  e <- hard@edges
  # The dimensions actually used: those not empty for every molecule.
  used <- colnames(EE@values)[colSums(!is.na(EE@values)) > 0]
  # SOFT weights, per arm, for every HARD edge (and only for those).
  w_EE <- .edge_dots(EE, e$a, e$b); w_RE <- .edge_dots(RE, e$a, e$b)
  out <- data.frame(a = e$a, b = e$b, edge_type = edge_type, source = e$source, evidence = e$evidence,
                    w_EE = w_EE, w_RE = w_RE, w_diff = w_EE - w_RE, n_dims = length(used), stringsAsFactors = FALSE)
  new("WeightedNetwork", edges = out, dims = setNames(list(used), edge_type))
}

#' Combine several weighted networks (e.g. protein-protein, metabolite-metabolite, metabolite-protein) into one.
#' @param ... `WeightedNetwork` objects.
#' @return one validated `WeightedNetwork` (its validity rules are re-checked on the union).
combine_networks <- function(...) {
  parts <- list(...)
  if (!length(parts) || !all(vapply(parts, is, logical(1), "WeightedNetwork"))) exnet_abort("combine_networks() takes WeightedNetwork objects")
  new("WeightedNetwork", edges = do.call(rbind, lapply(parts, function(p) p@edges)), dims = do.call(c, lapply(parts, function(p) p@dims)))
}
