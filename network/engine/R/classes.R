# =====================================================================================================
# classes.R — the three typed S4 classes of the engine, each with a validity method
# =====================================================================================================
# Embedding       one exercise arm's normalised response vectors: one row per molecule, one column per
#                 response dimension (tissue x ome x time). The SOFT information.
# PhysicalEdges   pairs of molecules linked by a physical-interaction database (STRING, Rhea). The HARD
#                 information: it decides WHETHER an edge exists.
# WeightedNetwork the result: every physical edge with its weight in each arm and their difference.
# S4 slots are typed, and validity() runs whenever an object is created or modified with validObject(), so an
# object that exists is guaranteed to satisfy these rules.
# =====================================================================================================

#' @slot values numeric matrix, rows = molecules (named), columns = response dimensions (named).
#' @slot arm `character(1)`, "EE" (endurance) or "RE" (resistance).
#' @slot kind `character(1)`, what the rows are, e.g. "gene" or "metabolite".
setClass("Embedding", slots = c(values = "matrix", arm = "character", kind = "character"))

# Validity rules for an Embedding.
setValidity("Embedding", function(object) {
  v <- object@values; msgs <- character()
  # The matrix must hold numbers.
  if (!is.numeric(v)) msgs <- c(msgs, "values must be numeric")
  # Rows and columns must be named, uniquely, so molecules and dimensions can be matched by name.
  if (is.null(rownames(v)) || anyDuplicated(rownames(v)) || any(rownames(v) == "")) msgs <- c(msgs, "rows must have unique, non-empty names (molecule IDs)")
  if (is.null(colnames(v)) || anyDuplicated(colnames(v))) msgs <- c(msgs, "columns must have unique names (response dimensions)")
  # The arm label must be one of the two exercise arms.
  if (length(object@arm) != 1L || !object@arm %in% c("EE", "RE")) msgs <- c(msgs, "arm must be 'EE' or 'RE'")
  # Infinite values would make every dot product meaningless.
  if (is.numeric(v) && any(is.infinite(v))) msgs <- c(msgs, "values must not be infinite")
  # Missing values are allowed ONLY as whole empty columns (a dimension no molecule has, e.g. adipose
  # protein at 0.5 h); a gap inside a column would silently change some dot products and not others.
  if (is.numeric(v)) { part <- colSums(is.na(v)) > 0 & colSums(is.na(v)) < nrow(v)
    if (any(part)) msgs <- c(msgs, sprintf("missing values allowed only as whole empty columns; partly missing: %s", paste(colnames(v)[part], collapse = ", "))) }
  if (length(msgs)) msgs else TRUE
})

#' @slot edges data.frame with columns a, b (molecule IDs), source (database) and evidence (free text).
setClass("PhysicalEdges", slots = c(edges = "data.frame"))

# Validity rules for PhysicalEdges.
setValidity("PhysicalEdges", function(object) {
  e <- object@edges; msgs <- character()
  # The required columns must be present.
  need <- c("a", "b", "source", "evidence"); miss <- setdiff(need, names(e))
  if (length(miss)) return(sprintf("edges is missing column(s): %s", paste(miss, collapse = ", ")))
  # A molecule cannot interact with itself in this network.
  if (any(e$a == e$b)) msgs <- c(msgs, "self-loops are not allowed")
  # Edges are undirected: A-B and B-A are the same edge and may appear only once.
  key <- ifelse(e$a < e$b, paste(e$a, e$b, sep = "\r"), paste(e$b, e$a, sep = "\r"))
  if (anyDuplicated(key)) msgs <- c(msgs, "duplicate undirected edges (A-B listed twice or as B-A)")
  # Every edge must name the database that supports it.
  if (any(is.na(e$source) | e$source == "")) msgs <- c(msgs, "every edge needs a source database")
  if (length(msgs)) msgs else TRUE
})

#' @slot edges data.frame: a, b, edge_type, source, w_EE, w_RE, w_diff, n_dims.
#' @slot dims list of the dimension names used by each edge type.
setClass("WeightedNetwork", slots = c(edges = "data.frame", dims = "list"))

# Validity rules for a WeightedNetwork.
setValidity("WeightedNetwork", function(object) {
  e <- object@edges; msgs <- character()
  need <- c("a", "b", "edge_type", "source", "w_EE", "w_RE", "w_diff", "n_dims"); miss <- setdiff(need, names(e))
  if (length(miss)) return(sprintf("edges is missing column(s): %s", paste(miss, collapse = ", ")))
  # Every weight must be a finite number.
  if (any(!is.finite(c(e$w_EE, e$w_RE)))) msgs <- c(msgs, "weights must be finite")
  # The difference must be exactly EE minus RE (guards against columns being swapped).
  if (any(abs(e$w_diff - (e$w_EE - e$w_RE)) > 1e-12)) msgs <- c(msgs, "w_diff must equal w_EE - w_RE")
  if (length(msgs)) msgs else TRUE
})
