#!/usr/bin/env Rscript
# =====================================================================================================
# 12_normalization_comparison.R — STEP 12: HOW DOES THE CHOICE OF NORMALISATION CHANGE THE NETWORKS?
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The edge weights are dot products of normalised log fold changes (steps 1, 1b, 3, 6), so the
#   normalisation decides how much each ome, tissue and arm counts. This script rebuilds every edge weight
#   under FOUR normalisations and shows the endurance-minus-resistance edge differences side by side, so
#   the team can choose one. It reads the unnormalised log fold changes and changes nothing else.
#
#   The four options (the divisor is always computed separately for each ome: RNA, protein, metabolites):
#     1. max, pooled   : divide by the maximum |logFC| over all tissues, timepoints and BOTH arms
#                        (the current step 1 / 1b normalisation)
#     2. max, per arm  : divide each arm by the maximum |logFC| of THAT arm (endurance and resistance get
#                        their own divisor)
#     3. mean, pooled  : divide by the mean |logFC| over all tissues, timepoints and both arms
#     4. mean, per arm : divide each arm by the mean |logFC| of that arm
#   "mean" is the mean ABSOLUTE log fold change (the plain mean of signed changes is close to zero).
#
#   What to keep in mind when comparing:
#     - Pooled options use ONE divisor for both arms, so an arm that genuinely changes more keeps larger
#       values. Per-arm options rescale each arm to its own size, which removes any overall difference in
#       response size between the arms (differences in pattern remain).
#     - For the metabolite network (one ome), options 1 and 3 differ only by a single constant, so their
#       edge-difference pattern is identical; they only differ for genes, where RNA and protein get
#       different relative weights under max vs mean.
#
# OUTPUTS
#   $HACK_OUT/12_normalization_divisors.csv  every divisor: option, ome, arm, value
#   $HACK_OUT/12_normalization_summary.csv   per option and network: correlation between the arms' edge
#       weights, edges changing sign, share of edges higher in endurance, rank correlation of the edge
#       differences with option 1, and overlap of the 20 largest differences with option 1
#   $HACK_OUT/12_threshold_robustness.csv    can a weight cut-off survive the normalisation choice? For
#       each network, weight (w_EE, w_RE, w_diff) and cut-off "keep the top X% of edges by |weight|"
#       (X = 1 ... 50): k, edges kept under all four options, smallest and mean pairwise Jaccard overlap
#       of the kept sets, sign agreement of the edges kept in all four, the cut-off in option 1 units and
#       as a fraction of option 1's largest |weight|. A fixed ABSOLUTE cut-off cannot be compared across
#       options (max- and mean-normalised weights differ ~100-fold), so cut-offs are expressed as ranks.
#   $HACK_OUT/12_threshold_pairwise.csv      the same Jaccard overlap for every pair of options (top 10%, 20%)
#   $HACK_FIG/12a_gene_network_normalization_comparison.png        2 x 2 panels, genes
#   $HACK_FIG/12b_metabolite_network_normalization_comparison.png  2 x 2 panels, metabolites
#   Each panel is the step 11 difference network (red = endurance weight higher, blue = resistance higher,
#   thin grey = same), in the same node positions as figures 10 / 11. Because the options put the weights
#   on different scales, each panel's colours and widths are relative to that panel's own 95th percentile
#   of |w_EE - w_RE|: compare PATTERNS across panels, not raw magnitudes.
#
# TECH STACK
#   R 4.4; data.table, igraph (layout), ggplot2, ggrepel (labels).
#
# INPUTS:  $HACK_OUT: 01_nodes_{EE,RE}_raw_logFC.csv, 01b_metab_nodes_{EE,RE}_raw_logFC.csv,
#          03_weighted_edges.csv, 06_metabolite_edges.csv, 01c_metabolite_ids.csv
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# The same fixed seed as steps 10 / 11, so the layout matches those figures.
SEED <- 20260926
# The four options, in panel order.
OPTIONS <- c("1. max, pooled across arms", "2. max, per arm", "3. mean, pooled across arms", "4. mean, per arm")
# Colours for the difference scale (as in step 11) and text colour.
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"; INK <- "#1A1A1A"

# ---- helpers ----------------------------------------------------------------------------------------
# Read one arm's raw log fold changes as a matrix (rows = nodes, columns = dimensions); id_cols = the
# number of leading identifier columns; columns empty for every node (adipose protein 0.5 / 24 h) dropped.
read_raw <- function(file, id_cols, id_name) {
  x <- fread(file.path(OUT, file))
  m <- as.matrix(x[, -(seq_len(id_cols)), with = FALSE]); rownames(m) <- x[[id_name]]
  m[, colSums(is.na(m)) < nrow(m), drop = FALSE]
}
# The ome of each dimension, from its column name ("muscle_rna_4h" -> "rna").
ome_of <- function(cols) sub("^[a-z]+_([a-z]+)_.*$", "\\1", cols)

# Normalise both arms under one option; returns list(EE = matrix, RE = matrix, div = divisor table).
normalise <- function(E, R, option) {
  # which summary statistic, and whether the two arms share a divisor
  stat <- if (grepl("max", option)) function(v) max(abs(v), na.rm = TRUE) else function(v) mean(abs(v), na.rm = TRUE)
  pooled <- grepl("pooled", option)
  # the ome of every column
  om <- ome_of(colnames(E))
  # one divisor per ome (and per arm, unless pooled)
  div <- rbindlist(lapply(unique(om), function(o) {
    # the values of this ome in each arm
    e <- E[, om == o]; r <- R[, om == o]
    # pooled: one number from both arms; per arm: one number from each arm
    if (pooled) data.table(option = option, ome = o, arm = c("EE", "RE"), divisor = stat(c(e, r)))
    else data.table(option = option, ome = o, arm = c("EE", "RE"), divisor = c(stat(e), stat(r)))
  }))
  # divide every column by its ome's divisor (for that arm)
  dE <- div[arm == "EE"][match(om, ome), divisor]; dR <- div[arm == "RE"][match(om, ome), divisor]
  list(EE = sweep(E, 2, dE, "/"), RE = sweep(R, 2, dR, "/"), div = div)
}

# Edge weights (dot products) for an edge list with node names a, b.
weights <- function(M, a, b) unname(rowSums(M[a, , drop = FALSE] * M[b, , drop = FALSE]))

# Run all four options for one network; returns edge differences per option, the summary and divisors.
compare <- function(E, R, a, b, network) {
  # normalise and weight under each option
  res <- lapply(OPTIONS, function(op) {
    n <- normalise(E, R, op)
    wE <- weights(n$EE, a, b); wR <- weights(n$RE, a, b)
    list(edges = data.table(option = op, a = a, b = b, w_EE = wE, w_RE = wR, w_diff = wE - wR), div = n$div)
  })
  # all edge differences in one long table
  ed <- rbindlist(lapply(res, `[[`, "edges"))
  # reference: option 1 (the current normalisation)
  ref <- ed[option == OPTIONS[1]]
  # the 20 edges with the largest |difference| in the reference
  top_ref <- ref[order(-abs(w_diff))][1:20, paste(a, b)]
  # per-option summary numbers
  summ <- ed[, .(network = network,
                 cor_w_EE_w_RE = round(cor(w_EE, w_RE), 3),
                 edges_sign_change = sum(sign(w_EE) != sign(w_RE)),
                 share_higher_in_EE = round(mean(w_diff > 0), 3),
                 spearman_w_diff_vs_option1 = round(cor(w_diff, ref$w_diff, method = "spearman"), 3),
                 top20_overlap_with_option1 = sum(paste(a, b)[order(-abs(w_diff))][1:20] %in% top_ref)),
             by = option]
  list(edges = ed, summary = summ, div = rbindlist(lapply(res, `[[`, "div"))[, network := network])
}

# Draw the 2 x 2 comparison for one network, in the same layout as steps 10 / 11.
draw_grid <- function(ed, summ, node_order, edge_order, groups, title, file, label_all, width = 15, height = 10) {
  # the network built exactly as in steps 10 / 11 (same vertex order, same edge order, same seed)
  g <- graph_from_data_frame(edge_order, directed = FALSE,
                             vertices = data.table(node = node_order[node_order %in% c(edge_order$a, edge_order$b)]))
  set.seed(SEED); L0 <- layout_with_fr(g)
  # node positions scaled to 0..1
  norm01 <- function(v) (v - min(v)) / diff(range(v))
  pos <- data.table(node = V(g)$name, x = norm01(L0[, 1]), y = norm01(L0[, 2]))
  # edge coordinates, and each panel's difference relative to its own 95th percentile of |w_diff|
  E <- copy(ed)
  E[, `:=`(x = pos$x[match(a, pos$node)], y = pos$y[match(a, pos$node)],
           xend = pos$x[match(b, pos$node)], yend = pos$y[match(b, pos$node)])]
  E[, rel := w_diff / quantile(abs(w_diff), 0.95), by = option]
  # draw small differences first so large ones sit on top
  E <- E[order(option, abs(rel))]
  # per-panel node strength difference (sum of |w_EE| minus sum of |w_RE|), relative like the edges
  both <- rbind(E[, .(option, node = a, w_EE, w_RE)], E[, .(option, node = b, w_EE, w_RE)])
  N <- both[, .(delta = sum(abs(w_EE)) - sum(abs(w_RE))), by = .(option, node)]
  N <- merge(N, pos, by = "node")
  N[, rel_delta := abs(delta) / max(abs(delta)), by = option]
  # labels: every node (metabolites) or the 6 largest strength differences per panel (genes)
  N[, rk := frank(-abs(delta), ties.method = "first"), by = option]
  N[, lab := fifelse(label_all | rk <= 6, node, NA_character_)]
  # class labels above each connected group (metabolites only)
  GL <- NULL
  if (!is.null(groups)) {
    comp <- components(g)$membership
    GL <- merge(N[, .(option, node, x, y)], data.table(node = names(comp), comp = comp), by = "node")[
      , .(x = mean(range(x)), y = max(y) + 0.07, lab = groups[node[1]]), by = .(option, comp)]
  }
  # panel titles with the key numbers
  summ <- copy(summ)[, strip := sprintf("%s\ncor(EE, RE) = %.2f  |  sign changes = %d  |  top-20 overlap with option 1 = %d/20",
                                       option, cor_w_EE_w_RE, edges_sign_change, top20_overlap_with_option1)]
  lab_of <- setNames(summ$strip, summ$option)
  E[, panel := factor(lab_of[option], levels = lab_of[OPTIONS])]
  N[, panel := factor(lab_of[option], levels = lab_of[OPTIONS])]
  if (!is.null(GL)) GL[, panel := factor(lab_of[option], levels = lab_of[OPTIONS])]
  # the plot
  p <- ggplot() +
    # edges: colour and width = difference relative to the panel's 95th percentile
    geom_segment(data = E, aes(x, y, xend = xend, yend = yend, colour = rel, linewidth = abs(rel)), lineend = "round") +
    # nodes: grey, size = relative |strength difference|
    geom_point(data = N, aes(x, y, size = rel_delta), shape = 21, fill = "grey80", colour = "grey30", stroke = 0.2) +
    # node labels that avoid each other
    geom_text_repel(data = N[!is.na(lab)], aes(x, y, label = lab), size = if (label_all) 1.7 else 2, colour = "grey15",
                    min.segment.length = 0.2, segment.size = 0.12, max.overlaps = Inf, seed = SEED) +
    # class labels (metabolites)
    (if (!is.null(GL)) geom_text(data = GL, aes(x, y, label = lab), size = 2.3, colour = "grey35", fontface = "bold.italic") else NULL) +
    # one panel per option
    facet_wrap(~ panel, ncol = 2) +
    # blue - grey - red, relative to each panel's 95th percentile
    scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                           breaks = c(-1, 0, 1), labels = c("higher in resistance", "same", "higher in endurance"),
                           name = "edge difference w_EE − w_RE\n(relative to each panel's 95th percentile)") +
    scale_linewidth(range = c(0.08, 1.4), guide = "none") +
    scale_size(range = c(0.5, 3.5), guide = "none") +
    guides(colour = guide_colourbar(barwidth = unit(7, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5)) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.1), clip = "off") +
    labs(title = title) +
    theme_classic(base_size = 8) %+replace% theme(
      axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
      strip.background = element_rect(fill = "#EDF0F2", colour = NA), strip.text = element_text(size = 7, face = "bold", colour = "#20262B", margin = margin(3, 3, 3, 3)),
      legend.position = "bottom", legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 9.5, hjust = 0, colour = INK, margin = margin(b = 5)),
      plot.background = element_rect(fill = "white", colour = NA))
  # save
  ggsave(file, p, width = width, height = height, dpi = 300, bg = "white")
  message("-> ", file)
}

# ---- genes --------------------------------------------------------------------------------------
# Raw log fold changes, both arms, rows named by gene symbol (16 observed dimensions).
GE <- read_raw("01_nodes_EE_raw_logFC.csv", 2, "gene_symbol"); GR <- read_raw("01_nodes_RE_raw_logFC.csv", 2, "gene_symbol")
# The gene edges, in the step 3 order (the order steps 10 / 11 use).
ge <- fread(file.path(OUT, "03_weighted_edges.csv"))[, .(a = symbol_a, b = symbol_b)]
# Compare the four options.
G <- compare(GE, GR, ge$a, ge$b, "genes")
# Safety check: option 1 reproduces the current step 3 weights.
stopifnot(isTRUE(all.equal(G$edges[option == OPTIONS[1], w_EE], fread(file.path(OUT, "03_weighted_edges.csv"))$w_EE)))
# Draw the gene grid.
draw_grid(G$edges, G$summary, fread(file.path(OUT, "01_nodes_EE.csv"))$gene_symbol, ge, NULL,
          "Gene network: endurance minus resistance edge weights under four normalisations (471 genes, STRING >= 700)",
          file.path(FIG, "12a_gene_network_normalization_comparison.png"), label_all = FALSE)

# ---- metabolites ---------------------------------------------------------------------------------
# Raw log fold changes, both arms (9 dimensions; all columns are "metab").
ME <- read_raw("01b_metab_nodes_EE_raw_logFC.csv", 1, "metabolite"); MR <- read_raw("01b_metab_nodes_RE_raw_logFC.csv", 1, "metabolite")
# The metabolite edges, in the step 6 order.
me <- fread(file.path(OUT, "06_metabolite_edges.csv"))[, .(a = metabolite_a, b = metabolite_b)]
# Compare the four options.
M <- compare(ME, MR, me$a, me$b, "metabolites")
# Safety check: option 1 reproduces the current step 6 weights.
stopifnot(isTRUE(all.equal(M$edges[option == OPTIONS[1], w_EE], fread(file.path(OUT, "06_metabolite_edges.csv"))$w_EE)))
# Each metabolite's super class, for the group labels.
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv")); cls <- setNames(ids$super_class, ids$metabolite)
# Draw the metabolite grid.
draw_grid(M$edges, M$summary, fread(file.path(OUT, "01b_metab_nodes_EE.csv"))$metabolite, me, cls,
          "Metabolite network: endurance minus resistance edge weights under four normalisations (Rhea + STRING-linked enzymes, same super class)",
          file.path(FIG, "12b_metabolite_network_normalization_comparison.png"), label_all = TRUE)

# ---- threshold robustness: which cut-off keeps the same edges under every normalisation? ------------
# Jaccard overlap of two edge sets (shared / all).
jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
# Cut-offs as "top X% of edges by |weight|" (ranks are comparable across options; raw values are not).
TOP_PCT <- c(1, 2, 5, 10, 15, 20, 25, 30, 40, 50)
# Helper: the edges kept by each option at the top-k cut-off for one weight column.
kept <- function(ed, col, k) lapply(OPTIONS, function(o) { x <- ed[option == o]; paste(x$a, x$b)[order(-abs(x[[col]]))][seq_len(k)] })
# Helper: one robustness row per cut-off for one network and weight column.
robust <- function(ed, net, col) rbindlist(lapply(TOP_PCT, function(pct) {
  # number of edges kept (at least 1)
  k <- max(1, round(pct / 100 * uniqueN(paste(ed$a, ed$b))))
  # kept sets under the four options, and the edges kept by all four
  sets <- kept(ed, col, k); common <- Reduce(intersect, sets)
  # smallest and mean overlap over the 6 pairs of options
  pj <- combn(length(OPTIONS), 2, function(ix) jac(sets[[ix[1]]], sets[[ix[2]]]))
  # of the edges kept by all four, the share whose sign is the same under all four
  sg <- ed[paste(a, b) %in% common, .(same = uniqueN(sign(get(col))) == 1), by = .(a, b)][, mean(same)]
  # option 1's |weight| at the cut-off (its own units) and as a share of option 1's largest |weight|
  r1 <- sort(abs(ed[option == OPTIONS[1]][[col]]), decreasing = TRUE)
  data.table(network = net, weight = col, top_pct = pct, k = k, kept_by_all_4 = length(common),
             min_jaccard = round(min(pj), 2), mean_jaccard = round(mean(pj), 2),
             sign_agreement = if (length(common)) round(sg, 2) else NA_real_,
             option1_cutoff = signif(r1[k], 2), option1_cutoff_share_of_max = round(r1[k] / r1[1], 2))
}))
# Every network x weight combination.
thr <- rbindlist(lapply(c("w_diff", "w_EE", "w_RE"), function(col) rbind(robust(G$edges, "genes", col), robust(M$edges, "metabolites", col))))
fwrite(thr, file.path(OUT, "12_threshold_robustness.csv"))
# Pairwise overlaps (which options disagree) at the top 10% and 20% of |w_diff|.
pw <- rbindlist(lapply(list(genes = G$edges, metabolites = M$edges), function(ed) rbindlist(lapply(c(10, 20), function(pct) {
  k <- round(pct / 100 * uniqueN(paste(ed$a, ed$b))); sets <- kept(ed, "w_diff", k)
  rbindlist(combn(length(OPTIONS), 2, function(ix) data.table(top_pct = pct, option_a = OPTIONS[ix[1]], option_b = OPTIONS[ix[2]],
                                                               jaccard = round(jac(sets[[ix[1]]], sets[[ix[2]]]), 2)), simplify = FALSE))
}))), idcol = "network")
fwrite(pw, file.path(OUT, "12_threshold_pairwise.csv"))
# Show the w_diff rows.
print(thr[weight == "w_diff"])

# ---- tables ---------------------------------------------------------------------------------------
# All divisors.
fwrite(rbind(G$div, M$div), file.path(OUT, "12_normalization_divisors.csv"))
# All summary numbers.
summ <- rbind(G$summary, M$summary)
fwrite(summ, file.path(OUT, "12_normalization_summary.csv"))
# Show them.
print(summ)
