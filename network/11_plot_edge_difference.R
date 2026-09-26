#!/usr/bin/env Rscript
# =====================================================================================================
# 11_plot_edge_difference.R — STEP 11: ONE NETWORK SHOWING WHERE ENDURANCE AND RESISTANCE DIFFER
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Step 10 draws the endurance and resistance networks one above the other. This script compresses the
#   two into ONE network per data type, in which each edge shows the DIFFERENCE between the arms:
#       w_diff = w_EE - w_RE   (endurance edge weight minus resistance edge weight, steps 3 / 6)
#     - positive (endurance weight higher)  -> red, and thicker the bigger the difference;
#     - negative (resistance weight higher) -> blue, and thicker the bigger the difference;
#     - near zero (the same in both arms)   -> thin and light grey.
#   READ WITH CARE: the difference is of SIGNED weights. For an edge that is positive in both arms, red means
#   the endurance co-response is stronger. For an edge that is NEGATIVE in both arms (the two nodes respond
#   in opposite directions), red means the endurance edge is closer to zero, i.e. the RESISTANCE edge is the
#   more strongly negative one (e.g. IL18-CCL5: -9.2 in endurance, -54.0 in resistance -> red). The legend
#   therefore says "higher", not "stronger"; the per-edge signs are in 04_edge_diff.csv / 04b_metab_edge_diff.csv
#   (columns w_EE, w_RE, higher_in, larger_abs_in).
#   Two figures:
#     11a_gene_network_edge_difference.png        genes (STRING edges, 471-gene universe)
#     11b_metabolite_network_edge_difference.png  metabolites (Rhea + super-class edges, step 6)
#   Node positions are the SAME as in figures 10a / 10b (same layout code, same fixed seed, same edges),
#   so the figures can be read side by side.
#
# HOW EACH ELEMENT IS DRAWN
#   - Edge colour: w_diff on a blue - grey - red scale centred at 0. The colour limits are symmetric at
#     the 95th percentile of |w_diff|, so a few extreme edges do not wash out the rest; edges beyond it are
#     drawn in the strongest colour (their exact values are in 04_edge_diff.csv / 04b_metab_edge_diff.csv).
#   - Edge width: |w_diff| (thin = little difference).
#   - Edge line type: solid if the difference passes the bootstrap test at uncorrected p < 0.05 (steps
#     4 / 4b), dotted and faded otherwise (so large differences that are within noise do not dominate). This separates differences larger than measurement noise from those that
#     are not; p < 0.05 without correction is hypothesis-level only (see README for the FDR results).
#   - Nodes: grey (colour is reserved for the edges), size = |difference in node strength| between arms
#     (step 4 / 4b), labels for the 10 genes with the largest strength difference / every metabolite.
#     Genes: triangle = strength differs at uncorrected p < 0.05. Metabolites: shape = RefMet super class.
#   - Only nodes with at least one edge are drawn (286 genes, 44 metabolites).
#   Titles are descriptive only; interpretation belongs in the report text.
#
# TECH STACK
#   R 4.4; data.table, igraph (layout), ggplot2, ggrepel (labels).
#
# INPUTS:  $HACK_OUT: 01_nodes_EE.csv, 01b_metab_nodes_EE.csv, 01c_metabolite_ids.csv,
#          04_edge_diff.csv, 04_node_diff.csv, 04b_metab_edge_diff.csv, 04b_metab_node_diff.csv
# OUTPUTS: $HACK_FIG (default ~/Desktop/output/hackathon; never inside the repo)
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# The same fixed seed as step 10, so the layout is identical to figures 10a / 10b.
SEED <- 20260926
# How many genes to label (those with the largest strength difference).
NLAB <- 10
# Text colour.
INK <- "#1A1A1A"
# Colours for the difference scale: blue (resistance stronger) - light grey (same) - red (endurance stronger).
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"

# The theme: small text, no axes (a network has no meaningful axes), legends at the bottom.
theme_net <- function(base = 8) {
  theme_classic(base_size = base) %+replace% theme(
    text = element_text(colour = INK, size = base),
    axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    legend.title = element_text(colour = INK, size = base * 0.8),
    legend.text = element_text(colour = INK, size = base * 0.74),
    legend.background = element_blank(), legend.key = element_blank(),
    plot.title = element_text(colour = INK, face = "bold", size = base * 1.18, hjust = 0,
                              margin = margin(b = base * 0.55)),
    plot.background = element_rect(fill = "white", colour = NA))
}
# Scale a coordinate to 0..1 (a single value goes to the middle).
norm01 <- function(v) if (diff(range(v)) < 1e-9) rep(0.5, length(v)) else (v - min(v)) / diff(range(v))

# Build and draw one difference network.
#   node_order: all node names in the same order step 10 used (so the layout matches 10a / 10b)
#   edges: data.table(a, b, w_diff, p_boot); nodes: data.table(node, delta, shape_key)
draw_diff <- function(node_order, edges, nodes, title, shape_values, shape_name, label_rule, file,
                      width = 10.5, height = 5.4) {
  # the network over the connected nodes, built exactly as in step 10 (same vertex order, same edge order)
  g <- graph_from_data_frame(edges[, .(a, b)], directed = FALSE,
                             vertices = data.table(node = node_order[node_order %in% c(edges$a, edges$b)]))
  # the same layout as step 10 (same seed)
  set.seed(SEED); L0 <- layout_with_fr(g)
  # node positions scaled to 0..1
  pos <- data.table(node = V(g)$name, x = norm01(L0[, 1]), y = norm01(L0[, 2]))
  # edge segments with their difference and test result
  E <- copy(edges)
  E[, `:=`(x = pos$x[match(a, pos$node)], y = pos$y[match(a, pos$node)],
           xend = pos$x[match(b, pos$node)], yend = pos$y[match(b, pos$node)])]
  # solid for differences beyond measurement noise at uncorrected p < 0.05, dotted otherwise
  E[, tested := factor(fifelse(p_boot < 0.05, "p < 0.05 (uncorrected)", "not significant"),
                       levels = c("p < 0.05 (uncorrected)", "not significant"))]
  # drawing order: non-significant edges first, then smallest to largest difference, so the significant
  # and large differences sit on top (base-R order: data.table's setorder cannot sort on abs())
  E <- E[order(-as.integer(tested), abs(w_diff))]
  # symmetric colour limit at the 95th percentile of |w_diff|
  lim <- as.numeric(quantile(abs(E$w_diff), 0.95))
  # node table with positions, the strength difference and the shape key
  N <- merge(pos, nodes, by = "node")
  # labels for the chosen nodes
  N[, rk := frank(-abs(delta), ties.method = "first")]
  N[, lab := fifelse(label_rule(rk), node, NA_character_)]
  # the plot, back to front
  p <- ggplot() +
    # edges: colour and width = the difference; line type and opacity = bootstrap result
    geom_segment(data = E, aes(x, y, xend = xend, yend = yend, colour = w_diff, linewidth = abs(w_diff), linetype = tested,
                               alpha = tested), lineend = "round") +
    # opacity: significant edges fully drawn, the rest faded
    scale_alpha_manual(values = c("p < 0.05 (uncorrected)" = 1, "not significant" = 0.4), guide = "none") +
    # nodes: grey, size = |strength difference|, shape = shape key
    geom_point(data = N, aes(x, y, size = abs(delta), shape = shape_key), fill = "grey80", colour = "grey30", stroke = 0.22) +
    # labels that avoid each other and the nodes
    geom_text_repel(data = N[!is.na(lab)], aes(x, y, label = lab), size = 2.3, colour = "grey15",
                    min.segment.length = 0.2, segment.size = 0.15, max.overlaps = Inf, seed = SEED) +
    # blue - grey - red scale centred at zero
    scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-lim, lim),
                           oob = scales::squish, name = "edge difference  w_EE − w_RE",
                           breaks = c(-lim, 0, lim),
                           labels = c("higher in resistance\n(w_RE > w_EE)", "same", "higher in endurance\n(w_EE > w_RE)")) +
    # widths: thin for no difference, thick for large differences
    scale_linewidth(range = c(0.1, 1.8), guide = "none") +
    # line types for the bootstrap result
    scale_linetype_manual(values = c("p < 0.05 (uncorrected)" = "solid", "not significant" = "12"),
                          name = "difference vs noise (bootstrap)", drop = FALSE) +
    # node sizes and shapes
    scale_size(range = c(0.8, 4.5), name = "node strength difference (absolute)") +
    scale_shape_manual(values = shape_values, name = shape_name) +
    # legends
    guides(colour = guide_colourbar(order = 1, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"),
                                    title.position = "top", title.hjust = 0.5),
           linetype = guide_legend(order = 2, override.aes = list(colour = "grey30", linewidth = 0.6)),
           shape = guide_legend(order = 3, override.aes = list(size = 2.6)), size = guide_legend(order = 4)) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02), clip = "off") +
    labs(title = title) + theme_net()
  # save with a white background at 300 dpi
  ggsave(file, p, width = width, height = height, dpi = 300, bg = "white")
  message(sprintf("-> %s (%d nodes, %d edges; colour limit +-%.1f)", file, vcount(g), ecount(g), lim))
}

# ---- figure 11a: genes ---------------------------------------------------------------------------
# Node order exactly as step 10 used it (gene symbols in the order of the step 1 node table).
gorder <- fread(file.path(OUT, "01_nodes_EE.csv"))$gene_symbol
# Edge differences and bootstrap p-values (step 4), in the step 3 edge order used by step 10.
w3 <- fread(file.path(OUT, "03_weighted_edges.csv"))[, .(a = symbol_a, b = symbol_b)]
# (the differences and p-values from step 4)
d4 <- fread(file.path(OUT, "04_edge_diff.csv"))[, .(a = symbol_a, b = symbol_b, w_diff, p_boot)]
# (merge back onto the step 3 order so the layout is computed on the same edge sequence as step 10)
ge <- d4[w3, on = .(a, b)]
# Safety check: every edge found its difference.
stopifnot(!anyNA(ge$w_diff))
# Per-gene strength difference and whether it passes p < 0.05 (step 4).
nd <- fread(file.path(OUT, "04_node_diff.csv"))
# (one row per gene: name, strength difference, shape key)
gn <- data.table(node = nd$gene_symbol, delta = nd$delta_strength,
                 shape_key = factor(fifelse(nd$p_boot < 0.05, "strength differs EE vs RE (p < 0.05, uncorrected)", "no nominal difference"),
                                    levels = c("no nominal difference", "strength differs EE vs RE (p < 0.05, uncorrected)")))
# Draw and save the gene figure.
draw_diff(gorder, ge, gn, "Gene network: endurance minus resistance edge weights (471 genes, STRING combined score >= 700)",
          shape_values = c("no nominal difference" = 22, "strength differs EE vs RE (p < 0.05, uncorrected)" = 24),
          shape_name = "node", label_rule = function(rk) rk <= NLAB,
          file = file.path(FIG, "11a_gene_network_edge_difference.png"))

# ---- figure 11b: metabolites -----------------------------------------------------------------------
# Node order exactly as step 10 used it (metabolites in the order of the step 1b node table).
morder <- fread(file.path(OUT, "01b_metab_nodes_EE.csv"))$metabolite
# Edge differences and bootstrap p-values (step 4b), in the step 6 edge order used by step 10.
w6 <- fread(file.path(OUT, "06_metabolite_edges.csv"))[, .(a = metabolite_a, b = metabolite_b)]
# (the differences and p-values from step 4b)
d4b <- fread(file.path(OUT, "04b_metab_edge_diff.csv"))[, .(a = metabolite_a, b = metabolite_b, w_diff, p_boot)]
# (merge back onto the step 6 order)
me <- d4b[w6, on = .(a, b)]
# Safety check: every edge found its difference.
stopifnot(!anyNA(me$w_diff))
# Per-metabolite strength difference (step 4b) and super class (step 1c).
nb <- fread(file.path(OUT, "04b_metab_node_diff.csv"))
# (the metabolite classes from step 1c)
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
# (one row per metabolite: name, strength difference, super class as shape)
mn <- data.table(node = nb$metabolite, delta = nb$delta_strength, shape_key = factor(ids$super_class[match(nb$metabolite, ids$metabolite)]))
# The super classes present get distinct shapes (same order as step 10).
cls <- sort(unique(as.character(mn$shape_key)))
# Draw and save the metabolite figure.
draw_diff(morder, me, mn, "Metabolite network: endurance minus resistance edge weights (shared Rhea enzyme among the 471 genes + same RefMet super class)",
          shape_values = setNames(c(21, 22, 24, 23, 25)[seq_along(cls)], cls), shape_name = "RefMet super class",
          label_rule = function(rk) rep(TRUE, length(rk)),
          file = file.path(FIG, "11b_metabolite_network_edge_difference.png"))
