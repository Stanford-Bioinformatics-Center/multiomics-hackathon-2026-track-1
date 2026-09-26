#!/usr/bin/env Rscript
# =====================================================================================================
# 10_plot_arm_networks.R — STEP 10: PICTURES OF THE ENDURANCE AND RESISTANCE NETWORKS, ONE ABOVE THE OTHER
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Draws two figures, in the layered style of week 5's figure 3.1 (El-Kebir et al. 2015, Figure 4), with
#   the two networks stacked:
#     10a_gene_networks_EE_vs_RE.png        the gene networks (STRING edges, step 3 weights)
#     10b_metabolite_networks_EE_vs_RE.png  the metabolite networks (Rhea + class rule, step 6 weights)
#   Top layer = endurance (EE), bottom layer = resistance (RE). Both layers have the SAME edges (the edges
#   come from databases, not from the exercise data); what differs is how strong each edge and node is.
#
# HOW EACH ELEMENT IS DRAWN
#   - Node position: ONE layout, computed from the shared edges (Fruchterman-Reingold, fixed seed), is
#     used for BOTH layers, so every node sits in the same place in the endurance and resistance layer.
#     Unlike figure 3.1 there are no lines joining the layers: figure 3.1 needs them because human and rat
#     networks contain different genes (orthology), whereas here both networks contain exactly the same
#     nodes and edges by construction, so such lines would carry no information. With identical positions,
#     the eye compares what actually differs between arms: edge widths and signs, node sizes.
#   - Node colour: the node's mean response across all its dimensions (scaled logFC from steps 1 / 1b;
#     violet = down, white = little change, orange = up; limits -2 to 2, values beyond are shown at the
#     limit), in the same palette as figure 3.1.
#   - Node size: the node's strength in that arm (sum of the sizes of its edge weights).
#   - Edge width: the size of the edge weight in that arm. Solid = positive weight (the two nodes respond
#     in the same direction on balance); dashed = negative weight (opposite directions).
#   - Node shape: genes: all squares. Metabolites: shape = RefMet super class.
#     (No significance marks: the bootstrap test of arm differences has been removed for now.)
#   - Labels: the 10 strongest nodes of each layer (genes); every node (metabolites, only 44). In the
#     metabolite figure each connected group also gets its RefMet super-class name above it (every group
#     is a single class, because metabolite edges require the same class).
#   - Display filter: only nodes with at least one edge are drawn (isolated nodes carry no network
#     information); the numbers of drawn nodes and edges are stated in each panel's side label.
#   Titles are descriptive only; interpretation belongs in the report text.
#
# TECH STACK
#   R 4.4; data.table, igraph (layouts), ggplot2, ggrepel (non-overlapping labels).
#
# INPUTS:  $HACK_OUT: 01_nodes_{EE,RE}.csv, 01b_metab_nodes_{EE,RE}.csv, 03_weighted_edges.csv,
#          06_metabolite_edges.csv, 01c_metabolite_ids.csv
# OUTPUTS: $HACK_FIG (default ~/Desktop/output/hackathon; figures are never written into the repo)
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Fixed random seed so the layout is identical on every run.
SEED <- 20260926
# How many of the strongest genes to label per layer.
NLAB <- 10

# ---- figure style (copied from week 5 figure 3.1 / week_4/R/fig.R theme_pub) ---------------------
# Text, axis, muted and strip colours.
INK <- "#1A1A1A"; HAIR <- "#3A3A3A"; MUT <- "#6B7278"; STRIP_BG <- "#EDF0F2"; STRIP_INK <- "#20262B"
# Node colour scale: violet (down) - white - orange (up), limits -2..2, as in figure 3.1.
sc_fill <- scale_fill_gradient2(low = "#6A3D9A", mid = "white", high = "#E66100", midpoint = 0,
                                limits = c(-2, 2), oob = scales::squish,
                                name = "mean response (scaled logFC)")
# The publication theme: small text, no axes (a network has no meaningful axes), legend at the bottom.
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

# ---- helpers ----------------------------------------------------------------------------------------
# Scale a coordinate to 0..1 (a single value goes to the middle).
norm01 <- function(v) if (diff(range(v)) < 1e-9) rep(0.5, length(v)) else (v - min(v)) / diff(range(v))

# Build the drawing data for one network pair (EE above, RE below).
#   nodes: data.table(node, resp_EE, resp_RE, shape_key); edges: data.table(a, b, w_EE, w_RE, sig_EE, sig_RE)
build <- function(nodes, edges) {
  # the shared network over the connected nodes (same edges in both arms)
  g <- graph_from_data_frame(edges[, .(a, b)], directed = FALSE, vertices = nodes[node %in% c(edges$a, edges$b), .(node)])
  # the single shared layout (fixed seed), used for both layers
  set.seed(SEED); L0 <- layout_with_fr(g)
  # scaled to 0..1
  x0 <- norm01(L0[, 1]); y0 <- norm01(L0[, 2])
  # the same positions in both layers: EE in the top band (0.55-0.91), RE in the bottom band (0.05-0.41);
  # the gap above each band leaves room for the metabolite class labels
  lay <- rbind(data.table(node = V(g)$name, arm = "EE", x = x0, y = 0.55 + 0.36 * y0),
               data.table(node = V(g)$name, arm = "RE", x = x0, y = 0.05 + 0.36 * y0))
  # node strength per arm = sum of the sizes of its edge weights
  both <- rbind(edges[, .(node = a, w_EE, w_RE)], edges[, .(node = b, w_EE, w_RE)])
  st <- both[, .(EE = sum(abs(w_EE)), RE = sum(abs(w_RE))), by = node]
  st <- melt(st, id.vars = "node", variable.name = "arm", value.name = "strength")[, arm := as.character(arm)]
  # attach strength, response (the arm's own column) and shape key to each drawn node
  N <- merge(lay, st, by = c("node", "arm"))
  N <- merge(N, nodes, by = "node")
  N[, resp := fifelse(arm == "EE", resp_EE, resp_RE)]
  # edge segments per arm, with that arm's weight
  E <- rbindlist(lapply(c("EE", "RE"), function(arm) {
    # this arm's node positions. The loop variable is copied to `ar` first: inside data.table's [ ],
    # a bare `arm` would mean the COLUMN arm (always equal to itself), silently selecting both layers.
    ar <- arm
    pa <- lay[arm == ar]
    data.table(arm = arm, w = edges[[paste0("w_", arm)]],
               x = pa$x[match(edges$a, pa$node)], y = pa$y[match(edges$a, pa$node)],
               xend = pa$x[match(edges$b, pa$node)], yend = pa$y[match(edges$b, pa$node)])
  }))
  # sign of each edge weight, for the line type
  E[, direction := factor(fifelse(w >= 0, "same direction (w > 0)", "opposite direction (w < 0)"),
                          levels = c("same direction (w > 0)", "opposite direction (w < 0)"))]
  # the drawing data and the counts shown in the side labels
  # which connected group each node belongs to (used for the metabolite class labels)
  N[, comp := components(g)$membership[node]]
  # the drawing data and the counts shown in the side labels
  list(N = N, E = E, n_nodes = vcount(g), n_edges = ecount(g))
}

# Draw one figure from build()'s output.
draw <- function(G, title, shape_values, shape_name, label_rule, file, width = 10.5, height = 6.2,
                 group_labels = FALSE) {
  # which nodes get a label in each layer
  N <- copy(G$N)
  N[, rk := frank(-strength, ties.method = "first"), by = arm]
  N[, lab := fifelse(label_rule(rk), node, NA_character_)]
  # class labels (metabolites): one per connected group per layer, centred above the group. Every group is a
  # single class because edges require the same class; the label is the shape key (the super class).
  GL <- N[, .(x = mean(range(x)), y = max(y) + 0.06, lab = as.character(shape_key[1])), by = .(arm, comp)]
  # side labels for the two layers, stating what is drawn
  SL <- data.table(x = -0.06, y = c(0.73, 0.23),
                   lab = sprintf(c("endurance vs control\n%d nodes · %d edges", "resistance vs control\n%d nodes · %d edges"),
                                 G$n_nodes, G$n_edges))
  # the plot, layer by layer (back to front)
  p <- ggplot() +
    # dashed rule between the two layers, as in figure 3.1
    geom_hline(yintercept = 0.5, colour = HAIR, linetype = "22", linewidth = 0.25) +
    # the edges: width = size of the weight in that arm; solid / dashed = sign
    geom_segment(data = G$E, aes(x, y, xend = xend, yend = yend, linewidth = abs(w), linetype = direction),
                 colour = "grey45", alpha = 0.55) +
    # the nodes: colour = mean response, size = strength in that arm, shape = shape key
    geom_point(data = N, aes(x, y, fill = resp, size = strength, shape = shape_key), colour = "grey25", stroke = 0.22) +
    # labels that avoid each other and the nodes
    geom_text_repel(data = N[!is.na(lab)], aes(x, y, label = lab), size = 2.3, colour = "grey15",
                    min.segment.length = 0.2, segment.size = 0.15, max.overlaps = Inf, seed = SEED) +
    # the class label above each group (only when requested, i.e. for metabolites)
    (if (group_labels) geom_text(data = GL, aes(x, y, label = lab), size = 2.9, colour = "grey35",
                                 fontface = "bold.italic") else NULL) +
    # the layer names on the left
    geom_text(data = SL, aes(x, y, label = lab), angle = 90, size = 2.6, colour = "grey30", fontface = "bold", lineheight = 0.9) +
    # scales: line types, widths, node sizes, shapes, colours
    scale_linetype_manual(values = c("same direction (w > 0)" = "solid", "opposite direction (w < 0)" = "22"),
                          name = "edge weight sign", drop = FALSE) +
    scale_linewidth(range = c(0.1, 1.2), guide = "none") +
    scale_size(range = c(0.8, 4.5), name = "node strength (sum |w|)") +
    scale_shape_manual(values = shape_values, name = shape_name) +
    sc_fill +
    # room on the left for the layer names
    coord_cartesian(xlim = c(-0.1, 1.02), ylim = c(-0.02, 1.02), clip = "off") +
    # legends: shapes drawn grey so they read as shape only
    guides(shape = guide_legend(override.aes = list(fill = "grey85", size = 2.6), order = 1),
           linetype = guide_legend(order = 2), size = guide_legend(order = 3),
           fill = guide_colourbar(order = 4, barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"))) +
    labs(title = title) + theme_net()
  # save with a white background at 300 dpi
  ggsave(file, p, width = width, height = height, dpi = 300, bg = "white")
  message("-> ", file)
}

# ---- figure 10a: gene networks -------------------------------------------------------------------
# Node vectors per arm (16 observed dimensions; the two empty adipose protein columns are ignored).
ge <- fread(file.path(OUT, "01_nodes_EE.csv")); gr <- fread(file.path(OUT, "01_nodes_RE.csv"))
# Each gene's mean response across its dimensions, per arm (missing columns skipped).
gn <- data.table(node = ge$gene_symbol,
                 resp_EE = rowMeans(as.matrix(ge[, -(1:2)]), na.rm = TRUE),
                 resp_RE = rowMeans(as.matrix(gr[, -(1:2)]), na.rm = TRUE))
# Shape key: every gene is drawn as a square.
gn[, shape_key := factor("gene")]
# Weighted gene edges (step 3), named by gene symbol.
gw <- fread(file.path(OUT, "03_weighted_edges.csv"))[, .(a = symbol_a, b = symbol_b, w_EE, w_RE, sig_EE, sig_RE)]
# Build and draw.
G1 <- build(gn, gw)
# Draw and save the gene figure.
draw(G1, "Gene networks: endurance vs resistance (471 genes, STRING combined score >= 700)",
     shape_values = c(gene = 22),
     shape_name = "node", label_rule = function(rk) rk <= NLAB,
     file = file.path(FIG, "10a_gene_networks_EE_vs_RE.png"))

# ---- figure 10b: metabolite networks ---------------------------------------------------------------
# Metabolite vectors per arm (9 dimensions).
me <- fread(file.path(OUT, "01b_metab_nodes_EE.csv")); mr <- fread(file.path(OUT, "01b_metab_nodes_RE.csv"))
# Each metabolite's mean response, per arm.
mn <- data.table(node = me$metabolite, resp_EE = rowMeans(as.matrix(me[, -1])), resp_RE = rowMeans(as.matrix(mr[, -1])))
# Shape key: the RefMet super class (from step 1c).
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
# (each metabolite's super class becomes its shape)
mn[, shape_key := factor(ids$super_class[match(node, ids$metabolite)])]
# Weighted metabolite edges (step 6).
mw <- fread(file.path(OUT, "06_metabolite_edges.csv"))[, .(a = metabolite_a, b = metabolite_b, w_EE, w_RE, sig_EE, sig_RE)]
# Build; the super classes that actually appear get distinct fillable shapes.
G2 <- build(mn, mw)
# (the super classes present among the drawn metabolites)
cls <- sort(unique(as.character(G2$N$shape_key)))
# Draw and save the metabolite figure.
draw(G2, "Metabolite networks: endurance vs resistance (shared Rhea enzyme among the 471 genes + same RefMet super class)",
     shape_values = setNames(c(21, 22, 24, 23, 25)[seq_along(cls)], cls), shape_name = "RefMet super class",
     label_rule = function(rk) rep(TRUE, length(rk)), group_labels = TRUE,
     file = file.path(FIG, "10b_metabolite_networks_EE_vs_RE.png"))
