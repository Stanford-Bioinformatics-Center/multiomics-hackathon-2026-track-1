#!/usr/bin/env Rscript
# =====================================================================================================
# 15_joint_network_classes.R — STEP 15: THE JOINT NETWORK WITH METABOLITE CLASSES SHOWN (TWO OPTIONS)
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   In the step 14 joint network every metabolite is a plain triangle, so it is hard to see which kinds
#   of metabolites (organic acids, nucleic acids, fatty acyls, ...) sit where. This step draws the same
#   network two ways so the team can choose how to show metabolite classes.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Nothing is recomputed: nodes, edges and weights are read from step 14. The metabolite class is the
#   RefMet super class (step 1c). Proteins stay circles, metabolites triangles. Two options, each drawn
#   in the style of 14a (endurance layer above resistance layer) and of 14b (edges = w_EE - w_RE):
#     option 1, "group + name" (15a): a new layout in which metabolites of the same class are pulled
#       together (extra links between same-class metabolites, used for the layout ONLY, never drawn and
#       never weighted by the data); each class gets a faint outline and its name above it.
#     option 2, "colour + legend" (15b): the unchanged step 14 layout; each metabolite triangle is filled
#       with its class colour (proteins keep their 14a / 14b fill), and the class legend sits on the
#       right-hand side (the other legends stay at the bottom).
#   Figures:
#     15a_joint_classes_grouped_EE_vs_RE.png          15a_joint_classes_grouped_edge_difference.png
#     15b_joint_classes_coloured_EE_vs_RE.png         15b_joint_classes_coloured_edge_difference.png
#
# HOW TO RUN
#   After step 14:   Rscript network/15_joint_network_classes.R   (about 30 seconds)
#
# DATA AND PROVENANCE
#   Step 14 tables (joint network built from MoTrPAC human pre-suspension results, STRING >= 700 and
#   Rhea release 142); classes from RefMet (step 1c). See steps 1-6 and 14 for upstream QC.
#
# TECH STACK
#   R 4.4; data.table, igraph (layout), ggplot2, ggrepel (labels), ggforce (class outlines),
#   ggnewscale (a second fill scale, so proteins keep the response colours while metabolites get class colours).
#
# INPUTS (files and columns)
#   $HACK_OUT/14_joint_edges.csv   node_a, node_b, edge_type, w_EE, w_RE, w_diff
#   $HACK_OUT/14_joint_nodes.csv   node, node_type, class, resp_EE, resp_RE, strength_EE, strength_RE, x, y
#
# OUTPUTS (files, locations)
#   $HACK_FIG/15a_*.png and 15b_*.png (four figures, listed above; never in the repo)
#   $HACK_OUT/15_class_layout.csv  node, node_type, class, x, y of the grouped layout (option 1)
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   364 nodes (60 metabolites in 10 super classes: organic acids 19, nucleic acids 16, fatty acyls 11,
#   sphingolipids 6, carbohydrates 2, glycerophospholipids 2, four classes with 1), 764 edges, as in step
#   14. The script stops if an edge end is not a node or a metabolite has no class. 99_validate_outputs.R
#   checks the grouped layout covers exactly the step 14 nodes.
#
# KNOWN LIMITS
#   Option 1 changes node positions, so 15a cannot be overlaid on 14a / 15b; the extra same-class links
#   shape the picture and carry no evidence. Single-metabolite classes get a small outline of their own.
#   Ten class colours are near the limit of what can be told apart; the palette avoids the response
#   (violet / orange) and difference (red / blue) colours, so it cannot be a standard qualitative palette.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel); library(ggforce); library(ggnewscale) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Fixed seed so layouts are identical on every run (same as step 14).
SEED <- 20260926
# How many of the strongest proteins and metabolites to label in each layer (as in step 14).
NLAB <- 10
# Weight of the layout-only links between same-class metabolites (real edges have weight 1).
CLASS_PULL <- 3
# Text colours (as in steps 10, 11, 14).
INK <- "#1A1A1A"; HAIR <- "#3A3A3A"
# Edge colours and line types by type (as in step 14).
TYPE_COL <- c("protein - protein" = "grey55", "metabolite - metabolite" = "#1B7837", "metabolite - protein" = "#8C510A")
TYPE_LTY <- c("protein - protein" = "solid", "metabolite - metabolite" = "42", "metabolite - protein" = "11")
# Node shapes: circles for proteins, triangles for metabolites.
SHAPES <- c(protein = 21, metabolite = 24)
# Difference colours (as in steps 11 and 14).
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"
# Ten distinguishable class colours, chosen to avoid violet and orange (the response fill) and red and blue
# (the difference edges), in order of class size.
PAL10 <- c("#1B9E77", "#222222", "#E6AB02", "#A6761D", "#E7298A", "#66A61E", "#17BECF", "#999999", "#BCBD22", "#8DD3C7")

# ---- read step 14 --------------------------------------------------------------------------------------
# Edges with their weights and type.
E <- fread(file.path(OUT, "14_joint_edges.csv"))[, .(a = node_a, b = node_b, edge_type, w_EE, w_RE, w_diff)]
# Nodes with type, class, response, strength and the step 14 layout.
nodes <- fread(file.path(OUT, "14_joint_nodes.csv"))
# Safety checks: every edge end is a node, and every metabolite has a class.
stopifnot(all(c(E$a, E$b) %in% nodes$node), !anyNA(nodes[node_type == "metabolite", class]))
# Fixed order of edge types (for legends).
E[, edge_type := factor(edge_type, levels = names(TYPE_COL))]
# Classes ordered from most to fewest metabolites, each with its colour.
cls <- nodes[node_type == "metabolite", .N, by = class][order(-N, class), class]
CLASS_COL <- setNames(PAL10[seq_along(cls)], cls)
# Legend text with the number of metabolites in each class, e.g. "Organic acids (19)".
cls_n <- nodes[node_type == "metabolite", .N, by = class]
CLASS_LAB <- setNames(sprintf("%s (%d)", cls, cls_n$N[match(cls, cls_n$class)]), cls)
# Show the classes.
print(cls_n[order(-N)])

# ---- option 1 layout: same-class metabolites pulled together ------------------------------------------
# Real edges (weight 1) plus layout-only links between every pair of metabolites in the same class.
met <- nodes[node_type == "metabolite", .(node, class)]
pull <- met[met, on = "class", allow.cartesian = TRUE][node < i.node, .(a = node, b = i.node, weight = CLASS_PULL)]
# The layout graph (these extra links are used here and nowhere else).
gl <- graph_from_data_frame(rbind(E[, .(a, b, weight = 1)], pull), directed = FALSE, vertices = nodes[, .(node)])
# Force-directed layout; a higher weight pulls two nodes closer.
set.seed(SEED); L1 <- layout_with_fr(gl, weights = E(gl)$weight)
# Positions scaled to 0..1 and stored as a second layout.
norm01 <- function(v) (v - min(v)) / diff(range(v))
grp <- copy(nodes)[, `:=`(x = norm01(L1[, 1])[match(node, V(gl)$name)], y = norm01(L1[, 2])[match(node, V(gl)$name)])]
# Save the grouped layout.
fwrite(grp[, .(node, node_type, class, x, y)], file.path(OUT, "15_class_layout.csv"))

# ---- shared drawing helpers --------------------------------------------------------------------------
# Figure theme (as in step 14); legends at the bottom (option 2 moves only the class legend to the side).
theme_net <- function(base = 8) {
  theme_classic(base_size = base) %+replace% theme(
    text = element_text(colour = INK, size = base),
    axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    legend.title = element_text(colour = INK, size = base * 0.8), legend.text = element_text(colour = INK, size = base * 0.74),
    legend.background = element_blank(), legend.key = element_blank(),
    plot.title = element_text(colour = INK, face = "bold", size = base * 1.18, hjust = 0, margin = margin(b = base * 0.55)),
    plot.background = element_rect(fill = "white", colour = NA))
}
# Helper: edge segment coordinates for a given node-position table.
seg <- function(pos) { S <- copy(E)
  S[, `:=`(x = pos$x[match(a, pos$node)], y = pos$y[match(a, pos$node)], xend = pos$x[match(b, pos$node)], yend = pos$y[match(b, pos$node)])]; S }
# Helper: labels for the NLAB strongest proteins and metabolites by a strength column.
top_labels <- function(N, col) { N <- copy(N); N[, rk := frank(-get(col), ties.method = "first"), by = node_type]
  N[, lab := fifelse(rk <= NLAB, node, NA_character_)]; N }
# Helper: one class name per class, centred above its members (option 1).
class_names <- function(P) P[node_type == "metabolite", .(x = mean(x), y = max(y) + 0.035), by = class]
# Helper: faint class outlines (option 1); expand keeps single metabolites visible.
hulls <- function(P) geom_mark_hull(data = P[node_type == "metabolite"], aes(x, y, group = class), colour = "grey45",
                                    fill = "grey60", alpha = 0.08, linewidth = 0.25, expand = unit(2.2, "mm"),
                                    radius = unit(2, "mm"), concavity = 3)

# ---- EE-over-RE figure (14a style) for either option ----------------------------------------------------
draw_layers <- function(P, option, file, title) {
  # Two layers with the same positions: EE in the top band, RE in the bottom band.
  lay <- rbind(P[, .(node, node_type, class, arm = "EE", x, y = 0.55 + 0.38 * y, resp = resp_EE, strength = strength_EE)],
               P[, .(node, node_type, class, arm = "RE", x, y = 0.05 + 0.38 * y, resp = resp_RE, strength = strength_RE)])
  # Edges of each layer with that arm's weight and its sign.
  S <- rbindlist(lapply(c("EE", "RE"), function(ar) { s <- seg(lay[lay$arm == ar]); s[, `:=`(arm = ar, w = if (ar == "EE") w_EE else w_RE)]; s }))
  S[, sign := factor(fifelse(w >= 0, "same direction (w > 0)", "opposite direction (w < 0)"), levels = c("same direction (w > 0)", "opposite direction (w < 0)"))]
  # Labels: strongest proteins and metabolites per layer (option 1 labels proteins only; classes are named instead).
  lab <- rbindlist(lapply(c("EE", "RE"), function(ar) top_labels(lay[lay$arm == ar], "strength")))
  if (option == 1) lab[node_type == "metabolite", lab := NA_character_]
  # Symmetric colour limit for the mean response (95th percentile of |mean response|, as in step 14).
  lim_r <- as.numeric(quantile(abs(lay$resp), 0.95))
  # Layer names at the left.
  SL <- data.table(x = -0.06, y = c(0.74, 0.24), lab = sprintf(c("endurance vs control\n%d nodes · %d edges",
                                                                 "resistance vs control\n%d nodes · %d edges"), nrow(P), nrow(E)))
  # Base: rule between layers, (option 1) class outlines, edges.
  p <- ggplot() + geom_hline(yintercept = 0.5, colour = HAIR, linetype = "22", linewidth = 0.25) +
    (if (option == 1) lapply(c("EE", "RE"), function(ar) hulls(lay[lay$arm == ar]))) +
    geom_segment(data = S, aes(x, y, xend = xend, yend = yend, colour = edge_type, linewidth = abs(w), linetype = sign), alpha = 0.6) +
    scale_colour_manual(values = TYPE_COL, name = "edge type") +
    scale_linetype_manual(values = c("same direction (w > 0)" = "solid", "opposite direction (w < 0)" = "22"), name = "edge weight sign", drop = FALSE) +
    scale_linewidth(range = c(0.08, 1.2), guide = "none") +
    # proteins (and, in option 1, metabolites): fill = mean response
    geom_point(data = if (option == 1) lay else lay[node_type == "protein"], aes(x, y, fill = resp, size = strength, shape = node_type),
               colour = "grey25", stroke = 0.22) +
    scale_fill_gradient2(low = "#6A3D9A", mid = "white", high = "#E66100", midpoint = 0, limits = c(-lim_r, lim_r),
                         oob = scales::squish, name = if (option == 1) "mean response (normalised logFC)" else "protein mean response\n(normalised logFC)") +
    scale_size(range = c(0.6, 4), name = "node strength (sum |w|)") +
    scale_shape_manual(values = SHAPES, name = "node", drop = FALSE)
  # Option 2: metabolites filled by class on a second fill scale, legend on the side.
  if (option == 2) p <- p + new_scale_fill() +
    geom_point(data = lay[node_type == "metabolite"], aes(x, y, fill = class, size = strength, shape = node_type), colour = "grey20", stroke = 0.25) +
    scale_fill_manual(values = CLASS_COL, labels = CLASS_LAB, breaks = cls, name = "metabolite class\n(RefMet super class)",
                      guide = guide_legend(order = 6, position = "right", override.aes = list(shape = 24, size = 2.6)))
  # Option 1: class names above each group in both layers.
  if (option == 1) p <- p + geom_text_repel(data = rbindlist(lapply(c("EE", "RE"), function(ar) class_names(lay[lay$arm == ar]))),
                                            aes(x, y, label = class), size = 2.3, colour = "grey30", fontface = "bold.italic",
                                            direction = "both", box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.15,
                                            segment.colour = "grey50", max.overlaps = Inf, seed = SEED)
  # Labels, layer names, legends, frame.
  p <- p + geom_text_repel(data = lab[!is.na(lab)], aes(x, y, label = lab), size = 2.1, colour = "grey15", min.segment.length = 0.2,
                           segment.size = 0.12, max.overlaps = Inf, seed = SEED) +
    geom_text(data = SL, aes(x, y, label = lab), angle = 90, size = 2.6, colour = "grey30", fontface = "bold", lineheight = 0.9) +
    guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey85", size = 2.6)),
           colour = guide_legend(order = 2, override.aes = list(linewidth = 0.8, alpha = 1)), linetype = guide_legend(order = 3),
           size = guide_legend(order = 4)) +
    coord_cartesian(xlim = c(-0.1, 1.02), ylim = c(-0.02, 1.0), clip = "off") + labs(title = title) +
    theme_net()
  # Save.
  ggsave(file.path(FIG, file), p, width = if (option == 1) 11 else 12.5, height = 7.5, dpi = 300, bg = "white")
  message("-> ", file.path(FIG, file))
}

# ---- difference figure (14b style) for either option ----------------------------------------------------
draw_diff <- function(P, option, file, title) {
  # Edges at this layout, smallest differences drawn first; colour limit at the 95th percentile of |w_diff|.
  D <- seg(P); lim_d <- as.numeric(quantile(abs(D$w_diff), 0.95)); D <- D[order(abs(w_diff))]
  # Per-node absolute strength difference and labels for the largest (option 1: proteins only).
  N <- copy(P)[, abs_delta := abs(strength_EE - strength_RE)]; N <- top_labels(N, "abs_delta")
  if (option == 1) N[node_type == "metabolite", lab := NA_character_]
  # Base: (option 1) class outlines, then edges coloured by difference.
  q <- ggplot() + (if (option == 1) hulls(N)) +
    geom_segment(data = D, aes(x, y, xend = xend, yend = yend, colour = w_diff, linewidth = abs(w_diff), linetype = edge_type), lineend = "round") +
    scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-lim_d, lim_d), oob = scales::squish,
                           breaks = c(-lim_d, 0, lim_d), labels = c("higher in resistance\n(w_RE > w_EE)", "same", "higher in endurance\n(w_EE > w_RE)"),
                           name = "edge difference  w_EE − w_RE") +
    scale_linetype_manual(values = TYPE_LTY, name = "edge type") +
    scale_linewidth(range = c(0.1, 1.6), guide = "none") +
    # nodes: grey (option 1) or proteins grey + metabolites by class (option 2)
    geom_point(data = N, aes(x, y, size = abs_delta, shape = node_type, fill = if (option == 1) "all" else fifelse(node_type == "protein", "protein", class)),
               colour = "grey25", stroke = 0.22) +
    scale_fill_manual(values = c(all = "grey80", protein = "grey80", CLASS_COL), breaks = if (option == 2) cls else NULL,
                      labels = CLASS_LAB, name = "metabolite class\n(RefMet super class)",
                      guide = if (option == 2) guide_legend(order = 5, position = "right", override.aes = list(shape = 24, size = 2.6)) else "none") +
    scale_size(range = c(0.6, 4), name = "node strength difference (absolute)") +
    scale_shape_manual(values = SHAPES, name = "node")
  # Option 1: class names above each group.
  if (option == 1) q <- q + geom_text_repel(data = class_names(N), aes(x, y, label = class), size = 2.4, colour = "grey30",
                                            fontface = "bold.italic", box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.15,
                                            segment.colour = "grey50", max.overlaps = Inf, seed = SEED)
  # Labels, legends, frame.
  q <- q + geom_text_repel(data = N[!is.na(lab)], aes(x, y, label = lab), size = 2.1, colour = "grey15", min.segment.length = 0.2,
                           segment.size = 0.12, max.overlaps = Inf, seed = SEED) +
    guides(colour = guide_colourbar(order = 1, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"),
                                    title.position = "top", title.hjust = 0.5),
           linetype = guide_legend(order = 2, override.aes = list(colour = "grey30", linewidth = 0.6)),
           shape = guide_legend(order = 3, override.aes = list(fill = "grey80", size = 2.6)), size = guide_legend(order = 4)) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.05), clip = "off") + labs(title = title) +
    theme_net()
  # Save.
  ggsave(file.path(FIG, file), q, width = if (option == 1) 11 else 12.5, height = 7, dpi = 300, bg = "white")
  message("-> ", file.path(FIG, file))
}

# ---- the four figures ---------------------------------------------------------------------------------
# Option 1: grouped layout, class outlines and names.
draw_layers(grp, 1, "15a_joint_classes_grouped_EE_vs_RE.png",
            "Joint protein-metabolite networks, metabolites grouped by RefMet super class: endurance vs resistance")
draw_diff(grp, 1, "15a_joint_classes_grouped_edge_difference.png",
          "Joint protein-metabolite network, metabolites grouped by RefMet super class: endurance minus resistance edge weights")
# Option 2: step 14 layout, metabolites coloured by class, legend on the side.
draw_layers(nodes, 2, "15b_joint_classes_coloured_EE_vs_RE.png",
            "Joint protein-metabolite networks, metabolites coloured by RefMet super class: endurance vs resistance")
draw_diff(nodes, 2, "15b_joint_classes_coloured_edge_difference.png",
          "Joint protein-metabolite network, metabolites coloured by RefMet super class: endurance minus resistance edge weights")
