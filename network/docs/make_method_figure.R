#!/usr/bin/env Rscript
# =====================================================================================================
# docs/make_method_figure.R — THE METHOD FIGURE: hard (physical) edges x soft (dot-product) weights
# =====================================================================================================
# PURPOSE: one picture that explains how every edge of our network is made, computed (not drawn by hand)
#   by the exnet engine on its tested toy example, so every number in the figure is the engine's output.
#   a  HARD layer: physical-interaction databases (STRING >= 700 for protein pairs, Rhea for enzyme-metabolite
#      pairs) decide WHETHER an edge exists. A and E respond almost identically but have no physical link, so
#      they never get an edge.
#   b  SOFT layer: each molecule's normalised exercise response, per arm (one value per tissue x ome x time);
#      the metabolite's values are doubled into the RNA and protein slots of the same tissue.
#   c  The weight of one edge = the dot product of its two vectors, worked column by column for A-B.
#   d  The result: one network per arm (edge width = |w|, dashed = negative) and the arm-specific edges
#      (strong in one arm only; threshold 0.3 here).
# HOW TO RUN:  Rscript network/docs/make_method_figure.R   (seconds; run_all.sh step M)
# OUTPUTS: network/docs/figures/fig_method_hard_soft_edges.png (committed, shown in the README) and .pdf;
#          a copy in $HACK_FIG.
# =====================================================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
# Folders: this script's own folder (docs/), the network folder above it, the figure output folder.
args <- commandArgs(trailingOnly = FALSE); DOCS <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
HERE <- dirname(DOCS); FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# The engine (loaded from source) and the shared figure style.
suppressMessages(pkgload::load_all(file.path(HERE, "engine"), quiet = TRUE, export_all = FALSE))
source(file.path(HERE, "R", "figure_style.R"))

# ---- the toy network, built by the engine ----------------------------------------------------------------
t <- toy_example()
pp <- build_network(t$hard_pp, t$gene_EE, t$gene_RE, "protein - protein")
dE <- bind_embeddings(double_embedding(t$metab_EE, t$gene_EE), t$gene_EE); dR <- bind_embeddings(double_embedding(t$metab_RE, t$gene_RE), t$gene_RE)
net <- combine_networks(pp, build_network(t$hard_mp, dE, dR, "metabolite - protein"))
E <- edge_table(net); S <- arm_specific_edges(net, tau = 0.3)
# Fixed positions for the six molecules (so all network panels line up).
POS <- data.table(node = c("A", "B", "C", "D", "E", "M"), x = c(0, 1.1, 1.0, 2.1, -0.9, 2.9), y = c(1.0, 1.55, 0.25, 0.55, 1.75, 1.45),
                  type = c(rep("protein", 5), "metabolite"))
seg <- function(e) merge(merge(e, POS[, .(a = node, x1 = x, y1 = y)], by = "a"), POS[, .(b = node, x2 = x, y2 = y)], by = "b")
nodes_layer <- function() list(geom_point(data = POS, aes(x, y, shape = type), size = 6, fill = "white", colour = "black", stroke = 0.6),
                               geom_text(data = POS, aes(x, y, label = node), size = 2.8, fontface = "bold"),
                               scale_shape_manual(values = c(protein = 21, metabolite = 24), guide = "none"))
LIM <- list(xlim(-1.3, 3.3), ylim(-0.05, 2.05))

# ---- a: the hard layer ---------------------------------------------------------------------------------------
ea <- seg(E[, .(a, b, source)])
pa <- ggplot() +
  geom_segment(data = ea, aes(x1, y1, xend = x2, yend = y2, linetype = source), colour = "grey30", linewidth = 0.6) +
  geom_segment(aes(x = -0.9, y = 1.75, xend = 0, yend = 1.0), colour = "#D7301F", linewidth = 0.5, linetype = "22") +
  annotate("point", x = -0.45, y = 1.38, shape = 4, size = 3.2, stroke = 1.3, colour = "#D7301F") +   # a cross (a plotted symbol, so every device draws it)
  annotate("text", x = -1.28, y = 0.38, hjust = 0, size = 2.3, colour = "#D7301F", label = "A and E respond alike\n(dot product 0.57)\nbut no physical link:\nNO edge") +
  nodes_layer() + scale_linetype_manual(values = c(STRING = "solid", Rhea = "11"), name = "physical database") + LIM +
  labs(title = "Hard layer: does an edge exist?", subtitle = "STRING protein-protein (score >= 700) or Rhea enzyme-metabolite links") +
  theme_motrpac_void() + theme(legend.position = "bottom")

# ---- b: the soft layer ---------------------------------------------------------------------------------------
hm <- rbind(data.table(arm = "endurance", node = rownames(dE@values), dE@values), data.table(arm = "resistance", node = rownames(dR@values), dR@values))
hm <- melt(hm, id.vars = c("arm", "node"), variable.name = "dim", value.name = "z")
hm[, dim := factor(gsub("_", " ", sub("_4h$", "", dim)), levels = c("muscle rna", "muscle prot", "blood rna", "blood prot"))]
hm[, node := factor(ifelse(node == "M", "M (doubled)", node), levels = rev(c("A", "B", "C", "D", "E", "M (doubled)")))]
pb <- ggplot(hm, aes(dim, node, fill = z)) + geom_tile(colour = "white", linewidth = 0.6) + geom_text(aes(label = sprintf("%.1f", z)), size = 2.1) +
  facet_wrap(~arm) + scale_fill_gradient2(low = DIV_PAL[["low"]], mid = DIV_PAL[["mid"]], high = DIV_PAL[["high"]], limits = c(-0.6, 0.6), name = "response z") +
  labs(title = "Soft layer: how does each molecule respond?", subtitle = "normalised log fold change per tissue x ome (x time in the real data: 16 values)", x = NULL, y = NULL) +
  theme_motrpac() + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_blank())

# ---- c: one worked dot product (A-B) ------------------------------------------------------------------------------
wk <- rbindlist(lapply(c("EE", "RE"), function(a) { m <- if (a == "EE") t$gene_EE@values else t$gene_RE@values
  data.table(arm = c(EE = "endurance", RE = "resistance")[[a]], dim = factor(gsub("_", " ", sub("_4h$", "", colnames(m))), levels = levels(hm$dim)),
             zA = m["A", ], zB = m["B", ], prod = m["A", ] * m["B", ]) }))
wkl <- melt(wk, id.vars = c("arm", "dim"), variable.name = "col", value.name = "v")[, col := factor(col, levels = c("zA", "zB", "prod"), labels = c("z(A)", "z(B)", "z(A) x z(B)"))]
tot <- wk[, .(w = sum(prod)), by = arm]
pc <- ggplot(wkl, aes(col, dim)) + geom_tile(aes(fill = v), colour = "white", linewidth = 0.6) + geom_text(aes(label = sprintf("%.2f", round(v, 2) + 0)), size = 2.2) +   # (+ 0 turns -0.00 into 0.00)
  geom_text(data = tot, aes(x = 3, y = 0.25, label = sprintf("sum = w = %.2f", w)), inherit.aes = FALSE, size = 2.5, fontface = "bold") +
  facet_wrap(~arm) + scale_y_discrete(limits = rev, expand = expansion(add = c(0.9, 0.3))) +
  scale_fill_gradient2(low = DIV_PAL[["low"]], mid = DIV_PAL[["mid"]], high = DIV_PAL[["high"]], limits = c(-0.6, 0.6), guide = "none") +
  labs(title = "Weight of edge A-B = dot product of the two vectors", subtitle = "w(A,B) = sum over dimensions of z(A) x z(B), computed separately in each arm", x = NULL, y = NULL) +
  theme_motrpac() + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_blank())

# ---- d: the weighted networks and arm-specific edges ---------------------------------------------------------------
wnet <- function(arm) { e <- seg(E[, .(a, b, w = if (arm == "EE") w_EE else w_RE)])
  ggplot() + geom_segment(data = e, aes(x1, y1, xend = x2, yend = y2, linewidth = abs(w), colour = w, linetype = w < 0)) +
    geom_label(data = e, aes((x1 + x2) / 2, (y1 + y2) / 2, label = sprintf("%.2f", w)), size = 2, label.size = 0, label.padding = grid::unit(0.6, "mm")) +
    nodes_layer() + LIM + scale_linewidth(range = c(0.3, 3), limits = c(0, 0.5), guide = "none") +
    scale_colour_gradient2(low = DIV_PAL[["low"]], mid = "grey70", high = DIV_PAL[["high"]], limits = c(-0.5, 0.5), guide = "none") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    labs(title = c(EE = "endurance: w_EE", RE = "resistance: w_RE")[[arm]]) + theme_motrpac_void() }
es <- seg(S[, .(a, b, specificity)])[, specificity := factor(specificity, levels = c("endurance-specific", "resistance-specific", "both", "neither"))]
pd3 <- ggplot() + geom_segment(data = es, aes(x1, y1, xend = x2, yend = y2, colour = specificity), linewidth = 1.4) + nodes_layer() + LIM +
  scale_colour_manual(values = c(`endurance-specific` = ARM_COL[["endurance"]], `resistance-specific` = ARM_COL[["resistance"]], both = ARM_COL[["both"]], neither = ARM_COL[["neither"]]),
                      name = NULL, drop = FALSE) +
  labs(title = "arm-specific (w >= 0.3)") + theme_motrpac_void() + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))

# ---- assemble --------------------------------------------------------------------------------------------------------
# The three small networks form ONE panel (d), with a shared title.
pd <- wrap_elements(full = (wnet("EE") | wnet("RE") | pd3) + plot_annotation(title = "    The result: one weighted network per arm, and its arm-specific edges",
                          theme = theme(plot.title = element_text(face = "bold", size = 9, family = "Helvetica"))))
f <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  title = "Fig. M | How an edge is made: a physical link decides WHETHER, the exercise responses decide HOW STRONG",
                  subtitle = "w_arm(u, v) = 1[u-v linked in STRING >= 700 or Rhea] x sum_d z_u,d(arm) z_v,d(arm)    (toy example; the real vectors have 16 tissue x ome x time values)",
                  theme = theme(plot.title = element_text(face = "bold", size = 10, family = "Helvetica"), plot.subtitle = element_text(size = 8, family = "Helvetica", colour = "grey25")))
dir.create(file.path(DOCS, "figures"), showWarnings = FALSE)
save_figure(f, file.path(DOCS, "figures", "fig_method_hard_soft_edges"), width = 250, height = 150)
file.copy(file.path(DOCS, "figures", "fig_method_hard_soft_edges.png"), file.path(FIG, "fig_method_hard_soft_edges.png"), overwrite = TRUE)
