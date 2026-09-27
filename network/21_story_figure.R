#!/usr/bin/env Rscript
# =====================================================================================================
# 21_story_figure.R — STEP 21 (FIGURE 1): THE WHOLE STORY ON ONE PAGE
# =====================================================================================================
#
# PURPOSE
#   One summary figure, in the style of the MoTrPAC landscape papers, that carries the project's three-part
#   story: (1) DISCOVERY — in muscle, endurance (not resistance) moves the proteins altered in type 2 diabetes
#   (T2D) back toward healthy, and those proteins form a connected subgraph of our physically-gated,
#   response-weighted network; (2) REPLICATION — the same endurance-minus-resistance gap in independent T2D
#   cohorts; (3) THE BLOOD SIDE — in plasma (UK Biobank), it is RESISTANCE that moves the proteins linked to
#   ageing and to future T2D back toward the young / low-risk state.
#
# WHAT THIS SCRIPT DOES (plain language)
#   It draws nothing new statistically; every number comes from steps 14, 19 and 20:
#   a  muscle, T2D-altered proteins (Öhman 2021): T2D z vs mean normalised muscle response, per arm, with the
#      reversal statistic and its permutation p (step 20 set "ohman_2021").
#   b  the connector subgraph of those proteins (step 19 story 2: T2D proteins + nodes linked to >= 2 of them),
#      its edges coloured by arm specificity (red = strong after endurance only, blue = resistance only; the
#      engine's rule, top 25% of |w| over both arms), node outline = direction in T2D.
#   c  every tissue-matched disease test of step 20 in one forest (reversal per arm and endurance - resistance),
#      grouped by the tissue of the readout: muscle vs blood.
#   d, e  UK Biobank plasma: proteins associated with age (Sun 2023) and with future T2D (Gadd 2024) vs the mean
#      normalised blood-protein response (same Olink platform), per arm.
#
# HOW TO RUN
#   After steps 19 and 20:   Rscript network/21_story_figure.R   (seconds; run_all.sh step 21)
#
# INPUTS:  $HACK_OUT/20_tests.csv, 20_disease_scores.csv, 19_t2d_story_nodes.csv, 14_joint_edges.csv,
#          01_nodes_EE.csv, 01_nodes_RE.csv; network/engine (exnet), network/R/figure_style.R
# OUTPUTS: $HACK_FIG/fig1_story.png / .pdf; network/docs/figures/fig1_story.png (committed, shown in the README)
#
# KNOWN LIMITS
#   Direction matches between acute exercise in healthy adults and disease (or age) are not evidence of treatment;
#   the blood result compares plasma disease / age associations with the acute exercise response in blood.
# =====================================================================================================
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
# Folders: the network folder (this script's folder), outputs and figures.
args <- commandArgs(trailingOnly = FALSE); HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# The engine (for the arm-specific rule) and the shared figure style.
suppressMessages(pkgload::load_all(file.path(HERE, "engine"), quiet = TRUE, export_all = FALSE))
source(file.path(HERE, "R", "figure_style.R"))
SEED <- 20260926; ALPHA <- 0.05

# ---- inputs ----------------------------------------------------------------------------------------------
TS <- fread(file.path(OUT, "20_tests.csv")); DS <- fread(file.path(OUT, "20_disease_scores.csv"))
EEg <- fread(file.path(OUT, "01_nodes_EE.csv")); REg <- fread(file.path(OUT, "01_nodes_RE.csv"))
# Mean normalised response per gene over one tissue's cells (the readout used by step 20).
tmean <- function(M, prefix) { cols <- grep(paste0("^", prefix), names(M), value = TRUE); v <- rowMeans(as.matrix(M[, ..cols]), na.rm = TRUE); v[is.nan(v)] <- NA; setNames(v, M$gene_symbol) }
R <- list(muscle = list(EE = tmean(EEg, "muscle_"), RE = tmean(REg, "muscle_")), blood_prot = list(EE = tmean(EEg, "blood_prot_"), RE = tmean(REg, "blood_prot_")))

# ---- a, d, e: disease z vs exercise response, per arm -----------------------------------------------------------
scatter <- function(set_name, tissue, alpha, xlab, ylab, title) {
  d <- DS[set == set_name & gene %in% EEg$gene_symbol & p < alpha]; t <- TS[set == set_name]
  x <- rbind(d[, .(gene, z, r = R[[tissue]]$EE[gene], arm = "endurance")], d[, .(gene, z, r = R[[tissue]]$RE[gene], arm = "resistance")])[is.finite(r)]
  fp <- function(p) if (p < 1e-3) "p < 0.001" else sprintf("p = %.2g", p)          # permutation p (10,000 draws: smallest possible 1e-4)
  ann <- data.table(arm = c("endurance", "resistance"), lab = c(sprintf("reversal %.2f\n%s", t$rev_EE, fp(t$p_EE)), sprintf("reversal %.2f\n%s", t$rev_RE, fp(t$p_RE))))
  ggplot(x, aes(z, r)) + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
    geom_point(aes(colour = z < 0), size = 0.9, alpha = 0.8) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
    scale_colour_manual(values = c(`TRUE` = DIR_COL[["lower"]], `FALSE` = DIR_COL[["higher"]]), guide = "none") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.15, size = 2.3, inherit.aes = FALSE) +
    facet_wrap(~arm) + labs(x = xlab, y = ylab, title = title) + theme_motrpac(7) }
pa <- scatter("ohman_2021", "muscle", ALPHA, "T2D z (Öhman 2021; < 0 = lower in T2D)", "muscle response", "Muscle: 68 proteins altered in T2D")
pd <- scatter("sun_2023_ukb_age", "blood_prot", 1.7e-5, "age z (UK Biobank; < 0 = lower with age)", "blood protein response", "Blood: 293 plasma proteins associated with age")
pe <- scatter("gadd_2024_ukb_incident_T2D", "blood_prot", 3.1e-6, "future-T2D z (UK Biobank; < 0 = lower risk)", "blood protein response", "Blood: 297 plasma proteins linked to future T2D")

# ---- b: the connector subgraph with arm-specific edges ----------------------------------------------------------------
# Rebuild the whole joint network with the engine (it reproduces step 14 exactly; tested), so the arm-specific
# threshold is the same one the step 17 pages use: top 25% of |w| over all edges and both arms.
rd <- function(f, id) { x <- fread(file.path(OUT, f)); m <- as.matrix(x[, setdiff(names(x), c("entrez_gene", "gene_symbol", "metabolite")), with = FALSE]); rownames(m) <- x[[id]]; m }
G <- list(EE = rd("01_nodes_EE.csv", "gene_symbol"), RE = rd("01_nodes_RE.csv", "gene_symbol"))
M <- list(EE = rd("01b_metab_nodes_EE.csv", "metabolite"), RE = rd("01b_metab_nodes_RE.csv", "metabolite"))
J <- fread(file.path(OUT, "14_joint_edges.csv"))
layer <- function(type, EE, RE) { j <- J[edge_type == type]; build_network(physical_edges(j$node_a, j$node_b, type), EE, RE, type) }
gE <- embedding(G$EE, "EE"); gR <- embedding(G$RE, "RE")
net <- combine_networks(layer("protein - protein", gE, gR),
                        layer("metabolite - metabolite", embedding(M$EE, "EE", "metabolite"), embedding(M$RE, "RE", "metabolite")),
                        layer("metabolite - protein", bind_embeddings(double_embedding(embedding(M$EE, "EE", "metabolite"), gE), gE),
                                                      bind_embeddings(double_embedding(embedding(M$RE, "RE", "metabolite"), gR), gR)))
SP <- arm_specific_edges(net, q = 0.75)
# The connector subgraph (step 19, story 2) and a fixed-seed layout of it.
SN <- fread(file.path(OUT, "19_t2d_story_nodes.csv"))[story == "2 connector subgraph"]
es <- SP[a %in% SN$node & b %in% SN$node][, specificity := factor(specificity, levels = c("endurance-specific", "resistance-specific", "both", "neither"))]
g <- igraph::graph_from_data_frame(es[, .(a, b)], directed = FALSE, vertices = SN[, .(node)])
set.seed(SEED); L <- igraph::layout_with_kk(g); lay <- data.table(node = igraph::V(g)$name, x = L[, 1], y = L[, 2])
es <- merge(merge(es, lay[, .(a = node, x1 = x, y1 = y)], by = "a"), lay[, .(b = node, x2 = x, y2 = y)], by = "b")
nd <- merge(lay, SN[, .(node, role, node_type, ohman_z)], by = "node")
nd[, outline := fifelse(role == "connector", "#222222", fifelse(ohman_z < 0, DIR_COL[["lower"]], DIR_COL[["higher"]]))]
pb <- ggplot() + geom_segment(data = es[specificity %in% c("both", "neither")], aes(x1, y1, xend = x2, yend = y2), colour = "grey80", linewidth = 0.35) +
  geom_segment(data = es[specificity %in% c("endurance-specific", "resistance-specific")], aes(x1, y1, xend = x2, yend = y2, colour = specificity), linewidth = 1.1) +
  geom_point(data = nd, aes(x, y, shape = node_type), fill = "white", colour = nd$outline, size = 2.2, stroke = 0.9) +
  geom_text_repel(data = nd, aes(x, y, label = node), size = 1.9, seed = SEED, max.time = 60, max.iter = 1e4, box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.15) +
  scale_colour_manual(values = c(`endurance-specific` = ARM_COL[["endurance"]], `resistance-specific` = ARM_COL[["resistance"]]),
                      labels = c("strong after endurance only", "strong after resistance only"), name = NULL) +
  scale_shape_manual(values = c(protein = 21, metabolite = 24), guide = "none") +
  labs(title = sprintf("T2D proteins + connectors: %d-node subgraph (connectivity p 0.039)", nrow(SN))) +
  theme_motrpac_void(7) + theme(legend.position = "bottom")

# ---- c: every tissue-matched test, grouped by readout tissue -------------------------------------------------------------
tl <- c(muscle = "muscle", adipose = "adipose", blood_prot = "blood (plasma protein)")
fr <- TS[!is.na(p_diff) & primary == TRUE][, grp := factor(tl[tissue], levels = tl)]
fr[, lab := sub(" \\(post hoc\\)", "*", label)]
fr <- rbind(fr[, .(grp, lab, arm = "endurance", v = rev_EE, p = p_EE)], fr[, .(grp, lab, arm = "resistance", v = rev_RE, p = p_RE)], fr[, .(grp, lab, arm = "endurance - resistance", v = diff, p = p_diff)])
fr[, arm := factor(arm, levels = c("endurance", "resistance", "endurance - resistance"))]
fr[, sig := p < ALPHA]
pc <- ggplot(fr, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(shape = sig), size = 1.8, position = position_dodge(0.65), stroke = 0.6) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p < 0.05", `FALSE` = "p >= 0.05"), name = NULL) +
  scale_colour_manual(values = c(endurance = ARM_COL[["endurance"]], resistance = ARM_COL[["resistance"]], `endurance - resistance` = "black"), name = NULL) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = "reversal (> 0: exercise moves the altered proteins back toward healthy / young)", y = NULL,
       title = "Every tissue-matched disease test (step 20; * = pooled cohorts, post hoc)") +
  theme_motrpac(7) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0), legend.box = "vertical")

# ---- assemble ------------------------------------------------------------------------------------------------------------------
# free() stops the long labels of the forest (panel e) from pushing the panels above to the right
f <- (pa | pb) / (pd | pe) / free(pc) + plot_layout(heights = c(1, 0.8, 1.05)) +
  plot_annotation(tag_levels = "a", title = "Fig. 1 | Endurance and resistance move disease-linked proteins back in different tissues",
                  subtitle = "Physically gated, response-weighted multi-omic network (MoTrPAC acute exercise, adults) tested against human T2D, insulin-resistance and ageing proteomes",
                  theme = theme(plot.title = element_text(face = "bold", size = 10, family = "Helvetica"), plot.subtitle = element_text(size = 7.5, family = "Helvetica", colour = "grey25")))
save_figure(f, file.path(FIG, "fig1_story"), width = 220, height = 260)
dir.create(file.path(HERE, "docs", "figures"), showWarnings = FALSE)
file.copy(file.path(FIG, "fig1_story.png"), file.path(HERE, "docs", "figures", "fig1_story.png"), overwrite = TRUE)
