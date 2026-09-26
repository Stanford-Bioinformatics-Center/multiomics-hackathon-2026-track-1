#!/usr/bin/env Rscript
# =====================================================================================================
# 14_joint_network.R — STEP 14: ONE JOINT PROTEIN + METABOLITE NETWORK (ENDURANCE AND RESISTANCE)
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Genes/proteins and metabolites have so far been studied in separate networks. This step joins them
#   into one network per exercise arm, so relationships that cross between the two layers (an enzyme and
#   the molecules it acts on) can be seen next to the within-layer relationships, and compared between
#   endurance and resistance.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Nodes: the 471 genes/proteins and the 450 metabolites (only those with at least one edge are drawn).
#   Edges (on/off gates that never use the exercise data), three kinds:
#     protein - protein      : STRING, combined score >= 700 (step 2)
#     metabolite - metabolite: shared or STRING-interacting Rhea enzymes + same RefMet super class (step 6)
#     metabolite - protein   : Rhea (step 5): the metabolite is a substrate or product of a reaction that
#                              the protein (one of our 471 genes) catalyses
#   Edge weights, per arm (dot products of the normalised response vectors, steps 1 / 1b):
#     protein - protein and metabolite - metabolite: as in steps 3 and 6;
#     metabolite - protein: the metabolite's 9-value embedding is DOUBLED to 18 values (each tissue x time
#       value written into both the RNA slot and the protein slot of that tissue x time, in the gene
#       embedding's column order), then ONE dot product is taken with the gene's 18-value embedding:
#         w = metab18 . gene18
#       so the metabolite multiplies RNA and protein with equal weight (the team's design). The two empty
#       gene dimensions (adipose protein at 0.5 h and 24 h) are skipped, leaving 16 terms.
#   Then two figures:
#     14a_joint_network_EE_vs_RE.png        endurance (top) and resistance (bottom) layers, one identical
#                                           layout (like figures 10a / 10b)
#     14b_joint_network_edge_difference.png one network whose edges show w_EE - w_RE (like 11a / 11b)
#
# HOW TO RUN
#   After steps 1, 1b, 2, 3, 5 and 6:   Rscript network/14_joint_network.R   (about 20 seconds)
#
# DATA AND PROVENANCE
#   MoTrPAC human pre-suspension study (MotrpacHumanPreSuspensionAnalysis v0.2.4); curated STRING file
#   (combined score >= 700); Rhea release 142 (2026-09-02). Upstream QC and the normalisation (each ome
#   divided by its maximum |logFC| across tissues, times and both arms) are documented in steps 1 / 1b and
#   the README ("Critical QC step").
#
# TECH STACK
#   R 4.4; data.table, igraph (layout), ggplot2, ggrepel (labels).
#
# INPUTS (files and columns)
#   $HACK_OUT/01_nodes_{EE,RE}.csv         entrez_gene, gene_symbol, 18 normalised columns
#   $HACK_OUT/01b_metab_nodes_{EE,RE}.csv  metabolite, 9 normalised columns
#   $HACK_OUT/03_weighted_edges.csv        protein-protein edges: symbol_a, symbol_b, w_EE, w_RE
#   $HACK_OUT/06_metabolite_edges.csv      metabolite-metabolite edges: metabolite_a, metabolite_b, w_EE, w_RE
#   $HACK_OUT/05_metabolite_protein_links.csv  metabolite-protein links: metabolite, gene_symbol, n_reactions
#   $HACK_OUT/01c_metabolite_ids.csv       metabolite classes
#
# OUTPUTS (files, locations, columns)
#   $HACK_OUT/14_joint_edges.csv    node_a, node_b, edge_type, w_EE, w_RE, w_diff
#   $HACK_OUT/14_joint_nodes.csv    node, node_type, class, degree per edge type, strength per arm
#   $HACK_OUT/14_joint_summary.csv  counts per edge type, correlation between the arms per edge type
#   $HACK_FIG/14a_joint_network_EE_vs_RE.png, 14b_joint_network_edge_difference.png (never in the repo)
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   431 protein-protein, 147 metabolite-metabolite and 186 metabolite-protein edges (764 in total). The
#   script stops if the protein-protein or metabolite-metabolite weights differ from steps 3 / 6, or if a
#   metabolite-protein weight differs from a hand-computed doubled-embedding dot product.
#
# KNOWN LIMITS
#   The three edge types have different numbers of terms in their dot products (16, 9 and 16), so their
#   weights are on somewhat different scales; the figures draw them on one scale and edge colour (14a) /
#   line type (14b) shows the type. Differences between the arms are not tested. Metabolite-protein links
#   cover only the 60 metabolites that Rhea links to our 471 genes.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Fixed seed so the layout is identical on every run.
SEED <- 20260926
# How many of the strongest proteins and metabolites to label in each layer.
NLAB <- 10
# Text and figure-furniture colours (as in steps 10 / 11).
INK <- "#1A1A1A"; HAIR <- "#3A3A3A"
# Edge colours by type in figure 14a (neutral hues that do not clash with the violet-orange node fill).
TYPE_COL <- c("protein - protein" = "grey55", "metabolite - metabolite" = "#1B7837", "metabolite - protein" = "#8C510A")
# Edge line types by type in figure 14b (colour is used there for the difference).
TYPE_LTY <- c("protein - protein" = "solid", "metabolite - metabolite" = "42", "metabolite - protein" = "11")
# Node shapes: circles for proteins, triangles for metabolites (the team's choice).
SHAPES <- c(protein = 21, metabolite = 24)
# Difference colours (as in step 11).
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"

# ---- node vectors ---------------------------------------------------------------------------------
# Helper: read a normalised node table as a matrix (rows = node names).
mat <- function(file, id) { x <- fread(file.path(OUT, file)); keep <- setdiff(names(x), c("entrez_gene", "gene_symbol", "metabolite"))
  m <- as.matrix(x[, ..keep]); rownames(m) <- x[[id]]; m }
# Genes: 18 columns (adipose protein 0.5 / 24 h empty), rows named by gene symbol, both arms.
GE <- mat("01_nodes_EE.csv", "gene_symbol"); GR <- mat("01_nodes_RE.csv", "gene_symbol")
# Metabolites: 9 columns, both arms.
ME <- mat("01b_metab_nodes_EE.csv", "metabolite"); MR <- mat("01b_metab_nodes_RE.csv", "metabolite")
# The 9 tissue x time keys (e.g. "blood_4h") shared by the three layers.
key_m <- sub("_metab_", "_", colnames(ME))

# ---- metabolite - protein weights (metabolite embedding doubled to 18 values) -----------------------------
# For each of the 18 gene columns (e.g. "blood_rna_4h", "blood_prot_4h"), find the metabolite column with the
# same tissue and time ("blood_metab_4h"); each metabolite column is therefore used twice (RNA slot, protein slot).
dbl_idx <- match(sub("_(rna|prot)_", "_", colnames(GE)), key_m)
# Safety check: every gene column found its metabolite column, and each metabolite column is used exactly twice.
stopifnot(!anyNA(dbl_idx), all(tabulate(dbl_idx, length(key_m)) == 2))
# The doubled metabolite embeddings: 18 columns in exactly the gene column order, both arms.
ME18 <- ME[, dbl_idx]; colnames(ME18) <- colnames(GE)
MR18 <- MR[, dbl_idx]; colnames(MR18) <- colnames(GR)
# The Rhea links (step 5): one row per metabolite - gene pair.
mp <- unique(fread(file.path(OUT, "05_metabolite_protein_links.csv"))[, .(a = metabolite, b = gene_symbol, n_reactions)])
# One dot product per pair and arm: the metabolite's 18 values times the gene's 18 values, summed. The two
# empty gene dimensions (adipose protein at 0.5 h and 24 h) are skipped (na.rm), leaving 16 terms.
dot18 <- function(M18, G, m, g) unname(rowSums(M18[m, , drop = FALSE] * G[g, , drop = FALSE], na.rm = TRUE))
mp[, `:=`(edge_type = "metabolite - protein", w_EE = dot18(ME18, GE, a, b), w_RE = dot18(MR18, GR, a, b))]
# Safety check: the first pair's weight, recomputed by hand over the 16 observed dimensions, matches.
obs <- !is.na(GE[mp$b[1], ])
stopifnot(sum(obs) == 16, isTRUE(all.equal(sum(ME18[mp$a[1], obs] * GE[mp$b[1], obs]), mp$w_EE[1])))

# ---- the other two edge types (weights from steps 3 and 6) ---------------------------------------------
# Protein - protein edges with their weights.
pp <- fread(file.path(OUT, "03_weighted_edges.csv"))[, .(a = symbol_a, b = symbol_b, w_EE, w_RE)][, edge_type := "protein - protein"]
# Metabolite - metabolite edges with their weights.
mm <- fread(file.path(OUT, "06_metabolite_edges.csv"))[, .(a = metabolite_a, b = metabolite_b, w_EE, w_RE)][, edge_type := "metabolite - metabolite"]
# Safety check: recomputing a protein-protein weight from the vectors reproduces step 3 (first edge).
gcols <- colnames(GE)[colSums(is.na(GE)) < nrow(GE)]
stopifnot(isTRUE(all.equal(sum(GE[pp$a[1], gcols] * GE[pp$b[1], gcols]), pp$w_EE[1])))
# All edges in one table.
E <- rbind(pp, mm, mp[, .(a, b, w_EE, w_RE, edge_type)])
# The difference between the arms (positive = the endurance weight is higher).
E[, w_diff := w_EE - w_RE]
# Fixed order of edge types (for legends and drawing).
E[, edge_type := factor(edge_type, levels = names(TYPE_COL))]

# ---- nodes and layout --------------------------------------------------------------------------------
# Node table: every node with at least one edge, its type and (for metabolites) its super class.
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
nodes <- data.table(node = unique(c(E$a, E$b)))
nodes[, node_type := fifelse(node %in% rownames(ME), "metabolite", "protein")]
nodes[, class := fifelse(node_type == "metabolite", ids$super_class[match(node, ids$metabolite)], "protein")]
# Safety check: no name is both a gene symbol and a metabolite name.
stopifnot(!any(rownames(ME) %in% rownames(GE)))
# Each node's mean response per arm (across its observed dimensions), for the node colour in 14a.
resp <- function(G, M) c(rowMeans(G, na.rm = TRUE), rowMeans(M))
rE <- resp(GE, ME); rR <- resp(GR, MR)
nodes[, `:=`(resp_EE = rE[node], resp_RE = rR[node])]
# The joint graph (edges in table order) and one layout for everything (fixed seed).
g <- graph_from_data_frame(E[, .(a, b)], directed = FALSE, vertices = nodes[, .(node)])
set.seed(SEED); L0 <- layout_with_fr(g)
# Positions scaled to 0..1.
norm01 <- function(v) (v - min(v)) / diff(range(v))
nodes[, `:=`(x = norm01(L0[, 1])[match(node, V(g)$name)], y = norm01(L0[, 2])[match(node, V(g)$name)])]
# Node strength per arm (sum of |w| over the node's edges) and degree per edge type.
both <- rbind(E[, .(node = a, edge_type, w_EE, w_RE)], E[, .(node = b, edge_type, w_EE, w_RE)])
st <- both[, .(strength_EE = sum(abs(w_EE)), strength_RE = sum(abs(w_RE)), degree = .N,
               degree_pp = sum(edge_type == "protein - protein"), degree_mm = sum(edge_type == "metabolite - metabolite"),
               degree_mp = sum(edge_type == "metabolite - protein")), by = node]
nodes <- st[nodes, on = "node"]

# ---- save the tables -----------------------------------------------------------------------------
fwrite(E[, .(node_a = a, node_b = b, edge_type, w_EE, w_RE, w_diff)], file.path(OUT, "14_joint_edges.csv"))
# Node table.
fwrite(nodes, file.path(OUT, "14_joint_nodes.csv"))
# Summary per edge type: counts, correlation between the arms, sign changes.
summ <- E[, .(edges = .N, cor_w_EE_w_RE = round(cor(w_EE, w_RE), 3), sign_changes = sum(sign(w_EE) != sign(w_RE)),
              median_abs_w = signif(median(abs(c(w_EE, w_RE))), 3)), by = edge_type]
summ <- rbind(summ, data.table(edge_type = "all", edges = nrow(E), cor_w_EE_w_RE = round(cor(E$w_EE, E$w_RE), 3),
                               sign_changes = sum(sign(E$w_EE) != sign(E$w_RE)), median_abs_w = signif(median(abs(c(E$w_EE, E$w_RE))), 3)))
fwrite(summ, file.path(OUT, "14_joint_summary.csv"))
# Show it.
print(summ)

# ---- figure theme ------------------------------------------------------------------------------------
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
# Helper: labels for the NLAB strongest proteins and NLAB strongest metabolites by a strength column.
top_labels <- function(N, col) { N <- copy(N); N[, rk := frank(-get(col), ties.method = "first"), by = node_type]
  N[, lab := fifelse(rk <= NLAB, node, NA_character_)]; N }

# ---- figure 14a: endurance (top) and resistance (bottom), one identical layout ---------------------
# Two layers with the same positions: EE in the top band (0.55-0.95), RE in the bottom band (0.05-0.45).
lay <- rbind(nodes[, .(node, node_type, arm = "EE", x, y = 0.55 + 0.40 * y, resp = resp_EE, strength = strength_EE)],
             nodes[, .(node, node_type, arm = "RE", x, y = 0.05 + 0.40 * y, resp = resp_RE, strength = strength_RE)])
# Edges of each layer with that arm's weight.
S <- rbindlist(lapply(c("EE", "RE"), function(ar) { s <- seg(lay[lay$arm == ar]); s[, `:=`(arm = ar, w = if (ar == "EE") w_EE else w_RE)]; s }))
S[, sign := factor(fifelse(w >= 0, "same direction (w > 0)", "opposite direction (w < 0)"), levels = c("same direction (w > 0)", "opposite direction (w < 0)"))]
# Labels per layer.
lab <- rbindlist(lapply(c("EE", "RE"), function(ar) top_labels(lay[lay$arm == ar], "strength")))
# Colour limits for the mean response (95th percentile of |mean response|, symmetric).
lim_r <- as.numeric(quantile(abs(lay$resp), 0.95))
# Side labels for the two layers.
SL <- data.table(x = -0.06, y = c(0.75, 0.25),
                 lab = sprintf(c("endurance vs control\n%d nodes · %d edges", "resistance vs control\n%d nodes · %d edges"), nrow(nodes), nrow(E)))
# Draw.
p <- ggplot() +
  # dashed rule between the layers
  geom_hline(yintercept = 0.5, colour = HAIR, linetype = "22", linewidth = 0.25) +
  # edges: colour = type, width = |w|, line type = sign
  geom_segment(data = S, aes(x, y, xend = xend, yend = yend, colour = edge_type, linewidth = abs(w), linetype = sign), alpha = 0.6) +
  # nodes: fill = mean response, size = strength, shape = protein (circle) / metabolite (triangle)
  geom_point(data = lay, aes(x, y, fill = resp, size = strength, shape = node_type), colour = "grey25", stroke = 0.22) +
  # labels
  geom_text_repel(data = lab[!is.na(lab)], aes(x, y, label = lab), size = 2.1, colour = "grey15", min.segment.length = 0.2,
                  segment.size = 0.12, max.overlaps = Inf, seed = SEED) +
  # layer names
  geom_text(data = SL, aes(x, y, label = lab), angle = 90, size = 2.6, colour = "grey30", fontface = "bold", lineheight = 0.9) +
  # scales
  scale_colour_manual(values = TYPE_COL, name = "edge type") +
  scale_linetype_manual(values = c("same direction (w > 0)" = "solid", "opposite direction (w < 0)" = "22"), name = "edge weight sign", drop = FALSE) +
  scale_linewidth(range = c(0.08, 1.2), guide = "none") +
  scale_size(range = c(0.6, 4), name = "node strength (sum |w|)") +
  scale_shape_manual(values = SHAPES, name = "node") +
  scale_fill_gradient2(low = "#6A3D9A", mid = "white", high = "#E66100", midpoint = 0, limits = c(-lim_r, lim_r),
                       oob = scales::squish, name = "mean response (normalised logFC)") +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey85", size = 2.6)),
         colour = guide_legend(order = 2, override.aes = list(linewidth = 0.8, alpha = 1)), linetype = guide_legend(order = 3),
         size = guide_legend(order = 4), fill = guide_colourbar(order = 5, barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"))) +
  coord_cartesian(xlim = c(-0.1, 1.02), ylim = c(-0.02, 1.0), clip = "off") +
  labs(title = "Joint protein-metabolite networks: endurance vs resistance (STRING, Rhea and class-rule edges)") + theme_net()
# Save.
ggsave(file.path(FIG, "14a_joint_network_EE_vs_RE.png"), p, width = 11, height = 7, dpi = 300, bg = "white")
message("-> ", file.path(FIG, "14a_joint_network_EE_vs_RE.png"))

# ---- figure 14b: one network, edges = w_EE - w_RE -----------------------------------------------------
# Edges at the single layout.
D <- seg(nodes)
# Symmetric colour limit at the 95th percentile of |w_diff|; smallest differences drawn first.
lim_d <- as.numeric(quantile(abs(D$w_diff), 0.95)); D <- D[order(abs(w_diff))]
# Per-node strength difference (sum |w_EE| - sum |w_RE|) and labels for the largest.
N <- copy(nodes)[, delta := strength_EE - strength_RE][, abs_delta := abs(delta)]
N <- top_labels(N, "abs_delta")
# Draw.
q <- ggplot() +
  # edges: colour and width = difference, line type = edge type
  geom_segment(data = D, aes(x, y, xend = xend, yend = yend, colour = w_diff, linewidth = abs(w_diff), linetype = edge_type), lineend = "round") +
  # nodes: grey, size = |strength difference|, shape = protein (circle) / metabolite (triangle)
  geom_point(data = N, aes(x, y, size = abs_delta, shape = node_type), fill = "grey80", colour = "grey30", stroke = 0.22) +
  # labels
  geom_text_repel(data = N[!is.na(lab)], aes(x, y, label = lab), size = 2.1, colour = "grey15", min.segment.length = 0.2,
                  segment.size = 0.12, max.overlaps = Inf, seed = SEED) +
  # scales
  scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-lim_d, lim_d), oob = scales::squish,
                         breaks = c(-lim_d, 0, lim_d), labels = c("higher in resistance\n(w_RE > w_EE)", "same", "higher in endurance\n(w_EE > w_RE)"),
                         name = "edge difference  w_EE − w_RE") +
  scale_linetype_manual(values = TYPE_LTY, name = "edge type") +
  scale_linewidth(range = c(0.1, 1.6), guide = "none") +
  scale_size(range = c(0.6, 4), name = "node strength difference (absolute)") +
  scale_shape_manual(values = SHAPES, name = "node") +
  guides(colour = guide_colourbar(order = 1, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5),
         linetype = guide_legend(order = 2, override.aes = list(colour = "grey30", linewidth = 0.6)),
         shape = guide_legend(order = 3, override.aes = list(size = 2.6)), size = guide_legend(order = 4)) +
  coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02), clip = "off") +
  labs(title = "Joint protein-metabolite network: endurance minus resistance edge weights") + theme_net()
# Save.
ggsave(file.path(FIG, "14b_joint_network_edge_difference.png"), q, width = 11, height = 6.5, dpi = 300, bg = "white")
message("-> ", file.path(FIG, "14b_joint_network_edge_difference.png"))
