#!/usr/bin/env Rscript
# =====================================================================================================
# 19_t2d_stories.R — STEP 19 (FIGURES 19a / 19b / 19c): ENDURANCE VS RESISTANCE IN THE CONTEXT OF T2D
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Our key question is how endurance (EE) and resistance (RE) exercise differ in the context of preventing
#   disease. The disease data come from Amar et al. 2024 (Cell Metab 36:1411; doi 10.1016/j.cmet.2023.12.021).
#   Of its 9 disease sets only the two type 2 diabetes (T2D) sets are from a tissue we measure (skeletal muscle;
#   the others are heart or liver). We keep ONLY the strongest T2D signals: proteins validated in BOTH T2D
#   datasets in the same direction (the consensus set of step 17: Öhman 2021 p < 0.05 and listed by Chae 2018,
#   same sign). Three stories, each comparing EE with RE, each with a permutation test:
#     STORY 1 — small subgraphs: each consensus protein with its direct network neighbours.
#     STORY 2 — larger subgraph: the smallest tree in the network that connects the consensus proteins.
#     STORY 3 — protein by protein: the consensus proteins' muscle response per cell and arm, with the set-level test
#               per tissue (muscle primary; blood / adipose secondary).
#
# WHAT THIS SCRIPT DOES (plain language)
#   T2D direction per protein: Öhman 2021 muscle proteome (full table) and Chae 2018 (published significant
#   proteins, FDR < 0.1), as processed in the Venus project week 6; signed z from the p-value (sign from the log2
#   change), exactly as step 18. Consensus = Öhman p < 0.05 AND listed by Chae AND the same sign: 85 proteins
#   genome-wide, 8 among our 471, 6 of them nodes of the joint network.
#   Exercise response per node: the mean of its normalised log fold changes (the numbers the embeddings and edge
#   weights are built from) over the cells of one tissue. MUSCLE is primary (the T2D data are muscle).
#   REVERSAL of a protein (per arm) = -sign(T2D change) x response: positive = exercise moves it opposite to T2D
#   (up if lower in T2D, down if higher). Set statistic = the mean over the consensus proteins.
#   Tests (10,000 random draws unless exact; seed 20260926):
#     - reversal per arm and EE - RE over the consensus proteins: exact sign-flip test (all 2^8 patterns).
#     - story 1: mean w_EE - w_RE over the edges touching the consensus proteins vs the same for random protein
#       sets with the same degree make-up (5 degree bins, drawn from the network proteins Öhman measured).
#     - story 2: Steiner tree (Kou-Markowsky-Berman: minimum spanning tree over the terminals' shortest-path
#       distances, paths expanded, spanning tree, non-terminal leaves pruned). Its size (are the consensus proteins
#       closer than chance?) and mean w_EE - w_RE over its edges vs trees built the same way on degree-matched
#       random terminal sets from the same connected part of the network.
#     - sensitivity (single-study sets, NOT the primary question): the same reversal on the 68 Öhman-only
#       proteins and the 18 Chae-listed proteins among the 471 (Spearman, T2D z permuted).
#
# HOW TO RUN
#   After steps 14 and 17s:   Rscript network/19_t2d_stories.R   (about 2 minutes; run_all.sh step 19)
#
# DATA AND PROVENANCE
#   Disease: $DISEASE_SCORES (default ~/Desktop/output/week_6/_shared/disease_scores.csv.gz; Venus week 6 build of
#   the Amar et al. 2024 inputs; fingerprinted by step 0). Exercise: MotrpacHumanPreSuspensionAnalysis v0.2.4 via
#   steps 1 / 1b (normalised responses) and 17s (per-cell adj. p). Network: step 14 (mnet, STRING v12 >= 700, Rhea).
#
# TECH STACK:  R 4.4; data.table, igraph, ggplot2, ggrepel, patchwork.
#
# INPUTS:  $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 01_nodes_EE.csv, 01_nodes_RE.csv,
#          01b_metab_nodes_EE.csv, 01b_metab_nodes_RE.csv, 17_node_cell_stats.csv; $DISEASE_SCORES
# OUTPUTS: $HACK_OUT/19_t2d_consensus_tests.csv  tissue: consensus reversal EE / RE / difference, sign-flip p
#          $HACK_OUT/19_t2d_node_tests.csv       sensitivity: single-study sets x node set x tissue (Spearman)
#          $HACK_OUT/19_t2d_subgraph_tests.csv   story x test: observed, null mean, p
#          $HACK_OUT/19_t2d_ego_networks.csv     story 1: per consensus protein: neighbours, edge weights, reversal
#          $HACK_OUT/19_t2d_story_nodes.csv      every node drawn in a story: T2D values, responses per tissue, cells
#          $HACK_OUT/19_t2d_consensus_cells.csv  story 3: consensus proteins x muscle cell: logFC, adj. p
#          $HACK_OUT/reports/19_t2d_stories.md   the numbers of the three stories
#          $HACK_FIG/19a_t2d_small_subgraphs.png, 19b_t2d_steiner_subgraph.png, 19c_t2d_consensus_proteins.png
#
# KNOWN LIMITS
#   One disease (T2D) because only its data are in a tissue we measure; only 8 proteins pass the two-study rule
#   (most of the 85 genome-wide consensus proteins are mitochondrial and not among the 471), so every test has
#   little power (the smallest possible sign-flip p with 8 proteins is 2/256 = 0.0078). Two consensus proteins
#   (DECR1, PEBP1) have no network edge and one (BLVRB) sits in a separate part of the network. The exercise data
#   are acute responses in healthy adults; a direction opposite to T2D is a match of directions, not evidence that
#   exercise treats or prevents T2D. Blood and adipose readouts compare a MUSCLE disease signature with other
#   tissues and are secondary.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel); library(patchwork) })

# Folders and inputs (override with environment variables).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz"))
dir.create(file.path(OUT, "reports"), recursive = TRUE, showWarnings = FALSE)
# Fixed seed, number of random draws, threshold, tissues (muscle first = primary) and colours.
SEED <- 20260926; N_PERM <- 10000L; ALPHA <- 0.05; TISSUES <- c("muscle", "blood", "adipose")
COL_LOW <- "#5E3C99"; COL_HIGH <- "#E66100"; COL_OTHER <- "#222222"   # lower in T2D / higher in T2D / other node
ARMLAB <- c(EE = "endurance vs control", RE = "resistance vs control")

# ---- 1. inputs --------------------------------------------------------------------------------------------
# T2D sets: one row per set and gene (smallest p kept), signed z from the p-value as in step 18.
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))[set %in% c("ohman_2021", "chae_2018") & !is.na(p), .(set, gene, logFC, p)]
ds <- ds[order(p)][!duplicated(paste(set, gene))]
ds[, z := sign(logFC) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)]
OH <- ds[set == "ohman_2021"]; CH <- ds[set == "chae_2018"]
# Consensus: significant in both studies (Öhman p < 0.05; Chae lists only significant proteins), same direction.
CONS <- merge(OH[, .(gene, oz = z, ologFC = logFC, op = p)], CH[, .(gene, cz = z, clogFC = logFC)], by = "gene")[op < ALPHA & sign(ologFC) == sign(clogFC)]
# The joint network (step 14): nodes and edges with their EE / RE dot-product weights.
E <- fread(file.path(OUT, "14_joint_edges.csv")); N <- fread(file.path(OUT, "14_joint_nodes.csv"))
gJ <- graph_from_data_frame(E[, .(node_a, node_b, edge_type, w_EE, w_RE, w_diff)], directed = FALSE, vertices = N[, .(node, node_type)])
PROT <- N[node_type == "protein", node]
# Normalised responses (steps 1 / 1b): genes and metabolites, one column per tissue x ome x time.
EEg <- fread(file.path(OUT, "01_nodes_EE.csv")); REg <- fread(file.path(OUT, "01_nodes_RE.csv"))
EEm <- fread(file.path(OUT, "01b_metab_nodes_EE.csv")); REm <- fread(file.path(OUT, "01b_metab_nodes_RE.csv"))
GENES471 <- EEg$gene_symbol
stopifnot(!anyDuplicated(GENES471), all(PROT %in% GENES471), all(N[node_type == "metabolite", node] %in% EEm$metabolite))
C471 <- CONS[gene %in% GENES471, gene]                           # consensus proteins among the 471 (story 3)
TERM <- intersect(C471, PROT)                                    # consensus proteins that are network nodes (stories 1, 2)
message(sprintf("consensus (both T2D studies, same direction): %d genome-wide, %d among the 471 (%s), %d network nodes",
                nrow(CONS), length(C471), paste(C471, collapse = ", "), length(TERM)))
# Mean normalised response of every node over the cells of one tissue (NA when the tissue has no cell for it).
tmean <- function(M, id, tis) { cols <- grep(paste0("^", tis, "_"), names(M), value = TRUE)
  v <- rowMeans(as.matrix(M[, ..cols]), na.rm = TRUE); v[is.nan(v)] <- NA; data.table(node = M[[id]], v = v) }
RESP <- rbindlist(lapply(TISSUES, function(t) {
  a <- rbind(tmean(EEg, "gene_symbol", t), tmean(EEm, "metabolite", t)); b <- rbind(tmean(REg, "gene_symbol", t), tmean(REm, "metabolite", t))
  data.table(tissue = t, node = a$node, rE = a$v, rR = b$v[match(a$node, b$node)]) }))
MUS <- RESP[tissue == "muscle"]; rE <- setNames(MUS$rE, MUS$node); rR <- setNames(MUS$rR, MUS$node)
Z <- setNames(OH$z, OH$gene)                                     # Öhman z (its sign = the consensus direction)
# Per-cell exercise statistics (step 17s): logFC and adj. p per node x tissue x ome x arm x time.
CELL <- fread(file.path(OUT, "17_node_cell_stats.csv"))
# PTMs for the figures (drawn, not tested): MoTrPAC muscle phosphosites (step 17s) and glycosylation sites (mnet / GlyGen).
PH <- fread(file.path(OUT, "17_phospho_site_stats.csv"))[tissue == "muscle" & arm %in% c("EE", "RE")]
GLY <- fread(file.path(OUT, "17_glygen_protein_annotation.csv"))[, .(protein, glyco_N_sites, glyco_O_sites, glyco_OGlcNAc_sites, crosstalk_residues)]

# ---- helpers ----------------------------------------------------------------------------------------------
# Permutation p-values with the +1 correction: two-sided around 0, one-sided "at most", two-sided around the
# null's own centre (for edge-weight means, whose null is not centred on 0).
# Per site: responds after EE / RE (any muscle time point, adj. p < 0.05), direction of its strongest response, label.
dir_of <- function(lfc, p, m) { i <- which(m & p < ALPHA); if (length(i)) sign(lfc[i][which.min(p[i])]) else 0 }
SITE <- PH[, .(ee = any(adj_p[arm == "EE"] < ALPHA, na.rm = TRUE), re = any(adj_p[arm == "RE"] < ALPHA, na.rm = TRUE),
               de = dir_of(logFC, adj_p, arm == "EE"), dr = dir_of(logFC, adj_p, arm == "RE"), minp = min(adj_p, na.rm = TRUE)), by = .(protein, site)]
SITE[, cls := fifelse(ee & re, "both", fifelse(ee, "EE", fifelse(re, "RE", "none")))]
ARW <- function(d) fifelse(d > 0, "\u2191", "\u2193")
SITE[, lab := fifelse(cls == "EE", paste0(site, ARW(de)), fifelse(cls == "RE", paste0(site, ARW(dr)),
               fifelse(cls == "both", fifelse(de == dr, paste0(site, ARW(de)), paste0(site, " E", ARW(de), " R", ARW(dr))), site)))]
p2 <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (1 + length(null)) }
p1lo <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(null <= obs + 1e-12)) / (1 + length(null)) }
p2c <- function(obs, null) { null <- null[is.finite(null)]; m <- mean(null); (1 + sum(abs(null - m) >= abs(obs - m) - 1e-12)) / (1 + length(null)) }
# Exact sign-flip test of mean(v) against 0 (all 2^n sign patterns when n <= 16, else N_PERM random ones).
signflip_p <- function(v) { v <- v[is.finite(v)]; n <- length(v); if (!n) return(NA_real_)
  S <- if (n <= 16) as.matrix(expand.grid(rep(list(c(-1, 1)), n))) else matrix(sample(c(-1, 1), n * N_PERM, TRUE), ncol = n)
  mean(abs(S %*% v) >= abs(sum(v)) - 1e-12) }
# Spearman reversal test (sensitivity sets): -Spearman(T2D z, response) per arm; the same shuffle of z for both arms.
reversal_test <- function(z, a, b) {
  rz <- rank(z); ra <- rank(a); rb <- rank(b); oA <- -cor(rz, ra); oB <- -cor(rz, rb)
  P <- replicate(N_PERM, sample(rz)); nA <- -drop(cor(ra, P)); nB <- -drop(cor(rb, P))
  data.table(n = length(z), rev_EE = oA, p_EE = p2(oA, nA), rev_RE = oB, p_RE = p2(oB, nB), diff = oA - oB, p_diff = p2(oA - oB, nA - nB)) }
# Degree-binned sampler: random sets with the same degree-bin make-up as `nodes`, drawn from `univ`.
deg_sampler <- function(g, univ, nodes) { d <- degree(g)[univ]; b <- setNames(cut(rank(d, ties.method = "first"), 5, labels = FALSE), univ)
  tb <- table(b[nodes]); function() unlist(lapply(names(tb), function(k) { pool <- univ[b == as.integer(k)]; pool[sample.int(length(pool), tb[[k]])] })) }
# Edges touching a node set (at least one end in it).
ego_edges <- function(v) E[node_a %in% v | node_b %in% v]
# Steiner tree (Kou-Markowsky-Berman approximation) connecting the terminals `term` in graph g.
steiner <- function(g, term) {
  D <- distances(g, term, term); kg <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE, diag = FALSE)
  el <- as_edgelist(mst(kg))
  nodes <- unique(unlist(lapply(seq_len(nrow(el)), function(i) shortest_paths(g, el[i, 1], el[i, 2], output = "vpath")$vpath[[1]]$name)))
  tr <- mst(induced_subgraph(g, nodes))                           # a spanning tree of the expanded paths
  repeat { lv <- V(tr)$name[degree(tr) == 1 & !V(tr)$name %in% term]; if (!length(lv)) break; tr <- delete_vertices(tr, lv) }
  tr }

# ---- 2. STORY 3 statistics first: consensus proteins, per tissue (exact sign flip) ---------------------------
set.seed(SEED)
CT <- rbindlist(lapply(TISSUES, function(t) { R <- RESP[tissue == t]; a <- setNames(R$rE, R$node)[C471]; b <- setNames(R$rR, R$node)[C471]
  s <- -sign(Z[C471]); ok <- is.finite(a) & is.finite(b)
  data.table(tissue = t, n = sum(ok), rev_EE = mean((s * a)[ok]), p_EE = signflip_p((s * a)[ok]), rev_RE = mean((s * b)[ok]), p_RE = signflip_p((s * b)[ok]),
             diff = mean((s * (a - b))[ok]), p_diff = signflip_p((s * (a - b))[ok]),
             n_rev_EE = sum((s * a)[ok] > 0), n_rev_RE = sum((s * b)[ok] > 0), n_EE_gt_RE = sum((s * (a - b))[ok] > 0)) }))
fwrite(CT, file.path(OUT, "19_t2d_consensus_tests.csv"))
print(CT)
# Sensitivity (single-study sets; not the primary question): Spearman reversal on the Öhman-only and Chae-listed proteins.
set.seed(SEED)
NT <- rbindlist(lapply(list(list("ohman_2021 only (p < 0.05)", OH[p < ALPHA & gene %in% GENES471]), list("chae_2018 only (listed)", CH[gene %in% GENES471])), function(s) {
  rbindlist(lapply(TISSUES, function(t) { y <- merge(s[[2]][, .(node = gene, z)], RESP[tissue == t], by = "node")[is.finite(rE) & is.finite(rR)]
    cbind(data.table(set = s[[1]], nodes = "all 471", tissue = t), reversal_test(y$z, y$rE, y$rR)) })) }))
fwrite(NT, file.path(OUT, "19_t2d_node_tests.csv"))

# ---- 3. STORY 1 — small subgraphs: each consensus protein and its direct neighbours --------------------------
set.seed(SEED)
UNIV <- intersect(OH$gene, PROT)                                 # network proteins Öhman measured (the draw pool)
samp <- deg_sampler(gJ, UNIV, TERM)
e1 <- ego_edges(TERM)
null1 <- t(replicate(N_PERM, { es <- ego_edges(samp()); c(wd = mean(es$w_diff), wE = mean(es$w_EE), wR = mean(es$w_RE)) }))
ST <- data.table(story = "1 small subgraphs", test = c("mean w_EE - w_RE, edges touching the consensus proteins", "mean w_EE (same edges)", "mean w_RE (same edges)"),
                 observed = c(mean(e1$w_diff), mean(e1$w_EE), mean(e1$w_RE)), null_mean = colMeans(null1),
                 p = c(p2c(mean(e1$w_diff), null1[, "wd"]), p2c(mean(e1$w_EE), null1[, "wE"]), p2c(mean(e1$w_RE), null1[, "wR"])),
                 side = "two-sided vs degree-matched random protein sets")
EGO <- rbindlist(lapply(TERM, function(v) { es <- ego_edges(v); nb <- setdiff(unique(c(es$node_a, es$node_b)), v)
  data.table(protein = v, t2d = fifelse(Z[v] < 0, "lower in T2D", "higher in T2D"), neighbours = paste(nb, collapse = ";"), n_edges = nrow(es),
             mean_w_EE = mean(es$w_EE), mean_w_RE = mean(es$w_RE), mean_w_diff = mean(es$w_diff),
             reversal_EE = -sign(Z[v]) * rE[v], reversal_RE = -sign(Z[v]) * rR[v]) }))
fwrite(EGO, file.path(OUT, "19_t2d_ego_networks.csv"))
S1 <- unique(c(TERM, unlist(strsplit(EGO$neighbours, ";"))))
print(EGO[, .(protein, t2d, neighbours, n_edges, w_EE = round(mean_w_EE, 3), w_RE = round(mean_w_RE, 3), rev_EE = round(reversal_EE, 3), rev_RE = round(reversal_RE, 3))])

# ---- 4. STORY 2 — larger subgraph: Steiner tree connecting the consensus proteins ----------------------------
set.seed(SEED)
memb <- components(gJ)$membership
big <- as.integer(names(which.max(table(memb[TERM]))))          # the connected part holding most consensus proteins
T2 <- TERM[memb[TERM] == big]; OUTSIDE <- setdiff(TERM, T2)
tree <- steiner(gJ, T2); S2 <- V(tree)$name
e2 <- as_data_frame(tree, what = "edges")
# Null: the same construction on degree-matched random terminal sets from the same connected part.
U2 <- UNIV[memb[UNIV] == big]; samp2 <- deg_sampler(gJ, U2, T2)
null2 <- t(replicate(N_PERM, { tr <- steiner(gJ, samp2()); w <- edge_attr(tr, "w_diff"); c(size = vcount(tr), wd = if (length(w)) mean(w) else NA) }))
ST <- rbind(ST, data.table(story = "2 Steiner tree",
  test = c("tree size (nodes) connecting the consensus proteins", "mean w_EE - w_RE over the tree edges", "mean w_EE over the tree edges", "mean w_RE over the tree edges"),
  observed = c(length(S2), mean(e2$w_diff), mean(e2$w_EE), mean(e2$w_RE)), null_mean = c(mean(null2[, "size"]), mean(null2[, "wd"], na.rm = TRUE), NA, NA),
  p = c(p1lo(length(S2), null2[, "size"]), p2c(mean(e2$w_diff), null2[, "wd"]), NA, NA),
  side = c("one-sided (smaller = closer than chance) vs trees on degree-matched random terminals", "two-sided vs the same random trees", "", "")))
message(sprintf("story 2: tree of %d nodes (%d consensus proteins; %s outside this part of the network): %s", length(S2), length(T2),
                paste(OUTSIDE, collapse = ", "), paste(S2, collapse = ", ")))

# ---- 5. STORY 3 tables: consensus proteins per muscle cell ----------------------------------------------------
ST <- rbind(ST, CT[tissue == "muscle", .(story = "3 consensus proteins", test = c("mean reversal endurance (muscle)", "mean reversal resistance (muscle)", "mean reversal endurance - resistance (muscle)"),
                                          observed = c(rev_EE, rev_RE, diff), null_mean = 0, p = c(p_EE, p_RE, p_diff), side = "exact sign-flip test (all 2^n patterns)")])
fwrite(ST, file.path(OUT, "19_t2d_subgraph_tests.csv"))
print(ST)
CC <- CELL[node %in% C471 & tissue == "muscle" & arm %in% c("EE", "RE")][, t2d := fifelse(Z[node] < 0, "lower in T2D", "higher in T2D")]
fwrite(CC, file.path(OUT, "19_t2d_consensus_cells.csv"))

# ---- 6. story node table ------------------------------------------------------------------------------------
cellsig <- CELL[tissue == "muscle" & arm %in% c("EE", "RE"), .(sig_cells = sum(adj_p < ALPHA, na.rm = TRUE),
  best = { i <- which.min(adj_p); if (length(i)) sprintf("%s %s %s logFC %.2f adj.p %.2g", ome[i], time[i], arm, logFC[i], adj_p[i]) else "" }), by = .(node, arm)]
wide <- dcast(RESP, node ~ tissue, value.var = c("rE", "rR"))
SN <- rbind(data.table(story = "1 small subgraphs", node = S1, role = fifelse(S1 %in% TERM, "consensus protein", "neighbour")),
            data.table(story = "2 Steiner tree", node = S2, role = fifelse(S2 %in% TERM, "consensus protein", "path node")),
            data.table(story = "3 consensus proteins", node = C471, role = "consensus protein"))
SN <- merge(SN, N[, .(node, node_type)], by = "node", all.x = TRUE)[is.na(node_type), node_type := "protein"]
SN <- merge(SN, OH[, .(node = gene, ohman_log2FC = logFC, ohman_p = p, ohman_z = z)], by = "node", all.x = TRUE)
SN <- merge(SN, CH[, .(node = gene, chae_log2FC = logFC)], by = "node", all.x = TRUE)
SN <- merge(SN, wide, by = "node", all.x = TRUE)
SN <- merge(SN, dcast(cellsig, node ~ arm, value.var = c("sig_cells", "best")), by = "node", all.x = TRUE)
setorder(SN, story, role, node); fwrite(SN, file.path(OUT, "19_t2d_story_nodes.csv"))

# ---- 7. figures -----------------------------------------------------------------------------------------------
# Both arms of a figure share one fill range and one edge-width range, so the two panels are directly comparable.
FILL <- function(lim) scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = "mean normalised\nmuscle response")
# One network panel: nodes filled by the muscle response of one arm, outlined by T2D direction (consensus proteins)
# or black (other nodes); edges drawn from `es` (node_a, node_b, w_EE, w_RE) with width = that arm's |weight|.
net_panel <- function(es, lay, arm, lim, wlim, title) {
  es <- copy(es)[, w := if (arm == "EE") w_EE else w_RE]
  es <- merge(merge(es, lay[, .(node_a = node, xa = x, ya = y)], by = "node_a"), lay[, .(node_b = node, xb = x, yb = y)], by = "node_b")
  nd <- copy(lay)[, r := if (arm == "EE") rE[node] else rR[node]]
  nd[, outline := fifelse(!node %in% C471, COL_OTHER, fifelse(Z[node] < 0, COL_LOW, COL_HIGH))]
  nd[, `:=`(shape = fifelse(node %in% PROT, 21L, 24L), sz = fifelse(node %in% C471, 6, 4.2), st = fifelse(node %in% C471, 2, 1))]
  ggplot() +
    geom_segment(data = es, aes(x = xa, y = ya, xend = xb, yend = yb, linewidth = abs(w), linetype = w < 0), colour = "grey45", alpha = 0.85) +
    tag_layers(tag_data(lay)) +
    geom_point(data = nd, aes(x, y, fill = r, colour = outline, shape = shape, size = sz, stroke = st)) +
    geom_text_repel(data = nd, aes(x, y, label = node, fontface = ifelse(node %in% C471, "bold", "plain")), size = 3, seed = SEED, max.time = 60, max.iter = 1e4,
                    box.padding = 0.35, min.segment.length = 0.2, segment.size = 0.2) +
    scale_colour_identity() + scale_shape_identity() + scale_size_identity() + FILL(lim) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    scale_linewidth(range = c(0.2, 2.6), limits = c(0, wlim), name = "|edge weight|") +
    coord_equal(clip = "off") + labs(title = title) + theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10), legend.position = "bottom")
}
# PTM tags on the proteins, drawn exactly as in the figure 17 pages (their default tags): stalks fan out clockwise from
# the upper right, one every 32 degrees; phosphosite pins = circle with a white "P" (red = responds after endurance
# only, blue = resistance only, purple = both; MoTrPAC muscle, adj. p < 0.05 at any time; one pin per responding site,
# strongest first, up to 6, then "+n"); glycosylation = SNFG squares (blue = N-linked GlcNAc, yellow = O-linked GalNAc,
# blue with a white dot = O-GlcNAc; number beyond = sites; database knowledge via mnet); gold star = crosstalk residue
# (a MoTrPAC phosphosite that is also an O-glycosylation site; number = residues). Same tags on both arm panels.
TAGCOL <- c(EE = "#E41A1C", RE = "#377EB8", both = "#984EA3")
MAX_PINS <- 6L
tag_data <- function(lay) {
  u <- max(diff(range(lay$x)), diff(range(lay$y)), 1e-6) * 0.0032      # one "pixel" of the page, in layout units
  rbindlist(lapply(lay$node, function(v) {
    st <- SITE[protein == v & cls != "none"][order(minp)]; g <- GLY[protein == v]
    it <- data.table(kind = rep("P", min(nrow(st), MAX_PINS)), col = TAGCOL[head(st$cls, MAX_PINS)], n = 1L)
    if (nrow(st) > MAX_PINS) it <- rbind(it, data.table(kind = "more", col = "#222222", n = nrow(st) - MAX_PINS))
    if (nrow(g)) { if (g$glyco_N_sites > 0) it <- rbind(it, data.table(kind = "sq", col = "#0072BC", n = g$glyco_N_sites))
                   if (g$glyco_O_sites > 0) it <- rbind(it, data.table(kind = "sq", col = "#FFD400", n = g$glyco_O_sites))
                   if (g$glyco_OGlcNAc_sites > 0) it <- rbind(it, data.table(kind = "sqo", col = "#0072BC", n = g$glyco_OGlcNAc_sites))
                   if (g$crosstalk_residues > 0) it <- rbind(it, data.table(kind = "star", col = "#FFD700", n = g$crosstalk_residues)) }
    if (!nrow(it)) return(NULL)
    xy <- lay[node == v]; r <- (if (v %in% C471) 12 else 9) * u
    a <- (70 - (seq_len(nrow(it)) - 1) * 32) * pi / 180                    # the page's -70 + 32k degrees, with y pointing up
    it[, `:=`(node = v, x0 = xy$x + r * cos(a), y0 = xy$y + r * sin(a), x1 = xy$x + (r + 13 * u) * cos(a), y1 = xy$y + (r + 13 * u) * sin(a),
              hx = xy$x + (r + 18.5 * u) * cos(a), hy = xy$y + (r + 18.5 * u) * sin(a), nx = xy$x + (r + 27 * u) * cos(a), ny = xy$y + (r + 27 * u) * sin(a))] })) }
tag_layers <- function(td) {
  if (is.null(td) || !nrow(td)) return(list())
  P <- td[kind == "P"]; Q <- td[kind %in% c("sq", "sqo")]; S <- td[kind == "star"]; M <- td[kind == "more"]; Nn <- td[n > 1 & kind != "more"]
  list(geom_segment(data = td, aes(x = x0, y = y0, xend = x1, yend = y1), colour = "#555555", linewidth = 0.4),
       geom_point(data = P, aes(hx, hy, colour = col), shape = 16, size = 3.4), geom_point(data = P, aes(hx, hy), shape = 1, size = 3.4, colour = "#222222", stroke = 0.35),
       geom_text(data = P, aes(hx, hy, label = "P"), colour = "white", fontface = "bold", size = 1.9),
       geom_point(data = Q, aes(hx, hy, colour = col), shape = 15, size = 3), geom_point(data = Q, aes(hx, hy), shape = 0, size = 3, colour = "#222222", stroke = 0.35),
       geom_point(data = Q[kind == "sqo"], aes(hx, hy), shape = 16, size = 1, colour = "white"),
       geom_text(data = S, aes(hx, hy, label = "\u2605"), colour = "#FFD700", size = 4.2),
       geom_text(data = M, aes(hx, hy, label = paste0("+", n)), colour = "#222222", fontface = "bold", size = 2.4),
       geom_text(data = Nn, aes(nx, ny, label = n), colour = "#222222", size = 2)) }
# Key for the PTM tags (as the figure 17 legend).
key_panel <- function() {
  k <- data.table(y = c(7.6, 6.8, 6, 3.9, 3.1, 2.3), lab = c("after endurance only", "after resistance only", "after both arms", "N-linked (GlcNAc)", "O-linked (GalNAc)", "O-GlcNAc"),
                  col = c(TAGCOL, "#0072BC", "#FFD400", "#0072BC"), kind = c("P", "P", "P", "sq", "sq", "sqo"))
  ggplot(k, aes(0, y)) +
    geom_point(data = k[kind == "P"], aes(colour = col), shape = 16, size = 4.2) + geom_point(data = k[kind == "P"], shape = 1, size = 4.2, colour = "#222222", stroke = 0.35) +
    geom_text(data = k[kind == "P"], label = "P", colour = "white", fontface = "bold", size = 2.3) +
    geom_point(data = k[kind != "P"], aes(colour = col), shape = 15, size = 3.8) + geom_point(data = k[kind != "P"], shape = 0, size = 3.8, colour = "#222222", stroke = 0.35) +
    geom_point(data = k[kind == "sqo"], shape = 16, size = 1.2, colour = "white") +
    annotate("text", 0, 1.2, label = "\u2605", colour = "#FFD700", size = 5) +
    geom_text(aes(0.3, y, label = lab), hjust = 0, size = 2.8) + annotate("text", 0.3, 1.2, label = "phospho = O-glyco residue", hjust = 0, size = 2.8) +
    annotate("text", -0.2, 8.2, label = "MoTrPAC phosphosites responding\n(muscle, adj. p < 0.05, any time;\none pin per site, up to 6, then +n)", hjust = 0, vjust = 0, size = 2.7, fontface = "bold") +
    annotate("text", -0.2, 4.5, label = "glycosylation, SNFG symbols\n(UniProt via mnet; number = sites)", hjust = 0, vjust = 0, size = 2.7, fontface = "bold") +
    scale_colour_identity() + coord_cartesian(xlim = c(-0.3, 2.6), ylim = c(0.6, 9.6), clip = "off") + theme_void() }
# Layout of a graph: Fruchterman-Reingold per connected piece, fixed seed.
layout_of <- function(g) { set.seed(SEED); L <- layout_components(g, layout = layout_with_kk); data.table(node = V(g)$name, x = L[, 1], y = L[, 2]) }
lim_of <- function(nodes) max(quantile(abs(c(rE[nodes], rR[nodes])), 0.95, na.rm = TRUE), 1e-6)
hist_panel <- function(v, obs, lab, title) ggplot(data.table(v = v[is.finite(v)]), aes(v)) + geom_histogram(bins = 40, fill = "grey75", colour = "grey55") +
  geom_vline(xintercept = obs, colour = "#D7301F", linewidth = 1) + labs(x = lab, y = "random draws", title = title) + theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9))
ARMCOL <- c(endurance = "#D7301F", resistance = "#2B8CBE")

# 19a — story 1: ego networks (edges touching a consensus protein only)
g1 <- graph_from_data_frame(e1[, .(node_a, node_b)], directed = FALSE); L1 <- layout_of(g1)
wl1 <- max(abs(c(e1$w_EE, e1$w_RE))); lim1 <- lim_of(S1)
eb <- melt(EGO[, .(protein = sprintf("%s (%s)", protein, t2d), endurance = mean_w_EE, resistance = mean_w_RE)], id.vars = "protein", variable.name = "arm", value.name = "w")
pe <- ggplot(eb, aes(w, protein, fill = arm)) + geom_col(position = position_dodge(0.75), width = 0.7) + geom_vline(xintercept = 0, colour = "grey50") +
  scale_fill_manual(values = ARMCOL, name = NULL) + labs(x = "mean weight of the edges touching the protein", y = NULL,
  title = sprintf("Edge weights around each protein; all edges: w_EE - w_RE %.3f vs degree-matched random sets, p %.2g", ST[1, observed], ST[1, p])) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
rb <- melt(EGO[, .(protein = sprintf("%s (%s)", protein, t2d), endurance = reversal_EE, resistance = reversal_RE)], id.vars = "protein", variable.name = "arm", value.name = "rev")
pr <- ggplot(rb, aes(rev, protein, colour = arm)) + geom_vline(xintercept = 0, colour = "grey50") + geom_point(size = 3.2, position = position_dodge(0.5)) +
  scale_colour_manual(values = ARMCOL, name = NULL) + labs(x = "reversal = -sign(T2D change) x muscle response (> 0: opposite to T2D)", y = NULL,
  title = "Each protein's own muscle response, oriented to its T2D direction") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
row19a <- (net_panel(e1, L1, "EE", lim1, wl1, "Endurance: node fill = muscle response, edge width = w_EE") |
         net_panel(e1, L1, "RE", lim1, wl1, "Resistance: node fill = muscle response, edge width = w_RE") | key_panel()) + plot_layout(widths = c(1, 1, 0.26), guides = "collect") & theme(legend.position = "bottom")
f19a <- row19a /
        (pe | pr) + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(title = sprintf("Figure 19a. The %d T2D proteins validated in both studies that are network nodes, each with its direct neighbours: endurance vs resistance", length(TERM)),
                  subtitle = "Consensus = Öhman 2021 p < 0.05 and listed by Chae 2018, same direction. Bold label + thick outline = consensus protein (purple = lower, orange = higher in T2D); dashed edge = negative weight.",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19a_t2d_small_subgraphs.png"), f19a, width = 20, height = 14, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19a_t2d_small_subgraphs.png"))

# 19b — story 2: Steiner tree
L2 <- layout_of(tree); wl2 <- max(abs(c(e2$w_EE, e2$w_RE))); lim2 <- lim_of(S2)
e2dt <- as.data.table(e2)[, .(node_a = from, node_b = to, w_EE, w_RE)]
st2 <- ST[story == "2 Steiner tree"]
row19b <- (net_panel(e2dt, L2, "EE", lim2, wl2, "Endurance: node fill = muscle response, edge width = w_EE") |
         net_panel(e2dt, L2, "RE", lim2, wl2, "Resistance: node fill = muscle response, edge width = w_RE") | key_panel()) + plot_layout(widths = c(1, 1, 0.26), guides = "collect") & theme(legend.position = "bottom")
f19b <- row19b /
        (hist_panel(null2[, "size"], length(S2), "tree size (nodes)", sprintf("Tree size vs trees on %s degree-matched random terminal sets: p %.2g (one-sided, smaller)", format(N_PERM, big.mark = ","), st2$p[1])) |
         hist_panel(null2[, "wd"], mean(e2$w_diff), "mean w_EE - w_RE over the tree edges", sprintf("Mean w_EE - w_RE (%.3f; w_EE %.3f, w_RE %.3f) vs the same random trees: p %.2g", st2$observed[2], st2$observed[3], st2$observed[4], st2$p[2]))) +
  plot_layout(heights = c(1.6, 1)) +
  plot_annotation(title = sprintf("Figure 19b. Smallest tree in the joint network connecting %d of the T2D proteins validated in both studies: endurance vs resistance", length(T2)),
                  subtitle = sprintf("Steiner tree (Kou-Markowsky-Berman) of %d nodes; bold = consensus protein (purple = lower, orange = higher in T2D), black outline = path node. %s is in a separate part of the network.",
                                     length(S2), paste(OUTSIDE, collapse = ", ")), theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19b_t2d_steiner_subgraph.png"), f19b, width = 20, height = 15, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19b_t2d_steiner_subgraph.png"))

# 19c — story 3: consensus set per tissue (primary) + single-study sensitivity, and the per-cell muscle response
fr <- rbind(CT[, .(lab = sprintf("consensus, %s (n = %d)", tissue, n), arm = "endurance", v = rev_EE, p = p_EE)], CT[, .(lab = sprintf("consensus, %s (n = %d)", tissue, n), arm = "resistance", v = rev_RE, p = p_RE)],
            CT[, .(lab = sprintf("consensus, %s (n = %d)", tissue, n), arm = "endurance - resistance", v = diff, p = p_diff)])
fr[, arm := factor(arm, levels = c("endurance", "resistance", "endurance - resistance"))][, lab := factor(lab, levels = rev(unique(lab)))]
pf <- ggplot(fr, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 3, position = position_dodge(0.6)) +
  geom_text(aes(label = sprintf("p %.2g", p)), position = position_dodge(0.6), vjust = -0.9, size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(ARMCOL, `endurance - resistance` = "#222222"), name = NULL) +
  labs(x = "mean reversal over the consensus proteins (> 0: opposite to T2D)", y = NULL, title = "Consensus proteins per tissue (exact sign-flip p; muscle = tissue of the T2D data)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
ns <- NT[tissue == "muscle"]
ns <- rbind(ns[, .(lab = sprintf("%s (n = %d)", set, n), arm = "endurance", v = rev_EE, p = p_EE)], ns[, .(lab = sprintf("%s (n = %d)", set, n), arm = "resistance", v = rev_RE, p = p_RE)],
            ns[, .(lab = sprintf("%s (n = %d)", set, n), arm = "endurance - resistance", v = diff, p = p_diff)])
ns[, arm := factor(arm, levels = c("endurance", "resistance", "endurance - resistance"))]
pn <- ggplot(ns, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 3, position = position_dodge(0.6)) +
  geom_text(aes(label = sprintf("p %.2g", p)), position = position_dodge(0.6), vjust = -0.9, size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(ARMCOL, `endurance - resistance` = "#222222"), name = NULL) +
  labs(x = "-Spearman(T2D z, muscle response)", y = NULL, title = "Sensitivity: each T2D study alone (muscle; T2D z permuted)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
CC[, cell := factor(paste(ome, time), levels = c("rna 0.5h", "rna 4h", "rna 24h", "prot 0.5h", "prot 4h", "prot 24h"))][, gl := sprintf("%s (%s)", node, t2d)][, armlab := ARMLAB[arm]]
pd <- ggplot(CC[!is.na(cell)], aes(cell, gl)) + geom_point(aes(fill = logFC, size = pmin(-log10(adj_p), 4)), shape = 21, colour = "grey30") +
  geom_point(data = CC[!is.na(cell) & adj_p < ALPHA], aes(cell, gl), shape = 21, size = 6.5, colour = "black", stroke = 1.1, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "logFC") + scale_size(range = c(1, 6), limits = c(0, 4), name = "-log10 adj. p\n(capped at 4)") +
  facet_wrap(~armlab) + labs(x = "muscle cell (ome, time)", y = NULL, title = sprintf("The %d consensus proteins among the 471: muscle response per cell (ring = adj. p < 0.05)", length(C471))) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), axis.text.x = element_text(angle = 30, hjust = 1))
# The consensus proteins' measured muscle phosphosites, per time point and arm (ring = adj. p < 0.05).
PC <- PH[protein %in% C471][, t2d := fifelse(Z[protein] < 0, "lower in T2D", "higher in T2D")][, gl := sprintf("%s %s (%s)", protein, site, t2d)][, armlab := ARMLAB[arm]]
PC[, time := factor(time, levels = c("0.5h", "4h", "24h"))]
pp <- ggplot(PC, aes(time, gl)) + geom_point(aes(fill = logFC, size = pmin(-log10(adj_p), 4)), shape = 21, colour = "grey30") +
  geom_point(data = PC[adj_p < ALPHA], aes(time, gl), shape = 21, size = 6.5, colour = "black", stroke = 1.1, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "logFC") + scale_size(range = c(1, 6), limits = c(0, 4), name = "-log10 adj. p\n(capped at 4)") +
  facet_wrap(~armlab) + labs(x = "muscle phosphoproteome, time after exercise", y = NULL,
  title = sprintf("Their measured muscle phosphosites (%d sites on %d proteins; %s have none): logFC per time (ring = adj. p < 0.05)", uniqueN(PC$gl), uniqueN(PC$protein),
                  paste(setdiff(C471, PC$protein), collapse = ", "))) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9))
f19c <- (pf | pn) / (pd | pp) + plot_layout(heights = c(0.8, 1.2)) +
  plot_annotation(title = "Figure 19c. The T2D proteins validated in both studies, one by one: endurance vs resistance",
                  subtitle = "Consensus = Öhman 2021 p < 0.05 and listed by Chae 2018 (FDR < 0.1), same direction. Reversal > 0 = exercise moves the protein opposite to T2D.",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19c_t2d_consensus_proteins.png"), f19c, width = 20, height = 12, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19c_t2d_consensus_proteins.png"))

# ---- 8. report ----------------------------------------------------------------------------------------------
f3 <- function(x) formatC(x, digits = 3, format = "fg"); fp <- function(x) formatC(x, digits = 2, format = "g")
row_md <- function(d) paste0("| ", apply(d, 1, paste, collapse = " | "), " |")
md <- c("# Step 19 — endurance vs resistance in the context of T2D (three stories)", "",
  "Disease data: Amar et al. 2024 (Cell Metab 36:1411), T2D muscle sets only. Only proteins validated in BOTH T2D datasets in the same direction",
  sprintf("are used (Öhman 2021 p < 0.05 and listed by Chae 2018): %d genome-wide, %d among the 471 (%s), %d network nodes.", nrow(CONS), length(C471), paste(C471, collapse = ", "), length(TERM)),
  "Reversal = -sign(T2D change) x mean normalised muscle response (> 0 = opposite to T2D). Random draws: 10,000; seed 20260926.", "",
  "## Consensus proteins per tissue (exact sign-flip test)", "", "| tissue | n | reversal EE | p | reversal RE | p | EE - RE | p | proteins reversed EE / RE / EE > RE |", "|---|---|---|---|---|---|---|---|---|",
  row_md(CT[, .(tissue, n, f3(rev_EE), fp(p_EE), f3(rev_RE), fp(p_RE), f3(diff), fp(p_diff), sprintf("%d / %d / %d", n_rev_EE, n_rev_RE, n_EE_gt_RE))]), "",
  "## Story 1 — each consensus protein with its neighbours (figure 19a)", "", "| protein | T2D | neighbours | edges | mean w_EE | mean w_RE | reversal EE | reversal RE |", "|---|---|---|---|---|---|---|---|",
  row_md(EGO[, .(protein, t2d, neighbours, n_edges, f3(mean_w_EE), f3(mean_w_RE), f3(reversal_EE), f3(reversal_RE))]), "",
  sprintf("## Story 2 — Steiner tree (figure 19b): %d nodes: %s (%s in a separate part of the network)", length(S2), paste(S2, collapse = ", "), paste(OUTSIDE, collapse = ", ")), "",
  "## Subgraph and set tests", "", "| story | test | observed | null mean | p | test type |", "|---|---|---|---|---|---|",
  row_md(ST[, .(story, test, f3(observed), ifelse(is.na(null_mean), "", f3(null_mean)), ifelse(is.na(p), "", fp(p)), side)]), "",
  "## Sensitivity: single-study sets (not the primary question)", "", "| set | tissue | n | reversal EE | p | reversal RE | p | EE - RE | p |", "|---|---|---|---|---|---|---|---|---|",
  row_md(NT[, .(set, tissue, n, f3(rev_EE), fp(p_EE), f3(rev_RE), fp(p_RE), f3(diff), fp(p_diff))]), "",
  "## Figures (caption skeletons)", "",
  "- **19a** — consensus proteins that are network nodes, each with its direct neighbours (EE vs RE panels); edge weights around each protein; each protein's reversal.",
  "- **19b** — Steiner tree connecting the consensus proteins (EE vs RE panels); tree size and edge-difference null distributions.",
  "- **19c** — consensus set per tissue (sign flip), single-study sensitivity, consensus proteins per muscle cell.")
writeLines(md, file.path(OUT, "reports", "19_t2d_stories.md")); message("-> ", file.path(OUT, "reports", "19_t2d_stories.md"))
