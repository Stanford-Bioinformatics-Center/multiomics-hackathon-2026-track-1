#!/usr/bin/env Rscript
# =====================================================================================================
# 19_t2d_stories.R — STEP 19 (FIGURES 19a / 19b / 19c): ENDURANCE VS RESISTANCE IN THE CONTEXT OF T2D
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Our key question is how endurance (EE) and resistance (RE) exercise differ in the context of preventing
#   disease. The disease data come from Amar et al. 2024 (Cell Metab 36:1411; doi 10.1016/j.cmet.2023.12.021).
#   Of its 9 disease sets only the two type 2 diabetes (T2D) sets are from a tissue we measure (skeletal muscle;
#   the others are heart or liver), so this step tells three T2D stories, each comparing EE with RE on the
#   network, and each with a permutation test:
#     STORY 1 — small subgraph: the T2D-altered proteins' own connected pieces of the network.
#     STORY 2 — larger subgraph: the T2D proteins plus the "connector" nodes (proteins or metabolites) that link
#               two or more of them.
#     STORY 3 — replication: does an independent T2D cohort (Chae 2018) reproduce the EE-vs-RE result, and what do
#               the proteins validated in BOTH studies (the consensus set of step 17) do?
#
# WHAT THIS SCRIPT DOES (plain language)
#   T2D direction per protein: Öhman 2021 muscle proteome (full table; log2 T2D / normal glucose tolerance, p) and
#   Chae 2018 (published significant proteins only), as processed in the Venus project week 6; signed z from the
#   p-value (sign from the log2 change), exactly as step 18. "T2D-altered" = Öhman p < 0.05 (Chae: listed).
#   Exercise response per node: the mean of its normalised log fold changes (the same numbers the embeddings and
#   edge weights are built from) over the cells of one tissue. MUSCLE is primary (the T2D data are muscle);
#   blood and adipose are reported as secondary readouts (a systemic view of the same proteins, not tissue-matched).
#   REVERSAL (per arm) = minus the Spearman correlation between T2D z and the exercise response across the
#   T2D-altered nodes: positive = exercise moves proteins that are higher in T2D down and those lower in T2D up.
#   Rank correlation is unaffected by an overall shift of all responses (the "most things go up" imbalance that
#   misled the sign counts in the step 17 pages). Tests:
#     - each arm vs 0 and the EE - RE difference: permutation of the T2D z across the same nodes (10,000; the
#       same permutation for both arms, so the difference has a proper null); sensitivity: per-node arm swap.
#     - subgraph connectivity (stories 1, 2): observed largest connected piece vs 10,000 random node sets of the
#       same size drawn from the Öhman-measured network proteins within 5 degree bins (Menche et al. 2015).
#     - EE vs RE edge weights inside a subgraph: mean w_EE - w_RE vs 10,000 random CONNECTED subgraphs of the
#       same size grown in the same network (two-sided).
#     - small pieces and the consensus set: exact sign-flip test on per-node (reversal EE - reversal RE).
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
# OUTPUTS: $HACK_OUT/19_t2d_node_tests.csv      set x node universe x tissue: reversal EE / RE, difference, p-values
#          $HACK_OUT/19_t2d_subgraph_tests.csv  story x test: observed, null mean, p
#          $HACK_OUT/19_t2d_components.csv      story 1: every connected piece: members, edge weights, tests
#          $HACK_OUT/19_t2d_story_nodes.csv     every node drawn in a story: T2D values, responses per tissue, cells
#          $HACK_OUT/19_t2d_consensus_cells.csv story 3: consensus proteins x muscle cell: logFC, adj. p
#          $HACK_OUT/reports/19_t2d_stories.md  the numbers of the three stories
#          $HACK_FIG/19a_t2d_small_subgraph.png, 19b_t2d_connector_subgraph.png, 19c_t2d_replication.png
#
# KNOWN LIMITS
#   One disease (T2D) because only its data are in a tissue we measure. Öhman's p is a 4-group ANOVA p; Chae lists
#   only significant proteins (FDR < 0.1), so its test uses fewer, pre-selected proteins. The exercise data are
#   acute responses in healthy adults; a direction opposite to T2D is a match of directions, not evidence that
#   exercise treats or prevents T2D. Blood and adipose readouts compare a MUSCLE disease signature with other
#   tissues and are secondary. Story 2's connectors carry no T2D value of their own.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel); library(patchwork) })

# Folders and inputs (override with environment variables).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz"))
dir.create(file.path(OUT, "reports"), recursive = TRUE, showWarnings = FALSE)
# Fixed seed, permutation count, thresholds, tissues (muscle first = primary) and colours.
SEED <- 20260926; N_PERM <- 10000L; ALPHA <- 0.05; TISSUES <- c("muscle", "blood", "adipose")
COL_LOW <- "#5E3C99"; COL_HIGH <- "#E66100"; COL_CONN <- "#222222"   # lower in T2D / higher in T2D / connector
ARMLAB <- c(EE = "endurance vs control", RE = "resistance vs control")

# ---- 1. inputs --------------------------------------------------------------------------------------------
# T2D sets: one row per set and gene (smallest p kept), signed z from the p-value as in step 18.
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))[set %in% c("ohman_2021", "chae_2018") & !is.na(p), .(set, gene, logFC, p)]
ds <- ds[order(p)][!duplicated(paste(set, gene))]
ds[, z := sign(logFC) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)]
OH <- ds[set == "ohman_2021"]; CH <- ds[set == "chae_2018"]
# The joint network (step 14): nodes and edges with their EE / RE dot-product weights.
E <- fread(file.path(OUT, "14_joint_edges.csv")); N <- fread(file.path(OUT, "14_joint_nodes.csv"))
gJ <- graph_from_data_frame(E[, .(node_a, node_b, edge_type, w_EE, w_RE, w_diff)], directed = FALSE, vertices = N[, .(node, node_type)])
PROT <- N[node_type == "protein", node]
# Normalised responses (steps 1 / 1b): genes and metabolites, one column per tissue x ome x time.
EEg <- fread(file.path(OUT, "01_nodes_EE.csv")); REg <- fread(file.path(OUT, "01_nodes_RE.csv"))
EEm <- fread(file.path(OUT, "01b_metab_nodes_EE.csv")); REm <- fread(file.path(OUT, "01b_metab_nodes_RE.csv"))
GENES471 <- EEg$gene_symbol
stopifnot(!anyDuplicated(GENES471), all(PROT %in% GENES471), all(N[node_type == "metabolite", node] %in% EEm$metabolite))
# Mean normalised response of every node over the cells of one tissue (NA when the tissue has no cell for it).
tmean <- function(M, id, tis) { cols <- grep(paste0("^", tis, "_"), names(M), value = TRUE)
  v <- rowMeans(as.matrix(M[, ..cols]), na.rm = TRUE); v[is.nan(v)] <- NA; data.table(node = M[[id]], v = v) }
RESP <- rbindlist(lapply(TISSUES, function(t) {
  a <- rbind(tmean(EEg, "gene_symbol", t), tmean(EEm, "metabolite", t)); b <- rbind(tmean(REg, "gene_symbol", t), tmean(REm, "metabolite", t))
  data.table(tissue = t, node = a$node, rE = a$v, rR = b$v[match(a$node, b$node)]) }))
# Per-cell exercise statistics (step 17s): logFC and adj. p per node x tissue x ome x arm x time.
CELL <- fread(file.path(OUT, "17_node_cell_stats.csv"))

# ---- helpers ----------------------------------------------------------------------------------------------
# Permutation p-values with the +1 correction (never exactly 0): two-sided around 0, one-sided "at least".
p2 <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (1 + length(null)) }
p1 <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(null >= obs - 1e-12)) / (1 + length(null)) }
# Two-sided p around the null's own centre (for edge-weight means, whose null is not centred on 0).
p2c <- function(obs, null) { null <- null[is.finite(null)]; m <- mean(null); (1 + sum(abs(null - m) >= abs(obs - m) - 1e-12)) / (1 + length(null)) }
# Reversal test: -Spearman(T2D z, response) per arm; null = the same shuffle of z for both arms; sensitivity =
# per-node arm swap (exchanging a node's EE and RE responses at random).
reversal_test <- function(z, a, b) {
  rz <- rank(z); ra <- rank(a); rb <- rank(b); n <- length(z)
  oA <- -cor(rz, ra); oB <- -cor(rz, rb)
  P <- replicate(N_PERM, sample(rz))                                      # n x N_PERM shuffled T2D ranks
  nA <- -drop(cor(ra, P)); nB <- -drop(cor(rb, P))
  S <- matrix(runif(n * N_PERM) < 0.5, n)                                 # which nodes swap arms
  sA <- -drop(cor(rz, apply(ifelse(S, b, a), 2, rank))); sB <- -drop(cor(rz, apply(ifelse(S, a, b), 2, rank)))
  data.table(n = n, rev_EE = oA, p_EE = p2(oA, nA), rev_RE = oB, p_RE = p2(oB, nB), diff = oA - oB,
             p_diff = p2(oA - oB, nA - nB), p_diff_armswap = p2(oA - oB, sA - sB)) }
# Exact sign-flip test of mean(v) against 0 (all 2^n sign patterns when n <= 16, else N_PERM random ones).
signflip_p <- function(v) { v <- v[is.finite(v)]; n <- length(v); if (!n) return(NA_real_)
  S <- if (n <= 16) as.matrix(expand.grid(rep(list(c(-1, 1)), n))) else matrix(sample(c(-1, 1), n * N_PERM, TRUE), ncol = n)
  mean(abs(S %*% v) >= abs(sum(v)) - 1e-12) }
# Random connected subgraph of k vertices: start at a random vertex, repeatedly add a random neighbour of the set.
nbr_list <- function(g) lapply(as_adj_list(g, mode = "all"), as.integer)
grow <- function(nb, k) { repeat { s <- sample.int(length(nb), 1); set <- s; fr <- setdiff(nb[[s]], set)
  while (length(set) < k && length(fr)) { x <- fr[sample.int(length(fr), 1)]; set <- c(set, x); fr <- setdiff(union(fr, nb[[x]]), set) }
  if (length(set) == k) return(set) } }
# Mean of an edge attribute inside the subgraph induced by vertex set v.
mean_edge <- function(g, v, attr = "w_diff") { w <- edge_attr(induced_subgraph(g, v), attr); if (length(w)) mean(w) else NA_real_ }
# Degree-binned sampler: random sets with the same degree-bin make-up as `nodes`, drawn from `univ`.
deg_sampler <- function(g, univ, nodes) { d <- degree(g)[univ]; b <- setNames(cut(rank(d, ties.method = "first"), 5, labels = FALSE), univ)
  tb <- table(b[nodes]); function() unlist(lapply(names(tb), function(k) { pool <- univ[b == as.integer(k)]; pool[sample.int(length(pool), tb[[k]])] })) }
lcc_size <- function(g, v) { if (!length(v)) return(0L); max(components(induced_subgraph(g, v))$csize) }

# ---- 2. node-level tests: reversal per arm and EE - RE, per tissue ------------------------------------------
# Node universes: the network's proteins (stories 1, 2) and all 471 proteins (story 3, more power).
set.seed(SEED)
NT <- rbindlist(lapply(list(list("ohman_2021", "network", OH[p < ALPHA & gene %in% PROT]), list("ohman_2021", "all 471", OH[p < ALPHA & gene %in% GENES471]),
                            list("chae_2018", "network", CH[gene %in% PROT]), list("chae_2018", "all 471", CH[gene %in% GENES471])), function(s) {
  rbindlist(lapply(TISSUES, function(t) { y <- merge(s[[3]][, .(node = gene, z)], RESP[tissue == t], by = "node")[is.finite(rE) & is.finite(rR)]
    if (nrow(y) < 5) return(NULL); cbind(data.table(set = s[[1]], nodes = s[[2]], tissue = t), reversal_test(y$z, y$rE, y$rR)) })) }))
# Benjamini-Hochberg over the three tissues within each set x universe (muscle is the pre-specified primary).
NT[, q_diff := p.adjust(p_diff, "BH"), by = .(set, nodes)]
fwrite(NT, file.path(OUT, "19_t2d_node_tests.csv"))
print(NT[, .(set, nodes, tissue, n, rev_EE = round(rev_EE, 3), p_EE, rev_RE = round(rev_RE, 3), p_RE, diff = round(diff, 3), p_diff, p_diff_armswap)])

# ---- 3. STORY 1 — small subgraph: the T2D proteins' own connected pieces ------------------------------------
set.seed(SEED)
UNIV <- intersect(OH$gene, PROT)                                   # network proteins measured by Öhman
DN <- OH[p < ALPHA & gene %in% PROT, gene]                         # T2D-altered network proteins
Z <- setNames(OH$z, OH$gene)
MUS <- RESP[tissue == "muscle"]; rE <- setNames(MUS$rE, MUS$node); rR <- setNames(MUS$rR, MUS$node)
sg1 <- induced_subgraph(gJ, DN); cm1 <- components(sg1)
# Whole-subgraph connectivity vs degree-matched random sets (largest piece, number of edges).
samp <- deg_sampler(gJ, UNIV, DN)
null1 <- t(replicate(N_PERM, { v <- samp(); c(lcc = lcc_size(gJ, v), ne = ecount(induced_subgraph(gJ, v))) }))
ST <- data.table(story = "1 small subgraph", test = c("largest connected piece (T2D proteins)", "edges among T2D proteins"),
                 observed = c(max(cm1$csize), ecount(sg1)), null_mean = colMeans(null1),
                 p = c(p1(max(cm1$csize), null1[, "lcc"]), p1(ecount(sg1), null1[, "ne"])), side = "one-sided (more than chance)")
# Every piece with >= 2 proteins: edge weights, size-matched connected null for w_EE - w_RE, sign-flip on nodes.
nbJ <- nbr_list(gJ); gP <- induced_subgraph(gJ, PROT); nbP <- nbr_list(gP)
null_wd <- list()   # cache: size -> null distribution of mean w_diff in random connected protein subgraphs
wd_null <- function(k) { key <- as.character(k); if (is.null(null_wd[[key]])) null_wd[[key]] <<- replicate(N_PERM, mean_edge(gP, grow(nbP, k))); null_wd[[key]] }
COMP <- rbindlist(lapply(which(cm1$csize >= 2), function(k) { v <- names(cm1$membership)[cm1$membership == k]; es <- E[node_a %in% v & node_b %in% v]
  vd <- -sign(Z[v]) * (rE[v] - rR[v])                                   # per node: reversal EE minus reversal RE (muscle)
  data.table(component = NA_integer_, size = length(v), members = paste(v, collapse = ";"), n_edges = nrow(es),
             mean_w_EE = mean(es$w_EE), mean_w_RE = mean(es$w_RE), mean_w_diff = mean(es$w_diff), p_w_diff = p2c(mean(es$w_diff), wd_null(length(v))),
             rev_EE = mean(-sign(Z[v]) * rE[v]), rev_RE = mean(-sign(Z[v]) * rR[v]), p_signflip = signflip_p(vd)) }))
COMP <- COMP[order(-size, mean_w_diff)][, component := .I]
fwrite(COMP, file.path(OUT, "19_t2d_components.csv"))
S1 <- unlist(strsplit(COMP$members, ";"))                            # story 1 nodes: all pieces with >= 2 proteins
# Pooled over all pieces: mean edge difference vs degree-matched random sets' induced edges.
wd1 <- mean_edge(gJ, S1); nd1 <- replicate(N_PERM, mean_edge(gJ, samp()))
ST <- rbind(ST, data.table(story = "1 small subgraph", test = "mean w_EE - w_RE over the pieces' edges", observed = wd1, null_mean = mean(nd1, na.rm = TRUE),
                           p = p2c(wd1, nd1), side = "two-sided vs degree-matched random sets"))
print(COMP[, .(component, size, members, n_edges, mean_w_EE = round(mean_w_EE, 3), mean_w_RE = round(mean_w_RE, 3), p_w_diff, rev_EE = round(rev_EE, 3), rev_RE = round(rev_RE, 3), p_signflip)])

# ---- 4. STORY 2 — larger subgraph: T2D proteins + connectors (linked to >= 2 of them) -------------------------
set.seed(SEED)
A <- as_adjacency_matrix(gJ, sparse = TRUE); A@x[] <- 1; VN <- rownames(A)
# Connectors of a seed set: nodes outside it adjacent to at least two seeds; the subgraph is seeds + connectors.
with_connectors <- function(seeds) { ind <- as.numeric(VN %in% seeds); cnt <- as.numeric(A %*% ind); c(seeds, VN[cnt >= 2 & ind == 0]) }
lcc_nodes <- function(v) { sg <- induced_subgraph(gJ, v); cm <- components(sg); names(cm$membership)[cm$membership == which.max(cm$csize)] }
V2all <- with_connectors(DN); S2 <- lcc_nodes(V2all)
CONN <- setdiff(S2, DN)
# Connectivity vs degree-matched random seed sets passed through the same connector rule.
null2 <- replicate(N_PERM, length(lcc_nodes(with_connectors(samp()))))
# EE vs RE edge weights inside the subgraph vs random connected subgraphs of the same size (whole joint network).
wd2 <- mean_edge(gJ, S2); nd2 <- replicate(N_PERM, mean_edge(gJ, grow(nbJ, length(S2))))
# Node-level reversal on the T2D proteins inside the subgraph (muscle).
y2 <- data.table(node = intersect(S2, DN))[, `:=`(z = Z[node], a = rE[node], b = rR[node])][is.finite(a) & is.finite(b)]
rt2 <- reversal_test(y2$z, y2$a, y2$b)
ST <- rbind(ST, data.table(story = "2 connector subgraph",
  test = c("largest connected piece (T2D proteins + connectors)", "mean w_EE - w_RE over its edges", "mean w_EE over its edges", "mean w_RE over its edges",
           "muscle reversal, endurance (T2D proteins inside)", "muscle reversal, resistance (T2D proteins inside)", "muscle reversal, endurance - resistance"),
  observed = c(length(S2), wd2, mean_edge(gJ, S2, "w_EE"), mean_edge(gJ, S2, "w_RE"), rt2$rev_EE, rt2$rev_RE, rt2$diff),
  null_mean = c(mean(null2), mean(nd2, na.rm = TRUE), NA, NA, 0, 0, 0),
  p = c(p1(length(S2), null2), p2c(wd2, nd2), NA, NA, rt2$p_EE, rt2$p_RE, rt2$p_diff),
  side = c("one-sided vs degree-matched seeds + same rule", "two-sided vs random connected subgraphs of the same size", "", "",
           "two-sided, T2D z permuted", "two-sided, T2D z permuted", "two-sided, T2D z permuted")))
# Connectors' edges to T2D proteins: which arm co-regulates them more strongly.
ce <- E[(node_a %in% CONN & node_b %in% DN) | (node_b %in% CONN & node_a %in% DN)]
ce <- ce[(node_a %in% S2) & (node_b %in% S2)]
message(sprintf("story 2: %d nodes (%d T2D proteins, %d connectors: %s); %d edges", length(S2), length(intersect(S2, DN)), length(CONN),
                paste(CONN, collapse = ", "), ecount(induced_subgraph(gJ, S2))))

# ---- 5. STORY 3 — replication (Chae 2018) and the proteins validated in both studies ---------------------------
CONS <- merge(OH[, .(gene, oz = z, ologFC = logFC, op = p)], CH[, .(gene, cz = z, clogFC = logFC)], by = "gene")[op < ALPHA & sign(ologFC) == sign(clogFC)]
CONS471 <- CONS[gene %in% GENES471]
vc <- -sign(CONS471$oz) * (rE[CONS471$gene] - rR[CONS471$gene])
ST <- rbind(ST, data.table(story = "3 replication",
  test = c("consensus proteins: mean reversal EE - RE (muscle)", "consensus proteins: mean reversal EE (muscle)", "consensus proteins: mean reversal RE (muscle)"),
  observed = c(mean(vc), mean(-sign(CONS471$oz) * rE[CONS471$gene]), mean(-sign(CONS471$oz) * rR[CONS471$gene])), null_mean = 0,
  p = c(signflip_p(vc), signflip_p(-sign(CONS471$oz) * rE[CONS471$gene]), signflip_p(-sign(CONS471$oz) * rR[CONS471$gene])),
  side = "exact sign-flip test (all 2^n patterns)"))
CC <- CELL[node %in% CONS471$gene & tissue == "muscle" & arm %in% c("EE", "RE")][, t2d := fifelse(Z[node] < 0, "lower in T2D", "higher in T2D")]
fwrite(CC, file.path(OUT, "19_t2d_consensus_cells.csv"))
fwrite(ST, file.path(OUT, "19_t2d_subgraph_tests.csv"))
print(ST)

# ---- 6. story node table ------------------------------------------------------------------------------------
cellsig <- CELL[tissue == "muscle" & arm %in% c("EE", "RE"), .(sig_cells = sum(adj_p < ALPHA, na.rm = TRUE),
  best = { i <- which.min(adj_p); if (length(i)) sprintf("%s %s %s logFC %.2f adj.p %.2g", ome[i], time[i], arm[i], logFC[i], adj_p[i]) else "" }), by = .(node, arm)]
wide <- dcast(RESP, node ~ tissue, value.var = c("rE", "rR"))
SN <- rbind(data.table(story = "1 small subgraph", node = S1, role = "T2D protein"),
            data.table(story = "2 connector subgraph", node = S2, role = fifelse(S2 %in% DN, "T2D protein", "connector")),
            data.table(story = "3 replication", node = CONS471$gene, role = "consensus (both studies)"))
SN <- merge(SN, N[, .(node, node_type)], by = "node", all.x = TRUE)[is.na(node_type), node_type := "protein"]
SN <- merge(SN, OH[, .(node = gene, ohman_log2FC = logFC, ohman_p = p, ohman_z = z)], by = "node", all.x = TRUE)
SN <- merge(SN, CH[, .(node = gene, chae_log2FC = logFC)], by = "node", all.x = TRUE)
SN <- merge(SN, wide, by = "node", all.x = TRUE)
SN <- merge(SN, dcast(cellsig, node ~ arm, value.var = c("sig_cells", "best")), by = "node", all.x = TRUE)
setorder(SN, story, role, node); fwrite(SN, file.path(OUT, "19_t2d_story_nodes.csv"))

# ---- 7. figures -----------------------------------------------------------------------------------------------
# One network panel: nodes filled by the muscle response of one arm, outlined by T2D direction; edges by that arm's weight.
# (both arms of a figure share one fill range and one edge-width range, so the two panels are directly comparable)
FILL <- function(lim) scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = "mean normalised\nmuscle response")
net_panel <- function(nodes, lay, arm, lim, wlim, title, label_all = TRUE) {
  es <- E[node_a %in% nodes & node_b %in% nodes]; es[, w := if (arm == "EE") w_EE else w_RE]
  es <- merge(merge(es, lay[, .(node_a = node, xa = x, ya = y)], by = "node_a"), lay[, .(node_b = node, xb = x, yb = y)], by = "node_b")
  nd <- copy(lay)[, r := if (arm == "EE") rE[node] else rR[node]]
  nd[, outline := fifelse(!node %in% DN, COL_CONN, fifelse(Z[node] < 0, COL_LOW, COL_HIGH))]
  nd[, shape := fifelse(node %in% PROT, 21L, 24L)]
  ggplot() +
    geom_segment(data = es, aes(x = xa, y = ya, xend = xb, yend = yb, linewidth = abs(w), linetype = w < 0), colour = "grey45", alpha = 0.8) +
    geom_point(data = nd, aes(x, y, fill = r, colour = outline, shape = shape), size = if (label_all) 5 else 3.6, stroke = 1.3) +
    geom_text_repel(data = nd, aes(x, y, label = node), size = if (label_all) 3 else 2.3, seed = SEED, max.time = 60, max.iter = 1e4,
                    box.padding = 0.3, min.segment.length = 0.2, segment.size = 0.2) +
    scale_colour_identity() + scale_shape_identity() + FILL(lim) + scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    scale_linewidth(range = c(0.2, 2.6), limits = c(0, wlim), name = "|edge weight|") +
    coord_equal(clip = "off") + labs(title = title) + theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10), legend.position = "right")
}
# Layout of a node set: Fruchterman-Reingold per connected piece, fixed seed.
layout_of <- function(nodes) { set.seed(SEED); sg <- induced_subgraph(gJ, nodes); L <- layout_components(sg, layout = layout_with_fr)
  data.table(node = V(sg)$name, x = L[, 1], y = L[, 2]) }
# Scatter of T2D z vs muscle response per arm, with the reversal statistic and its permutation p.
rev_scatter <- function(tab, zv, title, labs_nodes = character()) {
  d <- rbind(data.table(node = names(zv), z = zv, r = rE[names(zv)], arm = ARMLAB[["EE"]]), data.table(node = names(zv), z = zv, r = rR[names(zv)], arm = ARMLAB[["RE"]]))[is.finite(r)]
  ann <- data.table(arm = ARMLAB[c("EE", "RE")], lab = c(sprintf("reversal %.2f, p %.3g", tab$rev_EE, tab$p_EE), sprintf("reversal %.2f, p %.3g", tab$rev_RE, tab$p_RE)))
  ggplot(d, aes(z, r)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(aes(colour = z < 0), size = 1.8, alpha = 0.85) + scale_colour_manual(values = c(`TRUE` = COL_LOW, `FALSE` = COL_HIGH), guide = "none") +
    geom_text_repel(data = d[node %in% labs_nodes], aes(label = node), size = 2.4, seed = SEED, max.time = 60, max.iter = 1e4, min.segment.length = 0.1) +
    geom_label(data = ann, aes(x = -Inf, y = -Inf, label = lab), hjust = -0.05, vjust = -0.3, size = 3, label.size = 0, fill = "white", alpha = 0.85, inherit.aes = FALSE) +
    facet_wrap(~arm) + labs(x = "T2D z (signed; < 0 = lower in T2D)", y = "mean normalised muscle response", title = title) +
    theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10))
}
wlim_of <- function(nodes) { es <- E[node_a %in% nodes & node_b %in% nodes]; max(abs(c(es$w_EE, es$w_RE)), 1e-6) }
lim_of <- function(nodes) max(quantile(abs(c(rE[nodes], rR[nodes])), 0.95, na.rm = TRUE), 1e-6)
tab_m <- NT[set == "ohman_2021" & nodes == "network" & tissue == "muscle"]

# 19a — story 1
L1 <- layout_of(S1); lim1 <- lim_of(S1); wl1 <- wlim_of(S1)
cbar <- melt(COMP[, .(piece = sprintf("%d: %s", component, gsub(";", " · ", members)), `endurance (w_EE)` = mean_w_EE, `resistance (w_RE)` = mean_w_RE, p_w_diff)],
             id.vars = c("piece", "p_w_diff"), variable.name = "arm", value.name = "w")
cbar[, piece := factor(piece, levels = rev(unique(piece)))]
pc <- ggplot(cbar, aes(w, piece, fill = arm)) + geom_col(position = position_dodge(0.75), width = 0.7) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_text(data = unique(cbar[, .(piece, p_w_diff)]), aes(x = Inf, y = piece, label = sprintf("p %.2g", p_w_diff)), inherit.aes = FALSE, hjust = 1.05, size = 2.6) +
  scale_fill_manual(values = c(`endurance (w_EE)` = "#D7301F", `resistance (w_RE)` = "#2B8CBE"), name = NULL) +
  labs(x = "mean edge weight in the piece (dot product)", y = NULL, title = "Edge weights per piece (p: w_EE - w_RE vs random connected subgraphs of the same size)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
f19a <- (net_panel(S1, L1, "EE", lim1, wl1, "Endurance: node fill = muscle response, edge width = w_EE") |
         net_panel(S1, L1, "RE", lim1, wl1, "Resistance: node fill = muscle response, edge width = w_RE") + plot_layout(guides = "collect")) /
        (rev_scatter(tab_m, Z[DN], sprintf("All %d T2D-altered network proteins (Öhman p < 0.05): T2D z vs muscle response", length(DN)), S1) | pc) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(title = "Figure 19a. T2D-altered proteins (Öhman 2021 muscle) in the joint network: their connected pieces, endurance vs resistance",
                  subtitle = "Outline: purple = lower in T2D, orange = higher in T2D; dashed edge = negative weight. Pieces = connected components with >= 2 T2D-altered proteins.",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19a_t2d_small_subgraph.png"), f19a, width = 16, height = 11, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19a_t2d_small_subgraph.png"))

# 19b — story 2
L2 <- layout_of(S2); lim2 <- lim_of(S2); wl2 <- wlim_of(S2)
nd <- function(v, obs, lab, title) ggplot(data.table(v = v), aes(v)) + geom_histogram(bins = 40, fill = "grey75", colour = "grey55") +
  geom_vline(xintercept = obs, colour = "#D7301F", linewidth = 1) + labs(x = lab, y = "random sets", title = title) + theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9))
st2 <- ST[story == "2 connector subgraph"]
f19b <- (net_panel(S2, L2, "EE", lim2, wl2, "Endurance: node fill = muscle response, edge width = w_EE", FALSE) |
         net_panel(S2, L2, "RE", lim2, wl2, "Resistance: node fill = muscle response, edge width = w_RE", FALSE) + plot_layout(guides = "collect")) /
        (nd(null2, length(S2), "largest connected piece (nodes)", sprintf("Size vs %s degree-matched random seed sets (same connector rule): p %.3g", format(N_PERM, big.mark = ","), st2$p[1])) |
         nd(nd2, wd2, "mean w_EE - w_RE", sprintf("Mean w_EE - w_RE vs random connected subgraphs of %d nodes: p %.3g", length(S2), st2$p[2])) |
         rev_scatter(rt2[, .(rev_EE, p_EE, rev_RE, p_RE)], Z[intersect(S2, DN)], "T2D proteins inside the subgraph: T2D z vs muscle response")) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(title = "Figure 19b. T2D-altered proteins plus connector nodes (linked to >= 2 of them): the largest connected subgraph, endurance vs resistance",
                  subtitle = sprintf("%d nodes: %d T2D-altered proteins (outline purple = lower, orange = higher in T2D) and %d connectors (black outline); circles = proteins, triangles = metabolites.",
                                     length(S2), length(intersect(S2, DN)), length(CONN)), theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19b_t2d_connector_subgraph.png"), f19b, width = 17, height = 12, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19b_t2d_connector_subgraph.png"))

# 19c — story 3: forest of reversal tests, Chae replication scatter, consensus proteins per muscle cell
fr <- NT[tissue == "muscle"][, lab := sprintf("%s, %s (n = %d)", fifelse(set == "ohman_2021", "Öhman 2021", "Chae 2018"), nodes, n)]
fr <- rbind(fr[, .(lab, arm = "endurance", v = rev_EE, p = p_EE)], fr[, .(lab, arm = "resistance", v = rev_RE, p = p_RE)], fr[, .(lab, arm = "endurance - resistance", v = diff, p = p_diff)])
fr[, arm := factor(arm, levels = c("endurance", "resistance", "endurance - resistance"))][, lab := factor(lab, levels = rev(unique(lab)))]
pf <- ggplot(fr, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 3, position = position_dodge(0.6)) +
  geom_text(aes(label = sprintf("p %.2g", p)), position = position_dodge(0.6), vjust = -0.9, size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(endurance = "#D7301F", resistance = "#2B8CBE", `endurance - resistance` = "#222222"), name = NULL) +
  labs(x = "muscle reversal (-Spearman of T2D z vs response); difference = endurance - resistance", y = NULL, title = "Reversal per T2D study and node set (muscle; p: T2D z permuted)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
tab_c <- NT[set == "chae_2018" & nodes == "all 471" & tissue == "muscle"]
cz <- setNames(CH[gene %in% GENES471, z], CH[gene %in% GENES471, gene])
zsave <- Z; Z <- cz   # rev_scatter colours by the z it is given; use Chae's z for this panel
pcs <- rev_scatter(tab_c, cz, sprintf("Chae 2018 (independent cohort; %d listed proteins among the 471): T2D z vs muscle response", length(cz)), names(cz)); Z <- zsave
CC[, cell := paste(ome, time)][, cell := factor(cell, levels = c("rna 0.5h", "rna 4h", "rna 24h", "prot 0.5h", "prot 4h", "prot 24h"))]
CC[, gl := sprintf("%s (%s)", node, t2d)][, armlab := ARMLAB[arm]]
pd <- ggplot(CC[!is.na(cell)], aes(cell, gl)) + geom_point(aes(fill = logFC, size = pmin(-log10(adj_p), 4)), shape = 21, colour = "grey30") +
  geom_point(data = CC[!is.na(cell) & adj_p < ALPHA], aes(cell, gl), shape = 21, size = 6.5, colour = "black", stroke = 1.1, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "logFC") + scale_size(range = c(1, 6), name = "-log10 adj. p\n(capped at 4)") +
  facet_wrap(~armlab) + labs(x = "muscle cell (ome, time)", y = NULL, title = sprintf("The %d proteins validated in both T2D studies: muscle response per cell (ring = adj. p < 0.05)", nrow(CONS471))) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9), axis.text.x = element_text(angle = 30, hjust = 1))
f19c <- (pf | pcs) / pd + plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Figure 19c. Replication of the muscle reversal in a second T2D cohort (Chae 2018) and the proteins validated in both studies",
                  subtitle = "Öhman 2021 = primary T2D muscle proteome; Chae 2018 = independent cohort (published significant proteins only). Purple = lower in T2D, orange = higher in T2D.",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(FIG, "19c_t2d_replication.png"), f19c, width = 16, height = 11, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "19c_t2d_replication.png"))

# ---- 8. report ----------------------------------------------------------------------------------------------
f3 <- function(x) formatC(x, digits = 3, format = "fg"); fp <- function(x) formatC(x, digits = 2, format = "g")
row_md <- function(d) paste0("| ", apply(d, 1, paste, collapse = " | "), " |")
md <- c("# Step 19 — endurance vs resistance in the context of T2D (three stories)", "",
  "Disease data: Amar et al. 2024 (Cell Metab 36:1411), T2D muscle sets only (Öhman 2021 primary, Chae 2018 replication); the paper's heart and liver sets",
  "are not in a tissue we measure and are not used. Exercise response = mean normalised log fold change over muscle cells (primary); blood / adipose secondary.",
  "Reversal = -Spearman(T2D z, response): positive = exercise moves T2D-altered proteins opposite to T2D. Permutations: 10,000; seed 20260926.", "",
  "## Node-level tests (all tissues)", "", "| set | nodes | tissue | n | reversal EE | p | reversal RE | p | EE - RE | p (z permuted) | p (arm swap) | q (BH, 3 tissues) |", "|---|---|---|---|---|---|---|---|---|---|---|---|",
  row_md(NT[, .(set, nodes, tissue, n, f3(rev_EE), fp(p_EE), f3(rev_RE), fp(p_RE), f3(diff), fp(p_diff), fp(p_diff_armswap), fp(q_diff))]), "",
  "## Story 1 — small subgraph (figure 19a)", "", "| piece | size | members | edges | mean w_EE | mean w_RE | p (w_EE - w_RE) | reversal EE | reversal RE | p (sign flip) |", "|---|---|---|---|---|---|---|---|---|---|",
  row_md(COMP[, .(component, size, members, n_edges, f3(mean_w_EE), f3(mean_w_RE), fp(p_w_diff), f3(rev_EE), f3(rev_RE), fp(p_signflip))]), "",
  "## Subgraph tests (all stories)", "", "| story | test | observed | null mean | p | test type |", "|---|---|---|---|---|---|",
  row_md(ST[, .(story, test, f3(observed), ifelse(is.na(null_mean), "", f3(null_mean)), ifelse(is.na(p), "", fp(p)), side)]), "",
  sprintf("Story 2 connectors: %s.", paste(CONN, collapse = ", ")), "",
  sprintf("Story 3 consensus proteins among the 471: %s.", paste(sprintf("%s (%s)", CONS471$gene, fifelse(CONS471$oz < 0, "lower", "higher")), collapse = ", ")), "",
  "## Figures (caption skeletons)", "",
  "- **19a** — T2D-altered proteins in the joint network: connected pieces (EE vs RE panels), node-level reversal scatter, per-piece edge weights.",
  "- **19b** — T2D proteins + connectors: largest connected subgraph (EE vs RE panels), size and edge-difference null distributions, reversal inside.",
  "- **19c** — reversal per study and node set (forest), Chae 2018 replication scatter, consensus proteins per muscle cell.")
writeLines(md, file.path(OUT, "reports", "19_t2d_stories.md")); message("-> ", file.path(OUT, "reports", "19_t2d_stories.md"))
