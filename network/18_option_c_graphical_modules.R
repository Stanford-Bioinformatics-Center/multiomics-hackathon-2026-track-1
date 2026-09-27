#!/usr/bin/env Rscript
# =====================================================================================================
# 18_option_c_graphical_modules.R — STEP 18c (FIGURE 18c): OPTION C, MODULES AND DISEASE LINKS AS IN AMAR ET AL. 2024
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   A third way to find disease-connected modules, following the methodology of Amar et al. 2024 (Cell Metab
#   36:1411; the MoTrPAC "graphical analysis" it builds on) as closely as our data allow, WITHOUT reusing their
#   code. The network (nodes, edges, layout) is unchanged from our steps 14-15; only module identification and
#   the disease connection follow the paper. Compared with options A (disease filter) and B (overlay on
#   structural modules) from 18_disease_modules.R.
#
# WHAT THIS SCRIPT DOES (plain language), next to what the paper did
#   1. Selection. Paper: training-regulated features (F-test q < 0.05). Here: muscle features (RNA, protein,
#      metabolites; all features, as the paper worked per tissue) with adj. p < 0.05 in any of the 6 exercise cells
#      (endurance or resistance vs control at 0.5 / 4 / 24 h; MoTrPAC delta-delta contrasts). Sex cannot be split
#      (the human results are sex-adjusted), so the two ARMS take the place of the paper's two sexes.
#   2. States with repfdr. Paper: z-scores per sex x time -> repfdr (ztobins, piem) prior over all up / down / null
#      configurations; the most probable configurations kept (190, 86% of the prior mass); each feature assigned
#      its maximum-posterior configuration. Here: the same on the 6 cells (3^6 = 729 configurations), keeping the
#      most probable configurations up to 86% of the prior mass.
#   3. Graphical sets. Paper: node sets (features sharing a state at one time point, e.g. "8w_F1_M1" = up in both
#      sexes at 8 weeks) and edge sets (features sharing a transition between consecutive time points), selected by
#      size (>= 21 features in the saved result). Here: node sets "4h_EE1_RE0" (up after endurance, unchanged after
#      resistance, at 4 h) and edge sets "0.5h_EE1_RE1---4h_EE1_RE0", >= 10 features (smaller data). A module is a
#      set's features that belong to OUR network nodes; modules with >= 5 network nodes are kept.
#   4. Disease connection. Paper: direction concordance of cluster features with each tissue-matched disease set
#      (training up & disease up, training up & disease down, ... per sex), sign test. Here: per module, per arm
#      (the state at the set's time; for an edge set the later state if non-null, else the earlier), the module's
#      genes that are disease-significant (p < 0.05; all listed genes for the significant-only sets Chae, Yuan):
#      concordant (same direction) vs discordant ("exercise opposes disease"), binomial sign test vs 0.5, BH within
#      each disease set. Only T2D muscle (Öhman; Chae = replication) is tissue-matched to our muscle modules; heart
#      and liver sets are reported as exploratory. The overall concordance of all selected features is reported as
#      the background (our addition).
#   5. Network coherence (our addition, as in the Venus week-6 check of the paper's clusters): edges among a module's
#      nodes in our joint network vs 10,000 random node sets of the same size.
#   Figure 18c: module x disease concordance heatmap (per arm) + the top T2D modules drawn on our network.
#
# HOW TO RUN
#   After steps 14, 15, 17s, 18d:   Rscript network/18_option_c_graphical_modules.R   (about 2 minutes; run_all.sh 18c)
#
# DATA AND PROVENANCE
#   Exercise: MotrpacHumanPreSuspensionAnalysis v0.2.4 (MUSCLE_TRNSCRPT_DA, MUSCLE_PROT_PR_DA, MUSCLE_METAB_DA,
#   HUMAN_FEATURE_TO_GENE). Disease: $DISEASE_SCORES (Amar et al. 2024 disease sets, Venus week 6). repfdr 1.2.3.
#   The paper's method was read from its saved graphical-analysis objects (graphical_analysis_results_20220126.RData),
#   not its code.
#
# TECH STACK:  R 4.4; data.table, repfdr, igraph, ggplot2, ggrepel, patchwork, MotrpacHumanPreSuspensionAnalysis.
#
# INPUTS:  $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 15_class_layout.csv, 02_nodes_string.csv,
#          01c_metabolite_ids.csv, 18_approach_comparison.csv; $DISEASE_SCORES
# OUTPUTS: $HACK_OUT/18c_feature_states.csv      selected muscle features: z per cell, assigned configuration
#          $HACK_OUT/18c_modules.csv             module (graphical set), node, feature, node_type
#          $HACK_OUT/18c_module_disease.csv      module x disease set x arm: concordant, discordant, sign-test p, FDR
#          $HACK_OUT/18c_module_summary.csv      per module: name, sizes, network coherence, T2D results
#          $HACK_OUT/18_approach_comparison_ABC.csv  options A / At / B / B' / C side by side
#          $HACK_FIG/18c_graphical_modules.png   (never in the repo)
#
# KNOWN LIMITS
#   Acute human exercise (3 time points, 2 arms) replaces 8-week rat training (4 time points, 2 sexes); the
#   selection rule (any-cell adj. p) replaces the paper's F-test; set-size and prior-mass thresholds are scaled to
#   our smaller data. Modules are defined by shared response patterns, not network connectivity (as in the paper),
#   so they need not be connected on our graph. Direction concordance with healthy-adult exercise is not evidence
#   of treatment.
# =====================================================================================================

# Load packages quietly (the MoTrPAC package is attached: its data are loaded by name).
suppressMessages({ library(data.table); library(repfdr); library(igraph); library(ggplot2); library(ggrepel); library(patchwork)
  library(MotrpacHumanPreSuspensionAnalysis) })
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz"))
SEED <- 20260926; ALPHA <- 0.05; PRIOR_MASS <- 0.86; MIN_SET <- 10L; MIN_NODES <- 5L; N_PERM <- 10000L
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
CELLS <- c("EE_0.5h", "EE_4h", "EE_24h", "RE_0.5h", "RE_4h", "RE_24h")   # the paper's order: one "sex" (arm) at a time, time within

# ---- 1. selected muscle features and their z-scores per cell ----------------------------------------------------
data("HUMAN_FEATURE_TO_GENE", package = "MotrpacHumanPreSuspensionAnalysis")
f2g <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[, .(feature_id = as.character(feature_id), gene = as.character(gene_symbol))])[!is.na(gene)][!duplicated(feature_id)]
da <- rbindlist(lapply(c(rna = "MUSCLE_TRNSCRPT_DA", prot = "MUSCLE_PROT_PR_DA", metab = "MUSCLE_METAB_DA"), function(obj) {
  data(list = obj, package = "MotrpacHumanPreSuspensionAnalysis")
  x <- as.data.table(get(obj))[contrast_category %in% c("EE-CON", "RE-CON") & Timepoint %in% TPS]
  x[, .(feature_id = as.character(feature_id), platform = if ("platform" %in% names(x)) as.character(platform) else "",
        cell = paste0(substr(as.character(contrast_category), 1, 2), "_", names(TPS)[match(as.character(Timepoint), TPS)]), z = z.std, adj_p = adj_p_value)]
}), idcol = "ome")
da[, feature := paste(ome, feature_id, platform, sep = ":")]
Z <- dcast(da, feature + ome + feature_id ~ cell, value.var = "z"); P <- dcast(da, feature ~ cell, value.var = "adj_p")
Z <- Z[complete.cases(Z[, ..CELLS])]
pm <- as.matrix(P[, ..CELLS]); sel <- intersect(P$feature[apply(pm < ALPHA, 1, any, na.rm = TRUE)], Z$feature)
Z <- Z[feature %in% sel]
zmat <- as.matrix(Z[, ..CELLS]); rownames(zmat) <- Z$feature
message(sprintf("selected muscle features (adj. p < 0.05 in any cell): %d (RNA %d, protein %d, metabolites %d)", nrow(Z), Z[ome == "rna", .N], Z[ome == "prot", .N], Z[ome == "metab", .N]))

# ---- 2. repfdr: prior over configurations, most probable configurations, maximum-posterior assignment ---------------
bz <- ztobins(zmat, n.association.status = 3, plot.diagnostics = FALSE)
pe <- piem(bz$pdf.binned.z, bz$binned.z.mat)
Pi <- as.data.table(pe$last.iteration)
cfg_cols <- setdiff(names(Pi), "Pi"); stopifnot(length(cfg_cols) == length(CELLS))
Pi <- Pi[order(-Pi)][, cum := cumsum(Pi)]
keep <- Pi[seq_len(which(cum >= PRIOR_MASS)[1])]
H <- as.matrix(keep[, ..cfg_cols]); pri <- keep$Pi
message(sprintf("configurations kept: %d covering %.0f%% of the prior mass (paper: 190 / 86%%)", nrow(H), 100 * sum(pri)))
# the density of each feature's bin under each state, per cell: pdf.binned.z[cell, bin, state], states ordered -1, 0, 1
pdf <- bz$pdf.binned.z; B <- bz$binned.z.mat
stopifnot(dim(pdf)[1] == length(CELLS), dim(pdf)[3] == 3)
# sanity check of the state order: the "+1" density should sit on the high-z bins
hi_bin <- max(B[, 1]); stopifnot(pdf[1, hi_bin, 3] > pdf[1, hi_bin, 1])
logpost <- sapply(seq_len(nrow(H)), function(k) {
  lp <- log(pri[k]); for (j in seq_along(CELLS)) lp <- lp + log(pmax(pdf[j, B[, j], H[k, j] + 2], 1e-300)); lp })
best <- max.col(logpost, ties.method = "first")
S <- H[best, , drop = FALSE]; colnames(S) <- CELLS
Z[, config := apply(S, 1, paste, collapse = ",")]
Z <- cbind(Z, setnames(as.data.table(S), paste0("s_", CELLS)))

# ---- 3. graphical node sets and edge sets, restricted to our network ---------------------------------------------------
st <- function(t) sprintf("%s_EE%d_RE%d", t, Z[[paste0("s_EE_", t)]], Z[[paste0("s_RE_", t)]])
for (t in names(TPS)) Z[, paste0("node_", t) := st(t)]
Z[, `:=`(edge_1 = paste(node_0.5h, node_4h, sep = "---"), edge_2 = paste(node_4h, node_24h, sep = "---"))]
# map features to our network nodes: genes (RNA / protein) by symbol, metabolites by name
Nn <- fread(file.path(OUT, "14_joint_nodes.csv")); lay <- fread(file.path(OUT, "15_class_layout.csv"))
Z[, node := fifelse(ome == "metab", feature_id, f2g$gene[match(feature_id, f2g$feature_id)])]
Z[, in_network := node %in% Nn$node]
fwrite(Z[, c("feature", "ome", "feature_id", "node", "in_network", CELLS, "config", paste0("node_", names(TPS)), "edge_1", "edge_2"), with = FALSE],
       file.path(OUT, "18c_feature_states.csv"))
long <- rbind(melt(Z[, c("feature", "node", "ome", "in_network", paste0("node_", names(TPS))), with = FALSE], id.vars = c("feature", "node", "ome", "in_network"), value.name = "set"),
              melt(Z[, .(feature, node, ome, in_network, edge_1, edge_2)], id.vars = c("feature", "node", "ome", "in_network"), value.name = "set"))[, variable := NULL]
# sets selected by size on all selected features (as in the paper), excluding the all-null state at every time
sz <- long[, .(n_features = .N), by = set][n_features >= MIN_SET & !grepl("^[^-]*_EE0_RE0$", set) & !grepl("_EE0_RE0---.*_EE0_RE0$", set)]
MODc <- long[set %in% sz$set & in_network == TRUE][, .(features = paste(sort(unique(feature)), collapse = ";"), omes = paste(sort(unique(ome)), collapse = "+")), by = .(set, node)]
kept <- MODc[, .N, by = set][N >= MIN_NODES, set]
MODc <- MODc[set %in% kept][, module := paste0("C:", set)]
MODc[, node_type := Nn$node_type[match(node, Nn$node)]]
fwrite(MODc[, .(module, set, node, node_type, omes, features)], file.path(OUT, "18c_modules.csv"))
message(sprintf("graphical sets >= %d features: %d; with >= %d network nodes (modules): %d", MIN_SET, nrow(sz), MIN_NODES, length(kept)))

# ---- 4. disease connection: direction concordance per module, arm and disease set -----------------------------------------
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))
SETS <- c(ohman_2021 = "T2D muscle (Öhman)", chae_2018 = "T2D muscle (Chae)", coats_2018 = "HCM heart", havlenova_2021 = "HF heart (rat)",
          park_2019 = "MI heart (mouse)", niu_2022_nash = "NASH liver", niu_2022_cirrhosis = "Cirrhosis liver", yuan_2020 = "NAFLD liver", stocks_2022 = "ob/ob liver (mouse)")
MATCHED <- c("ohman_2021", "chae_2018")      # tissue-matched to our muscle modules
dsig <- ds[!is.na(p) & (p < ALPHA | set %in% c("chae_2018", "yuan_2020")), .(set, gene, dis_dir = sign(logFC))][dis_dir != 0][!duplicated(paste(set, gene))]
arm_state <- function(set_name, arm) {   # the module's state for an arm: node set -> that time; edge set -> later state if non-null, else earlier
  parts <- strsplit(set_name, "---")[[1]]
  s <- sapply(parts, function(p) as.integer(sub(paste0(".*_", arm, "(-?[0-9]).*"), "\\1", p)))
  if (length(s) == 2) { if (s[2] != 0) s[2] else s[1] } else s[1]
}
MDc <- rbindlist(lapply(kept, function(sn) { mem <- MODc[set == sn, node]
  rbindlist(lapply(c("EE", "RE"), function(arm) { a <- arm_state(sn, arm)
    rbindlist(lapply(names(SETS), function(s) { g <- dsig[set == s & gene %in% mem]
      if (a == 0 || !nrow(g)) return(data.table(module = paste0("C:", sn), arm = arm, arm_dir = a, set = s, n = nrow(g), concordant = NA_integer_, discordant = NA_integer_, p = NA_real_))
      conc <- sum(g$dis_dir == a); disc <- nrow(g) - conc
      data.table(module = paste0("C:", sn), arm = arm, arm_dir = a, set = s, n = nrow(g), concordant = conc, discordant = disc, p = binom.test(conc, nrow(g), 0.5)$p.value)
    })) })) }))
MDc[, fdr := p.adjust(p, "BH"), by = set]
MDc[, `:=`(tissue_matched = set %in% MATCHED, frac_concordant = concordant / n)]
fwrite(MDc, file.path(OUT, "18c_module_disease.csv"))
# background concordance of all selected network features (our addition)
bgc <- rbindlist(lapply(c("EE", "RE"), function(arm) { x <- Z[in_network == TRUE, .(node, s = rowSums(sign(as.matrix(.SD)))), .SDcols = paste0("s_", arm, "_", names(TPS))][s != 0]
  x <- unique(x[, .(node, dir = sign(s))])[dsig[set == "ohman_2021"], on = c(node = "gene"), nomatch = 0]
  data.table(arm = arm, n = nrow(x), frac_concordant = mean(x$dir == x$dis_dir)) }))

# ---- 5. names and network coherence -------------------------------------------------------------------------------------
E <- fread(file.path(OUT, "14_joint_edges.csv")); g_all <- graph_from_data_frame(E[, .(node_a, node_b)], directed = FALSE, vertices = Nn[, .(node)])
genes <- fread(file.path(OUT, "02_nodes_string.csv"))$gene_symbol; mets <- fread(file.path(OUT, "01c_metabolite_ids.csv"))$metabolite
DBS <- c("REACTOME", "KEGG_MEDICUS", "WP", "PID", "BIOCARTA", "GOBP", "MITOCARTA")
DBL <- c(REACTOME = "Reactome", KEGG_MEDICUS = "KEGG", WP = "WikiPathways", PID = "PID", BIOCARTA = "BioCarta", GOBP = "GO BP", MITOCARTA = "MitoCarta", REFMET = "RefMet")
STOP <- c("of", "by", "to", "the", "in", "and", "via", "for", "on", "a", "an", "or", "with", "from", "into", "at", "as")
pretty_set <- function(set, db) { x <- sub("^(REACTOME|KEGG_MEDICUS|WP|PID|BIOCARTA|GOBP|MITOCARTA|REFMET)_", "", set)
  if (db == "REFMET") return(if (x == "Cer") "Ceramides" else x)
  w <- strsplit(x, "_")[[1]]; keep <- grepl("[0-9]", w) | (nchar(w) <= 4 & !tolower(w) %in% STOP & !grepl("[AEIOU]", substr(w, 2, nchar(w)))) |
    w %in% c("ADME", "MHC", "NAD", "TNF", "RNA", "DNA", "II", "III", "IV", "ER", "ATP", "GTP", "TCA")
  w <- ifelse(keep, w, tolower(w)); o <- paste(w, collapse = " "); paste0(toupper(substr(o, 1, 1)), substring(o, 2)) }
ora1 <- function(m, bg, db) { m <- intersect(m, bg); if (length(m) < 3) return(NULL)
  r <- tryCatch(as.data.table(MotrpacHumanPreSuspensionAnalysis::run_ORA(input = m, background = bg, database = db, min_size = 5L, overlap_cutoff = 0)), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) NULL else r[set_size_in_input >= 2] }
set.seed(SEED)
SUMc <- rbindlist(lapply(kept, function(sn) { mem <- MODc[set == sn]
  r <- rbindlist(list(ora1(mem[node_type == "protein", node], genes, DBS), ora1(mem[node_type == "metabolite", node], mets, "REFMET")), fill = TRUE)
  nm <- if (nrow(r) && any(r$adj_p_value < ALPHA)) r[order(adj_p_value)][1, paste0(pretty_set(as.character(set), as.character(database)), " [", DBL[as.character(database)], "]")] else "no significant pathway"
  v <- intersect(mem$node, V(g_all)$name); obs <- ecount(induced_subgraph(g_all, v))
  null <- replicate(N_PERM, ecount(induced_subgraph(g_all, sample(V(g_all)$name, length(v)))))
  data.table(module = paste0("C:", sn), set = sn, name = nm, n = nrow(mem), n_prot = sum(mem$node_type == "protein"), n_met = sum(mem$node_type == "metabolite"),
             internal_edges = obs, expected_edges = mean(null), coherence_p = (1 + sum(null >= obs)) / (1 + N_PERM))
}))
t2 <- MDc[set == "ohman_2021" & !is.na(p), .(t2d_best_arm = arm[which.min(p)], t2d_n = n[which.min(p)], t2d_concordant = concordant[which.min(p)],
                                             t2d_discordant = discordant[which.min(p)], t2d_p = min(p), t2d_fdr = fdr[which.min(p)]), by = module]
ch <- MDc[set == "chae_2018" & !is.na(p), .(chae_n = sum(n), chae_min_p = min(p), chae_same_direction_as_ohman = NA), by = module]
SUMc <- ch[t2[SUMc, on = "module"], on = "module"]
# does Chae agree with Öhman in the same module and arm (replication of the relation)?
rel <- merge(MDc[set == "ohman_2021" & !is.na(p), .(module, arm, oh_frac = frac_concordant)], MDc[set == "chae_2018" & !is.na(p) & n >= 2, .(module, arm, ch_frac = frac_concordant)], by = c("module", "arm"))
SUMc[, chae_same_direction_as_ohman := sapply(module, function(mo) { r <- rel[module == mo]; if (!nrow(r)) NA else any(sign(r$oh_frac - 0.5) == sign(r$ch_frac - 0.5) & r$oh_frac != 0.5) })]
SUMc[, relation := fcase(!is.na(t2d_p) & t2d_p < ALPHA & t2d_discordant > t2d_concordant, "exercise opposes T2D",
                         !is.na(t2d_p) & t2d_p < ALPHA, "exercise moves with T2D", default = "no significant T2D link")]
fwrite(SUMc[order(t2d_p)], file.path(OUT, "18c_module_summary.csv"))
print(SUMc[order(t2d_p)][1:min(12, .N), .(module, name = substr(name, 1, 40), n, coh = signif(coherence_p, 2), t2d_arm = t2d_best_arm, conc = t2d_concordant, disc = t2d_discordant,
                                          t2d_p = signif(t2d_p, 2), t2d_fdr = signif(t2d_fdr, 2), chae_rep = chae_same_direction_as_ohman, relation)])
print(bgc)

# ---- comparison with options A / At / B / B' ---------------------------------------------------------------------------------
cmpAB <- fread(file.path(OUT, "18_approach_comparison.csv"))
C <- data.table(approach = "C", nodes = uniqueN(MODc$node), modules = length(kept), nodes_in_modules = uniqueN(MODc$node),
                t2d_measured_frac = round(mean(unique(MODc[node_type == "protein", node]) %in% ds[set == "ohman_2021", gene]), 2),
                modules_t2d_sig = SUMc[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N], modules_any_disease_sig = uniqueN(MDc[!is.na(fdr) & fdr < ALPHA, module]),
                modules_chae_replicated = SUMc[!is.na(t2d_fdr) & t2d_fdr < ALPHA & chae_same_direction_as_ohman %in% TRUE, .N],
                chae_direction_agreement = sprintf("%d/%d", rel[sign(oh_frac - .5) == sign(ch_frac - .5), .N], nrow(rel)),
                modules_named = SUMc[name != "no significant pathway", .N], modules_exercise_sig = length(kept),
                exercise_sig_and_muscle_specific = NA_integer_, t2d_sig_and_exercise_sig = SUMc[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N])
CMP <- rbind(cmpAB, C, fill = TRUE)
CMP[, network_coherent_modules := c(rep(NA_character_, nrow(cmpAB)), sprintf("%d/%d", SUMc[coherence_p < ALPHA, .N], nrow(SUMc)))]
fwrite(CMP, file.path(OUT, "18_approach_comparison_ABC.csv")); print(CMP)

# ---- figure 18c: concordance heatmap + the top T2D modules on our network ----------------------------------------------------
hm <- MDc[!is.na(p)][SUMc[, .(module, name)], on = "module", nomatch = 0]
top <- SUMc[!is.na(t2d_p)][order(t2d_p)][1:min(20, .N), module]
hm <- hm[module %in% top]
hm[, `:=`(row = factor(sprintf("%s  %s", sub("^C:", "", module), substr(name, 1, 34)), levels = rev(unique(sprintf("%s  %s", sub("^C:", "", SUMc[match(top, module), module]), substr(SUMc[match(top, module), name], 1, 34))))),
          col = factor(sprintf("%s · %s", SETS[set], arm), levels = as.vector(outer(SETS, c("EE", "RE"), paste, sep = " · "))),
          lab = fifelse(fdr < ALPHA, sprintf("%d/%d*", concordant, n), sprintf("%d/%d", concordant, n)))]
p1 <- ggplot(hm, aes(col, row, fill = frac_concordant)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 1.9) +
  scale_fill_gradient2(low = "#5E3C99", mid = "white", high = "#E66100", midpoint = 0.5, limits = c(0, 1), name = "fraction concordant\n(< 0.5: exercise opposes disease)") +
  # dividers: tissue-matched T2D columns (thin) in each arm block, and between the endurance and resistance blocks (thick)
  geom_vline(xintercept = c(2.5, 11.5), colour = "grey30", linewidth = 0.3) + geom_vline(xintercept = 9.5, colour = "black", linewidth = 0.9) +
  labs(x = NULL, y = NULL, title = "Direction concordance of option-C modules with the Amar et al. 2024 disease sets (per exercise arm)",
       caption = sprintf("Cells: concordant / disease-significant module genes; * sign-test FDR < 0.05 (within set). Thick line: endurance | resistance; left of each thin line: tissue-matched (T2D muscle). Top %d modules by T2D p. Background concordance with Öhman: EE %.2f, RE %.2f.",
                         length(top), bgc[arm == "EE", frac_concordant], bgc[arm == "RE", frac_concordant])) +
  theme_minimal(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 8.5),
                                       plot.caption = element_text(size = 5.5, hjust = 0), legend.position = "right")
# the top T2D modules drawn on the unchanged network (layout of figure 15)
L <- lay[Nn[, .(node, node_type)], on = "node"]; EE_ <- copy(E)[, `:=`(x = L$x[match(node_a, L$node)], y = L$y[match(node_a, L$node)], xend = L$x[match(node_b, L$node)], yend = L$y[match(node_b, L$node)])]
show <- head(SUMc[!is.na(t2d_p)][order(t2d_p), module], 4)
pan <- rbindlist(lapply(show, function(mo) { mem <- MODc[module == mo, node]
  L[, .(node, node_type, x, y, member = node %in% mem, panel = sprintf("%s\n%s (T2D sign-test p %s)", sub("^C:", "", mo), SUMc[module == mo, name], signif(SUMc[module == mo, t2d_p], 2)))] }))
edg <- rbindlist(lapply(unique(pan$panel), function(pn) cbind(EE_, panel = pn)))
p2 <- ggplot() + geom_segment(data = edg, aes(x, y, xend = xend, yend = yend), colour = "grey85", linewidth = 0.15) +
  geom_point(data = pan[member == FALSE], aes(x, y, shape = node_type), size = 0.5, colour = "grey75") +
  geom_point(data = pan[member == TRUE], aes(x, y, shape = node_type), size = 1.6, fill = "#E66100", colour = "grey10", stroke = 0.3) +
  geom_text_repel(data = pan[member == TRUE], aes(x, y, label = node), size = 1.6, seed = SEED, max.time = 60, max.iter = 1e4, max.overlaps = 25, segment.size = 0.1) +
  facet_wrap(~ panel, nrow = 1, labeller = label_wrap_gen(width = 48)) + scale_shape_manual(values = c(protein = 21, metabolite = 24), name = NULL) +
  labs(caption = "Module members (orange) on the unchanged joint network (layout of figure 15); option-C modules are response-pattern sets, so they need not be connected.") +
  theme_void(base_size = 7) + theme(strip.text = element_text(size = 6), legend.position = "none", plot.caption = element_text(size = 5.5, hjust = 0))
q <- p1 / p2 + plot_layout(heights = c(1.5, 1)) + plot_annotation(title = "C · Modules and disease links as in Amar et al. 2024 (repfdr states, graphical node / edge sets, direction concordance) on our network",
  subtitle = sprintf("%d selected muscle features · %d configurations · %d modules (>= %d network nodes) · T2D-linked (FDR < 0.05): %d · network-coherent: %d / %d",
                     nrow(Z), nrow(H), length(kept), MIN_NODES, SUMc[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N], SUMc[coherence_p < ALPHA, .N], nrow(SUMc)),
  theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.5, colour = "grey30")))
ggsave(file.path(FIG, "18c_graphical_modules.png"), q, width = 14, height = 10, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "18c_graphical_modules.png"))
