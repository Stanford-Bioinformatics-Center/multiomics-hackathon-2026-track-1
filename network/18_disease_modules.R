#!/usr/bin/env Rscript
# =====================================================================================================
# 18_disease_modules.R — STEP 18 (FIGURES 18a / 18b): DISEASE-FILTERED VS DISEASE-OVERLAID NETWORK MODULES
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   How should disease data enter the exercise network: (A) FILTER the network to disease-relevant nodes and
#   edges and find modules there, or (B) find modules on the whole network and OVERLAY disease on them? Both are
#   run with all disease datasets of Amar et al. 2024 (Cell Metab 36:1411; doi 10.1016/j.cmet.2023.12.021) and
#   judged against that paper's type 2 diabetes (T2D) muscle data, to propose which approach to use.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Disease data (the paper's 8 proteomics datasets, 9 disease sets, processed in the Venus project week 6:
#   directions median-centred, rodent genes as human orthologs): T2D muscle (Öhman 2021 full table; Chae 2018
#   significant proteins only), HCM heart (Coats 2018), NASH and cirrhosis liver (Niu 2022), NAFLD liver (Yuan
#   2020, significant only), ob/ob mouse liver (Stocks 2022), heart failure rat heart (Havlenova 2021), MI mouse
#   heart (Park 2019). A node's disease score per set: signed z from the p-value (sign from the log2 change).
#     A  FILTER : proteins significant (p < 0.05) in any FULL-distribution set (Öhman, Coats, Niu x2, Stocks,
#                 Havlenova, Park) + their direct neighbours; Louvain modules on that graph (>= 5 members).
#     B  OVERLAY: Louvain modules on the whole joint network (disease not used), then disease overlaid.
#     B' OVERLAY on the strength-trimmed graph (edge-strength sum >= median; the earlier "strength first" idea).
#   For every module: its disease score in each of the 9 sets (mean signed z of measured members; permutation p
#   vs 10,000 random sets of the same size from the same graph), a pathway name (ORA vs our universe), and its
#   muscle exercise response (MoTrPAC CAMERA-PR; EE / RE vs control, 0.5 / 4 / 24 h) with a tissue-swap control
#   (the same module in blood and adipose).
#   Evaluation against the paper's T2D data. Öhman 2021 scores T2D direction. Chae 2018 is HELD OUT: never used for
#   filtering, it is an independent replication (module members over-represented among Chae's T2D proteins, and
#   the direction agreement of Öhman and Chae for shared members). Comparison per approach: modules, coverage,
#   T2D-scored modules, Chae replication, exercise-responsive modules and their muscle specificity.
#   Figures (static): 18a approach A, 18b approach B: the network with named module bubbles (node fill = T2D
#   change, outline = number of disease sets significant) + a module x disease heatmap.
#
# HOW TO RUN
#   After steps 14, 15, 17s:   Rscript network/18_disease_modules.R   (about 2 minutes; run_all.sh step 18d)
#
# DATA AND PROVENANCE
#   Disease: $DISEASE_SCORES (default ~/Desktop/output/week_6/_shared/disease_scores.csv.gz, built by the Venus
#   project week-6 code from the Amar et al. 2024 disease inputs; fingerprinted by step 0). Exercise:
#   MotrpacHumanPreSuspensionAnalysis v0.2.4. Network: steps 14 / 15 (team mnet resource, STRING v12 >= 700).
#
# TECH STACK:  R 4.4; data.table, igraph, ggplot2, ggrepel, ggforce, patchwork, MotrpacHumanPreSuspensionAnalysis + TMSig.
#
# INPUTS:  $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 15_class_layout.csv, 02_nodes_string.csv,
#          01c_metabolite_ids.csv; $DISEASE_SCORES
# OUTPUTS: $HACK_OUT/18_disease_node_scores.csv  node x disease set: log2 change, p, signed z
#          $HACK_OUT/18_disease_modules.csv      approach, module, node, node_type
#          $HACK_OUT/18_module_disease.csv       module x disease set: members measured, mean z, permutation p
#          $HACK_OUT/18_module_summary.csv       per module: name, size, T2D (Öhman), Chae replication, exercise
#          $HACK_OUT/18_approach_comparison.csv  per approach: the evaluation numbers
#          $HACK_OUT/18_module_camera.csv, 18_module_ora.csv; $HACK_OUT/reports/18_disease_modules.md
#          $HACK_FIG/18a_disease_filter_modules.png, 18b_disease_overlay_modules.png
#
# KNOWN LIMITS
#   Only the T2D sets are muscle; heart and liver disease sets are tissue-mismatched to our exercise data (muscle,
#   adipose, blood), so they inform relevance, not exercise direction. Öhman's p is a 4-group ANOVA p. Few module
#   members are measured per set, so module tests are underpowered. Filtering on disease significance makes
#   approach A's own disease scores partly circular; the within-graph permutation and the held-out Chae set
#   are there to judge it fairly. Direction matches between healthy-adult exercise and disease are not evidence
#   of treatment.
# =====================================================================================================

# Load packages quietly (the MoTrPAC package is attached: its DA loader looks tables up by name).
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel); library(ggforce); library(patchwork)
  library(MotrpacHumanPreSuspensionAnalysis) })

# Folders and inputs (override with environment variables).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz"))
dir.create(file.path(OUT, "reports"), recursive = TRUE, showWarnings = FALSE)
# Fixed seeds and thresholds.
SEED <- 20260926; MIN_SIZE <- 5L; ALPHA <- 0.05; N_PERM <- 10000L
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
# Disease sets: short labels (heatmap columns), in a fixed order; which can filter (full distribution).
SETS <- c(ohman_2021 = "T2D muscle (Öhman)", chae_2018 = "T2D muscle (Chae, held out)", coats_2018 = "HCM heart", havlenova_2021 = "HF heart (rat)",
          park_2019 = "MI heart (mouse)", niu_2022_nash = "NASH liver", niu_2022_cirrhosis = "Cirrhosis liver", yuan_2020 = "NAFLD liver (sig. only)",
          stocks_2022 = "ob/ob liver (mouse)")
FILTER_SETS <- c("ohman_2021", "coats_2018", "havlenova_2021", "park_2019", "niu_2022_nash", "niu_2022_cirrhosis", "stocks_2022")

# ---- disease scores per node --------------------------------------------------------------------------------
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))
stopifnot(all(names(SETS) %in% ds$set))
ds <- ds[!is.na(p), .(set, gene, logFC, p)][order(p)][!duplicated(paste(set, gene))]
ds[, z := sign(logFC) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)]
E <- fread(file.path(OUT, "14_joint_edges.csv")); Nn <- fread(file.path(OUT, "14_joint_nodes.csv"))
N <- fread(file.path(OUT, "15_class_layout.csv"))[, .(node, x, y)][Nn, on = "node"]
N[, strength := pmax(strength_EE, strength_RE)]
D <- ds[gene %in% N$node]
fwrite(D[, .(node = gene, set, logFC, p, z)], file.path(OUT, "18_disease_node_scores.csv"))
# per node: T2D (Öhman) change for colouring; number of full-distribution sets in which it is significant
N <- D[set == "ohman_2021", .(node = gene, t2d_logFC = logFC)][N, on = "node"]
N <- D[set %in% FILTER_SETS & p < ALPHA, .(n_disease = uniqueN(set)), by = .(node = gene)][N, on = "node"]
N[is.na(n_disease), n_disease := 0L]
chae <- unique(ds[set == "chae_2018", gene])

# ---- the three graphs and their modules --------------------------------------------------------------------------
# filter seeds: significant in any full-distribution disease set (A), or in the T2D muscle set only (At, At-strict)
seed_p <- N[node_type == "protein" & n_disease > 0, node]
t2d_p <- D[set == "ohman_2021" & p < ALPHA & gene %in% N[node_type == "protein", node], gene]
with_nb <- function(s) unique(c(s, E[node_a %in% s, node_b], E[node_b %in% s, node_a]))
GRAPHS <- list(A = with_nb(seed_p), At = with_nb(t2d_p), As = t2d_p, B = N$node, Bp = N[strength >= median(N$strength), node])
ALAB <- c(A = "A  filter: significant in any disease set + neighbours", At = "At filter: T2D-significant (Öhman) + neighbours",
          As = "As filter: only edges between T2D-significant proteins", B = "B  overlay: whole joint network", Bp = "B' overlay: strength-trimmed network")
MOD <- rbindlist(lapply(names(GRAPHS), function(a) {
  e <- E[node_a %in% GRAPHS[[a]] & node_b %in% GRAPHS[[a]]]
  g <- graph_from_data_frame(e[, .(node_a, node_b)], directed = FALSE)
  set.seed(SEED); cl <- cluster_louvain(g)
  m <- data.table(approach = a, node = V(g)$name, community = as.integer(membership(cl)))
  sz <- m[, .N, by = community][N >= MIN_SIZE][order(-N)]
  m <- m[community %in% sz$community][, module := sprintf("%s_M%02d", a, match(community, sz$community))][, community := NULL]
  m[, node_type := N$node_type[match(node, N$node)]]
}))
stopifnot(all(MOD$node %in% N$node))
fwrite(MOD, file.path(OUT, "18_disease_modules.csv"))
message(paste(sprintf("%s: %d nodes, %d modules", names(GRAPHS), lengths(GRAPHS), sapply(names(GRAPHS), function(a) uniqueN(MOD[approach == a, module]))), collapse = "; "))

# ---- disease overlay: module x set (permutation within each graph) -----------------------------------------------------
set.seed(SEED)
MD <- rbindlist(lapply(split(MOD, MOD$module), function(m) {
  pool_nodes <- GRAPHS[[m$approach[1]]]
  rbindlist(lapply(names(SETS), function(s) {
    ps <- D[set == s & gene %in% pool_nodes]; z <- ps[gene %in% m$node, z]
    if (s %in% c("chae_2018", "yuan_2020") || length(z) < 3)   # significant-only sets: no score (hit counts instead)
      return(data.table(module = m$module[1], set = s, n_measured = length(z), mean_z = if (length(z)) mean(z) else NA_real_, perm_p = NA_real_))
    obs <- mean(z); null <- replicate(N_PERM, mean(sample(ps$z, length(z))))
    data.table(module = m$module[1], set = s, n_measured = length(z), mean_z = obs, perm_p = (1 + sum(abs(null) >= abs(obs))) / (1 + N_PERM))
  }))
}))
fwrite(MD, file.path(OUT, "18_module_disease.csv"))

# ---- held-out T2D replication (Chae 2018) --------------------------------------------------------------------------------
# enrichment: are the module's proteins over-represented among Chae's T2D proteins (hypergeometric, background =
# the graph's proteins)? agreement: for members in both Öhman and Chae, do the two T2D directions agree?
oh <- D[set == "ohman_2021"]; ch <- ds[set == "chae_2018"]
REP <- rbindlist(lapply(split(MOD, MOD$module), function(m) {
  bg <- intersect(GRAPHS[[m$approach[1]]], N[node_type == "protein", node]); mem <- intersect(m$node, bg)
  k <- sum(mem %in% chae); K <- sum(bg %in% chae)
  both <- intersect(intersect(mem, oh$gene), ch$gene)
  data.table(module = m$module[1], chae_hits = k, chae_expected = length(mem) * K / length(bg),
             chae_enrich_p = phyper(k - 1, K, length(bg) - K, length(mem), lower.tail = FALSE),
             chae_ohman_shared = length(both), chae_ohman_agree = sum(sign(oh$logFC[match(both, oh$gene)]) == sign(ch$logFC[match(both, ch$gene)])))
}))

# ---- module names (ORA vs our universe) ----------------------------------------------------------------------------------
genes <- fread(file.path(OUT, "02_nodes_string.csv"))$gene_symbol; mets <- fread(file.path(OUT, "01c_metabolite_ids.csv"))$metabolite
DBS <- c("REACTOME", "KEGG_MEDICUS", "WP", "PID", "BIOCARTA", "GOBP", "MITOCARTA")
DB_LABEL <- c(REACTOME = "Reactome", KEGG_MEDICUS = "KEGG", WP = "WikiPathways", PID = "PID", BIOCARTA = "BioCarta", GOBP = "GO BP", MITOCARTA = "MitoCarta", REFMET = "RefMet")
STOP <- c("of", "by", "to", "the", "in", "and", "via", "for", "on", "a", "an", "or", "with", "from", "into", "at", "as")
pretty_set <- function(set, db) {
  x <- sub("^(REACTOME|KEGG_MEDICUS|WP|PID|BIOCARTA|GOBP|MITOCARTA|REFMET)_", "", set)
  if (db == "REFMET") return(if (x == "Cer") "Ceramides" else x)
  w <- strsplit(x, "_")[[1]]
  keep <- grepl("[0-9]", w) | (nchar(w) <= 4 & !tolower(w) %in% STOP & !grepl("[AEIOU]", substr(w, 2, nchar(w)))) |
    w %in% c("ADME", "MHC", "NAD", "TNF", "RNA", "DNA", "II", "III", "IV", "ER", "ATP", "GTP", "TCA")
  w <- ifelse(keep, w, tolower(w)); out <- paste(w, collapse = " "); paste0(toupper(substr(out, 1, 1)), substring(out, 2))
}
ora_one <- function(members, background, db) {
  members <- intersect(members, background); if (length(members) < 3) return(NULL)
  r <- tryCatch(as.data.table(MotrpacHumanPreSuspensionAnalysis::run_ORA(input = members, background = background, database = db, min_size = MIN_SIZE, overlap_cutoff = 0)),
                error = function(e) NULL)
  if (is.null(r) || !nrow(r)) NULL else r[set_size_in_input >= 2]
}
ORA <- rbindlist(lapply(split(MOD, MOD$module), function(m) {
  r <- rbindlist(list(ora_one(m[node_type == "protein", node], genes, DBS), ora_one(m[node_type == "metabolite", node], mets, "REFMET")), fill = TRUE)
  if (!nrow(r)) NULL else r[, module := m$module[1]]
}), fill = TRUE)
ORA[, label := paste0(mapply(pretty_set, as.character(set), as.character(database)), " [", DB_LABEL[as.character(database)], "]")]
fwrite(ORA[order(module, adj_p_value), .(module, database = as.character(database), set = as.character(set), label, overlap_n = set_size_in_input, p_value, adj_p_value)],
       file.path(OUT, "18_module_ora.csv"))

# ---- exercise response per module (MoTrPAC CAMERA-PR, genes), muscle + tissue-swap control ------------------------------
gmt <- file.path(OUT, "18_disease_modules.gmt")
gp <- MOD[node_type == "protein", .(node = list(node)), by = module][lengths(node) >= MIN_SIZE]
writeLines(sapply(seq_len(nrow(gp)), function(i) paste(c(gp$module[i], "disease module (genes)", gp$node[[i]]), collapse = "\t")), gmt)
cam <- as.data.table(MotrpacHumanPreSuspensionAnalysis::run_cameraPR(selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ol"),
  selected_tissues = "all", path_to_gmt = gmt, min_size = MIN_SIZE, overlap_cutoff = 0.7))[contrast_type == "exercise_with_controls"]
cam[, cs := as.character(contrast_short)]
cam[, `:=`(arm = fifelse(grepl("^Endur", cs), "EE", "RE"), time = names(TPS)[sapply(cs, function(s) which(sapply(TPS, grepl, x = s))[1])],
           ome = fifelse(assay == "transcript-rna-seq", "rna", "prot"), module = as.character(set), tissue = as.character(tissue))]
cam[, approach := sub("_M[0-9]+$", "", module)]   # A, At, As, B, Bp
cam[, fdr := p.adjust(p_value, "BH"), by = .(approach, tissue, ome, arm, time)]
CAM <- cam[, .(module, tissue, ome, arm, time, n = set_size, direction = as.character(direction), z = z.std, p = p_value, fdr)]
fwrite(CAM, file.path(OUT, "18_module_camera.csv"))
EX <- CAM[, .(muscle_sig = sum(tissue == "muscle" & fdr < ALPHA), muscle_cells = sum(tissue == "muscle"),
              swap_sig = sum(tissue != "muscle" & fdr < ALPHA), swap_cells = sum(tissue != "muscle"),
              best_muscle = .SD[tissue == "muscle"][order(fdr)][1, sprintf("%s %s %s %s (FDR %.2g)", arm, ome, time, fifelse(z > 0, "up", "down"), fdr)],
              best_muscle_z = .SD[tissue == "muscle"][order(fdr)][1, z]), by = module]

# ---- per-module summary and the approach comparison --------------------------------------------------------------------------
nm <- ORA[adj_p_value < ALPHA][order(adj_p_value)][, .SD[1], by = module][, .(module, name = label)]
SUM <- MOD[, .(approach = approach[1], n = .N, n_prot = sum(node_type == "protein"), n_met = sum(node_type == "metabolite")), by = module]
SUM <- Reduce(function(a, b) b[a, on = "module"], list(SUM, nm, MD[set == "ohman_2021", .(module, t2d_n = n_measured, t2d_mean_z = mean_z, t2d_perm_p = perm_p)], REP, EX))
SUM[is.na(name), name := "no significant pathway"]
SUM[, n_disease_sets_sig := sapply(module, function(mo) MD[module == mo & !is.na(perm_p) & perm_p < ALPHA, .N])]
SUM[, `:=`(t2d_sig = !is.na(t2d_perm_p) & t2d_perm_p < ALPHA, chae_rep = chae_hits >= 2 & chae_enrich_p < ALPHA,
           ex_sig = !is.na(muscle_sig) & muscle_sig > 0, muscle_specific = !is.na(muscle_sig) & muscle_sig / pmax(muscle_cells, 1) > swap_sig / pmax(swap_cells, 1))]
fwrite(SUM[order(approach, module)], file.path(OUT, "18_module_summary.csv"))
CMP <- SUM[, .(nodes = sapply(approach[1], function(a) length(GRAPHS[[a]])), modules = .N, nodes_in_modules = sum(n),
               t2d_measured_frac = round(sum(t2d_n, na.rm = TRUE) / sum(n_prot), 2),
               modules_t2d_sig = sum(t2d_sig), modules_any_disease_sig = sum(n_disease_sets_sig > 0),
               modules_chae_replicated = sum(chae_rep), chae_direction_agreement = sprintf("%d/%d", sum(chae_ohman_agree), sum(chae_ohman_shared)),
               modules_named = sum(name != "no significant pathway"), modules_exercise_sig = sum(ex_sig),
               exercise_sig_and_muscle_specific = sum(ex_sig & muscle_specific),
               t2d_sig_and_exercise_sig = sum(t2d_sig & ex_sig)), by = approach][order(approach)]
fwrite(CMP, file.path(OUT, "18_approach_comparison.csv"))
print(CMP)

# ---- figures 18a (filter) and 18b (overlay): network + module x disease heatmap ----------------------------------------------
lim_t <- as.numeric(quantile(abs(N$t2d_logFC), 0.95, na.rm = TRUE))
fig <- function(a, file, title) {
  nodes <- N[node %in% GRAPHS[[a]]]; e <- E[node_a %in% nodes$node & node_b %in% nodes$node]
  e[, `:=`(x = nodes$x[match(node_a, nodes$node)], y = nodes$y[match(node_a, nodes$node)], xend = nodes$x[match(node_b, nodes$node)],
           yend = nodes$y[match(node_b, nodes$node)], w = pmax(abs(w_EE), abs(w_RE)))]
  m <- MOD[approach == a][nodes[, .(node, x, y)], on = "node", nomatch = 0]
  ml <- m[, .(x = mean(x), y = max(y) + 0.03), by = module][SUM[, .(module, name)], on = "module", nomatch = 0][, lab := sub("^[A-Za-z]+_", "", module)]
  net <- ggplot() +
    geom_mark_hull(data = m, aes(x, y, group = module), colour = "grey45", fill = "grey60", alpha = 0.07, linewidth = 0.25, expand = unit(2.2, "mm"), radius = unit(2, "mm"), concavity = 3) +
    geom_segment(data = e, aes(x, y, xend = xend, yend = yend, linewidth = w), colour = "grey60", alpha = 0.7) +
    geom_point(data = nodes, aes(x, y, fill = t2d_logFC, shape = node_type, size = pmin(n_disease, 4)), colour = "grey10", stroke = 0.3) +
    geom_text_repel(data = ml, aes(x, y, label = lab), size = 2.4, fontface = "bold", colour = "grey20", seed = SEED, max.time = 60, max.iter = 1e4, box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.15) +
    scale_fill_gradient2(low = "#5E3C99", mid = "white", high = "#E66100", midpoint = 0, limits = c(-lim_t, lim_t), oob = scales::squish, na.value = "grey88",
                         name = "T2D vs NGT, log2 (Öhman; grey = not measured)", guide = guide_colourbar(barwidth = unit(4.5, "cm"), barheight = unit(0.22, "cm"), title.position = "top")) +
    scale_shape_manual(values = c(protein = 21, metabolite = 24), name = NULL) +
    scale_size(range = c(1, 3.6), breaks = 0:4, labels = c("0", "1", "2", "3", "4+"), name = "disease sets significant (p < 0.05)") +
    scale_linewidth(range = c(0.1, 1), guide = "none") + coord_cartesian(clip = "off") +
    theme_void(base_size = 8) + theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 6.5), legend.text = element_text(size = 6))
  # heatmap: modules (rows, with names) x disease sets; colour = mean signed z; star = permutation p < 0.05; Chae column = replication
  md <- MD[module %in% ml$module]
  md <- merge(md, SUM[, .(module, name, chae_enrich_p, chae_hits)], by = "module")
  md[, `:=`(row = sprintf("%s  %s", sub("^[A-Za-z]+_", "", module), substr(name, 1, 44)), col = factor(SETS[set], levels = SETS))]
  md[set %in% c("chae_2018", "yuan_2020"), mean_z := NA_real_]
  md[set == "chae_2018", perm_p := chae_enrich_p]
  md[, star := fifelse(!is.na(perm_p) & perm_p < ALPHA, "*", "")]
  md[set == "chae_2018", star := fifelse(chae_hits >= 2 & chae_enrich_p < ALPHA, paste0(chae_hits, "R"), fifelse(chae_hits > 0, as.character(chae_hits), ""))]
  md[set == "yuan_2020", star := fifelse(n_measured > 0, as.character(n_measured), "")]
  md[, row := factor(row, levels = rev(sort(unique(row))))]
  hm <- ggplot(md, aes(col, row, fill = mean_z)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 2.4) +
    scale_fill_gradient2(low = "#5E3C99", mid = "white", high = "#E66100", midpoint = 0, limits = c(-3, 3), oob = scales::squish, na.value = "grey93",
                         name = "module mean signed z\n(< 0 lower, > 0 higher in disease)") +
    labs(x = NULL, y = NULL, caption = paste0("* permutation p < 0.05 (vs same-size sets of the same graph). Significant-only sets (Chae, held out; NAFLD):\n",
                                         "number of module proteins in the published list; R = over-represented (hypergeometric p < 0.05).")) +
    theme_minimal(base_size = 7) + theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), legend.position = "bottom",
                                         plot.caption = element_text(size = 5.5, hjust = 0))
  cm <- CMP[approach == a]
  q <- net + hm + plot_layout(widths = c(1.25, 1)) + plot_annotation(title = title,
    subtitle = sprintf("%d nodes, %d modules · modules T2D-scored significant (Öhman): %d · replicated in held-out Chae: %d · exercise-responsive in muscle: %d (muscle-specific %d)",
                       cm$nodes, cm$modules, cm$modules_t2d_sig, cm$modules_chae_replicated, cm$modules_exercise_sig, cm$exercise_sig_and_muscle_specific),
    theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.5, colour = "grey30")))
  ggsave(file.path(FIG, file), q, width = 14, height = 8, dpi = 300, bg = "white"); message("-> ", file.path(FIG, file))
}
fig("At", "18a_disease_filter_modules.png", "A · Disease filter first: modules of the network trimmed to T2D-significant muscle proteins (Öhman 2021, in Amar et al. 2024) and their neighbours")
fig("B", "18b_disease_overlay_modules.png", "B · Modules first, disease overlaid: structural modules of the whole joint network scored against the Amar et al. 2024 disease sets")

# ---- report ------------------------------------------------------------------------------------------------------------------
fmt <- function(v) formatC(v, digits = 2, format = "g")
rep <- c("# Disease-filtered vs disease-overlaid modules (step 18, figures 18a / 18b)", "",
         sprintf("Generated by `network/18_disease_modules.R` on %s from the 8 disease datasets (9 sets) of Amar et al. 2024 (Cell Metab 36:1411).", format(Sys.Date())), "",
         "## Approach comparison", "", paste0("| metric | ", paste(ALAB[CMP$approach], collapse = " | "), " |"), paste0("|---|", strrep("---|", nrow(CMP))),
         sapply(setdiff(names(CMP), "approach"), function(k) sprintf("| %s | %s |", k, paste(CMP[[k]], collapse = " | "))), "",
         "## Modules", "", "| module | name | n | T2D n | T2D mean z | T2D perm p | Chae hits (exp.) | Chae p | disease sets sig | best muscle response | swap sig |", "|---|---|---|---|---|---|---|---|---|---|---|",
         SUM[order(approach, module), sprintf("| %s | %s | %d | %s | %s | %s | %d (%.1f) | %s | %d | %s | %d/%d |", module, name, n, t2d_n, fmt(t2d_mean_z), fmt(t2d_perm_p),
                                             chae_hits, chae_expected, fmt(chae_enrich_p), n_disease_sets_sig, best_muscle, swap_sig, swap_cells)])
writeLines(rep, file.path(OUT, "reports", "18_disease_modules.md")); message("-> ", file.path(OUT, "reports", "18_disease_modules.md"))
