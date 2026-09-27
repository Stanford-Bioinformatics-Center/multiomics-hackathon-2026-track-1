#!/usr/bin/env Rscript
# =====================================================================================================
# 18_option_bc_hybrid.R — STEP 18bc (FIGURE 18d): OPTION B MODULES + OPTION C DISEASE TEST
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   The proposed hybrid: modules defined from OUR network without disease data (option B: structural Louvain
#   modules of the whole joint network, as in figure 18b), connected to disease with the test of Amar et al. 2024
#   (option C: direction concordance of exercise responses with disease changes, per exercise arm, sign test).
#
# WHAT THIS SCRIPT DOES (plain language)
#   Modules: option B's modules (18_disease_modules.csv, approach "B"; Louvain on the whole joint network, seed fixed,
#   >= 5 members; disease not used). Exercise direction of each member feature (option C's repfdr states,
#   18c_feature_states.csv; muscle RNA / protein / metabolite features selected at adj. p < 0.05 in any cell): per
#   arm, "up" or "down" if its states at 0.5 / 4 / 24 h are non-null and agree in sign; features that go both ways
#   in an arm, or are null throughout, are left out. Disease connection (as the paper and option C): for each
#   module, arm and disease set, the module's directional features whose gene is disease-significant (p < 0.05; all
#   listed genes for the significant-only sets Chae, Yuan) are concordant (same direction) or discordant ("exercise
#   opposes disease"); binomial sign test vs 0.5; BH within each disease set. The overall concordance of all
#   directional network features is reported as the background (our addition). T2D muscle is the only
#   tissue-matched set; heart and liver sets are exploratory.
#   Figure 18d: the network with option-B module bubbles + a module x (disease set, arm) concordance heatmap.
#
# HOW TO RUN
#   After steps 18d and 18c:   Rscript network/18_option_bc_hybrid.R   (seconds; run_all.sh step 18bc)
#
# DATA AND PROVENANCE:  as options B (18_disease_modules.R) and C (18_option_c_graphical_modules.R).
# TECH STACK:  R 4.4; data.table, ggplot2, ggrepel, ggforce, patchwork.
# INPUTS:   $HACK_OUT/18_disease_modules.csv, 18_module_summary.csv, 18c_feature_states.csv, 14_joint_edges.csv,
#           14_joint_nodes.csv, 15_class_layout.csv, 18_approach_comparison_ABC.csv; $DISEASE_SCORES
# OUTPUTS:  $HACK_OUT/18bc_module_disease.csv  module x disease set x arm: concordant, discordant, sign-test p, FDR
#           $HACK_OUT/18bc_module_summary.csv  per module: name, directional features, T2D (Öhman) concordance, Chae
#           $HACK_OUT/18_approach_comparison_ABCD.csv  options A / At / B / B' / C / B+C side by side
#           $HACK_FIG/18d_structural_modules_concordance.png (never in the repo)
# KNOWN LIMITS: only features selected as exercise-responsive in muscle carry a direction, so each module's test
#   uses a subset of its members; direction concordance with healthy-adult exercise is not evidence of treatment.
# =====================================================================================================

suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(ggforce); library(patchwork) })
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz"))
SEED <- 20260926; ALPHA <- 0.05
SETS <- c(ohman_2021 = "T2D muscle (Öhman)", chae_2018 = "T2D muscle (Chae)", coats_2018 = "HCM heart", havlenova_2021 = "HF heart (rat)",
          park_2019 = "MI heart (mouse)", niu_2022_nash = "NASH liver", niu_2022_cirrhosis = "Cirrhosis liver", yuan_2020 = "NAFLD liver", stocks_2022 = "ob/ob liver (mouse)")
MATCHED <- c("ohman_2021", "chae_2018")

# ---- option B modules and option C feature directions -------------------------------------------------------------
MOD <- fread(file.path(OUT, "18_disease_modules.csv"))[approach == "B"]
SUMB <- fread(file.path(OUT, "18_module_summary.csv"))[approach == "B", .(module, name)]
FS <- fread(file.path(OUT, "18c_feature_states.csv"))[in_network == TRUE]
arm_dir <- function(arm) { s <- as.matrix(FS[, paste0(arm, "_", c("0.5h", "4h", "24h")), with = FALSE])
  apply(s, 1, function(v) { v <- v[v != 0]; if (!length(v) || length(unique(sign(v))) > 1) 0L else as.integer(sign(v[1])) }) }
FS[, `:=`(dir_EE = arm_dir("EE"), dir_RE = arm_dir("RE"))]
# the directional features of each module's members
MF <- FS[, .(feature, node, ome, dir_EE, dir_RE)][MOD[, .(module, node)], on = "node", nomatch = 0, allow.cartesian = TRUE]
message(sprintf("option B modules: %d; member nodes with a selected muscle feature: %d; directional features: EE %d, RE %d",
                uniqueN(MOD$module), uniqueN(MF$node), MF[dir_EE != 0, .N], MF[dir_RE != 0, .N]))

# ---- option C disease test on option B modules ----------------------------------------------------------------------
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))
dsig <- ds[!is.na(p) & (p < ALPHA | set %in% c("chae_2018", "yuan_2020")), .(set, gene, dis_dir = sign(logFC))][dis_dir != 0][!duplicated(paste(set, gene))]
MD <- rbindlist(lapply(unique(MOD$module), function(mo) rbindlist(lapply(c("EE", "RE"), function(arm) {
  f <- MF[module == mo][get(paste0("dir_", arm)) != 0, .(node, d = get(paste0("dir_", arm)))]
  rbindlist(lapply(names(SETS), function(s) { x <- f[dsig[set == s], on = c(node = "gene"), nomatch = 0]
    if (!nrow(x)) return(data.table(module = mo, arm = arm, set = s, n = 0L, concordant = NA_integer_, discordant = NA_integer_, p = NA_real_))
    conc <- sum(x$d == x$dis_dir)
    data.table(module = mo, arm = arm, set = s, n = nrow(x), concordant = conc, discordant = nrow(x) - conc, p = binom.test(conc, nrow(x), 0.5)$p.value) }))
}))))
MD[, fdr := p.adjust(p, "BH"), by = set]
MD[, `:=`(frac_concordant = concordant / n, tissue_matched = set %in% MATCHED)]
fwrite(MD, file.path(OUT, "18bc_module_disease.csv"))
# background: all directional network features vs Öhman (our addition)
bg <- rbindlist(lapply(c("EE", "RE"), function(arm) { x <- FS[get(paste0("dir_", arm)) != 0, .(node, d = get(paste0("dir_", arm)))][dsig[set == "ohman_2021"], on = c(node = "gene"), nomatch = 0]
  data.table(arm = arm, n = nrow(x), concordant = sum(x$d == x$dis_dir), frac = mean(x$d == x$dis_dir), p = binom.test(sum(x$d == x$dis_dir), nrow(x), 0.5)$p.value) }))
print(bg)

# ---- per-module summary and comparison ------------------------------------------------------------------------------------
t2 <- MD[set == "ohman_2021" & !is.na(p)][order(p), .SD[1], by = module][, .(module, t2d_arm = arm, t2d_n = n, t2d_concordant = concordant, t2d_discordant = discordant, t2d_p = p, t2d_fdr = fdr)]
chr <- merge(MD[set == "ohman_2021" & !is.na(p), .(module, arm, oh = frac_concordant)], MD[set == "chae_2018" & !is.na(p) & n >= 2, .(module, arm, ch = frac_concordant)], by = c("module", "arm"))
SUM <- MOD[, .(n = .N, n_prot = sum(node_type == "protein"), n_met = sum(node_type == "metabolite")), by = module]
SUM <- Reduce(function(a, b) b[a, on = "module"], list(SUM, SUMB, MF[, .(directional_features_EE = sum(dir_EE != 0), directional_features_RE = sum(dir_RE != 0)), by = module], t2))
SUM[, chae_same_direction := sapply(module, function(mo) { r <- chr[module == mo]; if (!nrow(r)) NA else any(sign(r$oh - .5) == sign(r$ch - .5) & r$oh != .5) })]
SUM[, n_disease_sets_fdr := sapply(module, function(mo) MD[module == mo & !is.na(fdr) & fdr < ALPHA, uniqueN(set)])]
SUM[, relation := fcase(!is.na(t2d_fdr) & t2d_fdr < ALPHA & t2d_discordant > t2d_concordant, "exercise opposes T2D",
                        !is.na(t2d_fdr) & t2d_fdr < ALPHA, "exercise moves with T2D",
                        !is.na(t2d_p) & t2d_p < ALPHA & t2d_discordant > t2d_concordant, "exercise opposes T2D (nominal p only)",
                        !is.na(t2d_p) & t2d_p < ALPHA, "exercise moves with T2D (nominal p only)", default = "no T2D link")]
fwrite(SUM[order(t2d_p)], file.path(OUT, "18bc_module_summary.csv"))
print(SUM[order(t2d_p), .(module, name = substr(name, 1, 40), n, dirEE = directional_features_EE, dirRE = directional_features_RE, arm = t2d_arm,
                          conc = t2d_concordant, disc = t2d_discordant, t2d_p = signif(t2d_p, 2), fdr = signif(t2d_fdr, 2), chae = chae_same_direction, sets_fdr = n_disease_sets_fdr, relation)])
cmp <- fread(file.path(OUT, "18_approach_comparison_ABC.csv"))
BC <- data.table(approach = "B+C", nodes = uniqueN(MOD$node), modules = uniqueN(MOD$module), nodes_in_modules = nrow(MOD),
                 t2d_measured_frac = cmp[approach == "B", t2d_measured_frac], modules_t2d_sig = SUM[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N],
                 modules_any_disease_sig = SUM[n_disease_sets_fdr > 0, .N], modules_chae_replicated = SUM[!is.na(t2d_fdr) & t2d_fdr < ALPHA & chae_same_direction %in% TRUE, .N],
                 chae_direction_agreement = sprintf("%d/%d", chr[sign(oh - .5) == sign(ch - .5), .N], nrow(chr)),
                 modules_named = SUM[name != "no significant pathway", .N], modules_exercise_sig = cmp[approach == "B", modules_exercise_sig],
                 exercise_sig_and_muscle_specific = cmp[approach == "B", exercise_sig_and_muscle_specific], t2d_sig_and_exercise_sig = SUM[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N],
                 network_coherent_modules = "by construction")
CMP <- rbind(cmp, BC, fill = TRUE); CMP[, modules_t2d_nominal := c(rep(NA_integer_, nrow(cmp)), SUM[!is.na(t2d_p) & t2d_p < ALPHA, .N])]
fwrite(CMP, file.path(OUT, "18_approach_comparison_ABCD.csv")); print(CMP[, .(approach, modules, modules_t2d_sig, modules_t2d_nominal, modules_any_disease_sig, modules_chae_replicated, chae_direction_agreement)])

# ---- figure 18d --------------------------------------------------------------------------------------------------------------
E <- fread(file.path(OUT, "14_joint_edges.csv")); Nn <- fread(file.path(OUT, "14_joint_nodes.csv"))
N <- fread(file.path(OUT, "15_class_layout.csv"))[, .(node, x, y)][Nn, on = "node"]
oh <- ds[set == "ohman_2021", .(node = gene, t2d_logFC = logFC)][!duplicated(node)]; N <- oh[N, on = "node"]
N[, has_dir := node %in% FS[dir_EE != 0 | dir_RE != 0, node]]
e <- copy(E)[, `:=`(x = N$x[match(node_a, N$node)], y = N$y[match(node_a, N$node)], xend = N$x[match(node_b, N$node)], yend = N$y[match(node_b, N$node)], w = pmax(abs(w_EE), abs(w_RE)))]
m <- MOD[N[, .(node, x, y)], on = "node", nomatch = 0]
ml <- m[, .(x = mean(x), y = max(y) + 0.03), by = module][, lab := sub("^B_", "", module)]
lim_t <- as.numeric(quantile(abs(N$t2d_logFC), 0.95, na.rm = TRUE))
net <- ggplot() +
  geom_mark_hull(data = m, aes(x, y, group = module), colour = "grey45", fill = "grey60", alpha = 0.07, linewidth = 0.25, expand = unit(2.2, "mm"), radius = unit(2, "mm"), concavity = 3) +
  geom_segment(data = e, aes(x, y, xend = xend, yend = yend, linewidth = w), colour = "grey60", alpha = 0.7) +
  geom_point(data = N, aes(x, y, fill = t2d_logFC, shape = node_type, size = has_dir), colour = "grey10", stroke = 0.3) +
  geom_text_repel(data = ml, aes(x, y, label = lab), size = 2.4, fontface = "bold", colour = "grey20", seed = SEED, max.time = 60, max.iter = 1e4, box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.15) +
  scale_fill_gradient2(low = "#5E3C99", mid = "white", high = "#E66100", midpoint = 0, limits = c(-lim_t, lim_t), oob = scales::squish, na.value = "grey88",
                       name = "T2D vs NGT, log2 (Öhman; grey = not measured)", guide = guide_colourbar(barwidth = unit(4.5, "cm"), barheight = unit(0.22, "cm"), title.position = "top")) +
  scale_shape_manual(values = c(protein = 21, metabolite = 24), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.6), labels = c(`FALSE` = "no directional exercise feature", `TRUE` = "directional exercise feature (option C states)"), name = NULL) +
  scale_linewidth(range = c(0.1, 1), guide = "none") + coord_cartesian(clip = "off") +
  theme_void(base_size = 8) + theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 6.5), legend.text = element_text(size = 6))
hm <- MD[SUM[, .(module, name)], on = "module"]
hm[, `:=`(row = sprintf("%s  %s", sub("^B_", "", module), substr(name, 1, 40)), col = factor(sprintf("%s · %s", SETS[set], arm), levels = as.vector(outer(SETS, c("EE", "RE"), paste, sep = " · "))),
          lab = fifelse(is.na(n) | n == 0, "", fifelse(!is.na(fdr) & fdr < ALPHA, sprintf("%d/%d*", concordant, n), sprintf("%d/%d", concordant, n))))]
hm[, row := factor(row, levels = rev(sort(unique(row))))]
hmp <- ggplot(hm[n > 0], aes(col, row, fill = frac_concordant)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 1.9) +
  geom_vline(xintercept = c(2.5, 11.5), colour = "grey30", linewidth = 0.3) + geom_vline(xintercept = 9.5, colour = "black", linewidth = 0.9) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient2(low = "#5E3C99", mid = "white", high = "#E66100", midpoint = 0.5, limits = c(0, 1), name = "fraction concordant\n(< 0.5: exercise opposes disease)") +
  labs(x = NULL, y = NULL, caption = sprintf(paste0("Cells: concordant / directional features of disease-significant module genes; * sign-test FDR < 0.05 (within set). Thick line: endurance | resistance;\n",
                                                    "left of each thin line: tissue-matched (T2D muscle). Background concordance with Öhman over all directional network features: EE %d/%d (%.2f), RE %d/%d (%.2f)."),
                                             bg[arm == "EE", concordant], bg[arm == "EE", n], bg[arm == "EE", frac], bg[arm == "RE", concordant], bg[arm == "RE", n], bg[arm == "RE", frac])) +
  theme_minimal(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), legend.position = "bottom", plot.caption = element_text(size = 5.5, hjust = 0))
q <- net + hmp + plot_layout(widths = c(1, 1.15)) + plot_annotation(
  title = "B + C · Structural modules of the whole joint network (option B), linked to disease with the Amar et al. 2024 direction-concordance test (option C)",
  subtitle = sprintf("%d modules · T2D-linked (Öhman, FDR < 0.05): %d (nominal p < 0.05: %d) · modules with any disease link at FDR < 0.05: %d · Chae agrees with Öhman in %s module-arms",
                     uniqueN(MOD$module), SUM[!is.na(t2d_fdr) & t2d_fdr < ALPHA, .N], SUM[!is.na(t2d_p) & t2d_p < ALPHA, .N], SUM[n_disease_sets_fdr > 0, .N], BC$chae_direction_agreement),
  theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.5, colour = "grey30")))
ggsave(file.path(FIG, "18d_structural_modules_concordance.png"), q, width = 15, height = 8.5, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "18d_structural_modules_concordance.png"))
