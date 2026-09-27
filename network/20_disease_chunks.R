#!/usr/bin/env Rscript
# =====================================================================================================
# 20_disease_chunks.R — STEP 20 (FIGURES 20a / 20b / 20c): THE BEST ENDURANCE-VS-RESISTANCE STORY PER DISEASE CHUNK
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   How do endurance (EE) and resistance (RE) exercise differ in the context of disease, using disease data from
#   tissues we measure (muscle, adipose)? Three chunks of published human data:
#     OLD T2D  — Amar et al. 2024 (Cell Metab 36:1411) T2D muscle sets: Öhman 2021 (full table), Chae 2018 (significant
#                proteins only). Its three stories are step 19; here it is re-run through the same battery for comparison.
#     NEW T2D  — Kjærgaard et al. 2025 (Cell 188:4106; T2D vs normal glucose tolerance, muscle proteome and
#                phosphoproteome, discovery 77 people + validation 46), Needham et al. 2024 (Cell Metab 36:2542;
#                insulin-resistant vs insulin-sensitive muscle proteome, 19 people), Larsen et al. 2023 (Sci Adv
#                9:eadi7548; subcutaneous adipose proteome, T2D vs lean, 48 men).
#                Plus UK Biobank plasma: Gadd et al. 2024 (Nat Aging 4:1616; 1,468 Olink proteins x incident T2D, 47,600 people).
#     AGEING   — Ubaida-Mohien et al. 2019 (eLife 8:e49874; vastus lateralis proteome, healthy adults aged 20-87;
#                proteins associated with age) and, scaled up, UK Biobank plasma: Sun et al. 2023 (Nature 622:329;
#                Olink Explore 3k, ~54,000 people; per-protein association with age).
#   For every disease set the SAME pre-specified tests are run, and each chunk's story is chosen by a fixed rule.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Disease direction per protein (or phosphosite): log2 change (disease vs healthy; age slope for ageing) and p,
#      turned into a signed z as in step 18. "Altered" = p < 0.05 (sets that publish only significant proteins: all
#      listed; UK Biobank sets, where p < 0.05 flags almost every protein at n ~50,000: the paper's own Bonferroni
#      threshold, 1.7e-5 for Sun 2023 and 3.1e-6 for Gadd 2024). Needham's phosphosite table gives no direction for insulin-resistant vs -sensitive, so only its
#      protein table is used. Larsen's adipose table is regenerated from the authors' own cleaned matrix (GitHub
#      fpm-cbmr/HIIT_adipose_project) with their limma model (group x time, subject correlation); their matrix is
#      already batch-corrected, so the batch term is dropped (same coefficients, slightly optimistic d.f.).
#   2. Exercise response per node: mean normalised log fold change over the cells of the disease set's tissue (muscle
#      or adipose; for the UK Biobank plasma sets, our blood PROTEIN cells, which are the same Olink platform). For phosphosites: the mean MoTrPAC muscle phosphosite logFC over 0.5 / 4 / 24 h, matched on
#      protein + residue.
#   3. Tests per set (10,000 random draws; seed 20260926):
#      A  protein / site level: reversal per arm = -Spearman(disease z, response) over the altered proteins / sites
#         (> 0 = exercise moves them opposite to the disease); EE - RE; T2D z permuted (same shuffle for both arms).
#      B  connector subgraph (protein sets): altered network proteins + nodes linked to >= 2 of them, largest
#         connected piece; its size vs degree-matched random seed sets through the same rule; mean w_EE - w_RE over
#         its edges vs random connected subgraphs of the same size.
#   4. Story per chunk (fixed rule): the test with the smallest EE - RE p among that chunk's primary sets.
#      POST HOC, disclosed: Kjærgaard's discovery result (p 0.039) was not reproduced by their validation cohort, so the
#      two cohorts are pooled per protein (Stouffer's z) and the pooled set is new T2D's primary protein set; a pool of
#      all three T2D muscle cohorts (Öhman + Kjærgaard x 2) is reported as a cross-chunk summary.
#
# HOW TO RUN
#   After steps 14, 17s and 19:   Rscript network/20_disease_chunks.R   (a few minutes; run_all.sh step 20)
#
# DATA AND PROVENANCE
#   $DISEASE_EXT (default ~/Desktop/output/hackathon-2026-track1/external/disease): kjaergaard_2025_cell/mmc1.xlsx
#   (Elsevier supplementary Table S1), needham_2024_cellmetab/mmc3.xlsx (Table S3C), larsen_2023_sciadv/data/
#   Exprs_adipose_clean.txt (authors' GitHub, commit in COMMIT_SHA). $UBAIDA (the Amar et al. repository's copy of
#   Ubaida-Mohien 2019 supplementary data). $DISEASE_SCORES for Öhman / Chae (as step 19). See network/README.md.
#
# TECH STACK:  R 4.4; data.table, readxl, limma, igraph, ggplot2, ggrepel, patchwork.
#
# INPUTS:  the files above; $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 01_nodes_EE.csv, 01_nodes_RE.csv,
#          17_phospho_site_stats.csv, 17_glygen_protein_annotation.csv
# OUTPUTS: $HACK_OUT/20_disease_sets.csv          every disease set: chunk, tissue, level, n measured / altered / in ours
#          $HACK_OUT/20_disease_scores.csv        set x protein / site: log2 change, p, z
#          $HACK_OUT/20_tests.csv                 set x test: reversal EE / RE, difference, p-values (all sets)
#          $HACK_OUT/20_stories.csv               the chosen story per chunk
#          $HACK_OUT/reports/20_disease_chunks.md; $HACK_FIG/20a_all_disease_tests.png, 20b_*.png, 20c_*.png
#
# KNOWN LIMITS
#   Direction matches between acute exercise in healthy adults and a disease (or age) are not evidence of treatment.
#   Disease sets differ in platform, cohort and definition (T2D vs NGT; insulin-resistant vs -sensitive; age slope).
#   Chae 2018 and Ubaida-Mohien 2019 publish only significant proteins, so their tests use pre-selected proteins.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(readxl); library(limma); library(igraph); library(ggplot2); library(ggrepel); library(patchwork) })
# The shared figure style (MoTrPAC landscape-paper look: small Helvetica type, bold lower-case panel tags).
.args <- commandArgs(trailingOnly = FALSE); .here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", .args, value = TRUE))))
source(file.path(.here, "R", "figure_style.R"))

# Folders and inputs (override with environment variables).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
EXT <- Sys.getenv("DISEASE_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/disease"))
UBAIDA <- Sys.getenv("UBAIDA", unset = file.path(Sys.getenv("AMAR_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/amar_2024")), "disease_datasets", "ubaida_mohien_ 2019_elife_stat.csv"))
DISEASE_SCORES <- Sys.getenv("DISEASE_SCORES", unset = file.path(Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network")), "16d_disease_scores.csv.gz"))
dir.create(file.path(OUT, "reports"), recursive = TRUE, showWarnings = FALSE)
SEED <- 20260926; N_PERM <- 10000L; ALPHA <- 0.05
zof <- function(lfc, p) sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)   # signed z, as step 18

# ---- 1. disease sets: one standard table (set, gene, site, logFC, p, z) ------------------------------------
# Old T2D (Amar et al. 2024 sets, step 16d build).
ds <- fread(cmd = sprintf("gzip -dc %s", shQuote(DISEASE_SCORES)))[set %in% c("ohman_2021", "chae_2018") & !is.na(p), .(set, gene, logFC, p)]
OLD <- ds[order(p)][!duplicated(paste(set, gene))][, site := NA_character_]
# Kjærgaard 2025: proteome and phosphoproteome, T2D - NGT, discovery and validation cohorts (Table S1 B-E).
K1 <- file.path(EXT, "kjaergaard_2025_cell", "mmc1.xlsx")
kprot <- function(sheet, set) { x <- as.data.table(read_excel(K1, sheet = sheet, .name_repair = "minimal"))
  gcol <- intersect(c("Gene_name", "Gene_ID"), names(x))[1]
  data.table(set = set, gene = sub(";.*$", "", x[[gcol]]), site = NA_character_, logFC = x[["Group (T2D - NGT) logFC"]], p = x[["Group (T2D - NGT) P.Value"]]) }
kphos <- function(sheet, set) { x <- as.data.table(read_excel(K1, sheet = sheet, .name_repair = "minimal"))
  lc <- grep("^Group.*logFC$", names(x), value = TRUE)[1]; pc <- grep("^Group.*P.Value$", names(x), value = TRUE)[1]
  key <- x$PTM_collapse_key                                        # e.g. "AAAS_S495_M1": gene _ residue _ multiplicity
  data.table(set = set, gene = sub(";.*$", "", x$PG.Genes), site = sub("^.*_([STY][0-9]+)_M[0-9]+$", "\\1", key), mult = sub("^.*_(M[0-9]+)$", "\\1", key),
             logFC = x[[lc]], p = x[[pc]])[mult == "M1"][, mult := NULL] }  # singly phosphorylated form, as MoTrPAC single sites
KJ <- rbind(kprot("C - Proteome DA - discovery", "kjaergaard_2025_prot_discovery"), kprot("E - Proteome DA - validation", "kjaergaard_2025_prot_validation"),
            kphos("B - Phospho DA - discovery", "kjaergaard_2025_phos_discovery"), kphos("D - Phospho DA - validation", "kjaergaard_2025_phos_validation"))
# Needham 2024: proteins, insulin-resistant vs insulin-sensitive (Table S3C).
n3 <- as.data.table(read_excel(file.path(EXT, "needham_2024_cellmetab", "mmc3.xlsx"), sheet = "Table S3C", skip = 1, .name_repair = "minimal"))
NE <- data.table(set = "needham_2024_prot_IR_vs_IS", gene = sub(";.*$", "", n3$GeneNames), site = NA_character_, logFC = n3[["Fold change (log2)"]], p = n3[["P-value"]])
# Larsen 2023: adipose proteome, T2D vs lean before training — the authors' limma model on their cleaned matrix.
LX <- fread(file.path(EXT, "larsen_2023_sciadv", "data", "Exprs_adipose_clean.txt"))
lid <- LX$ID; M <- as.matrix(LX[, !"ID"]); cn <- colnames(M)          # protein IDs are the last column
parts <- tstrsplit(cn, "_"); grp <- paste0(parts[[1]], parts[[2]]); subj <- parts[[4]]
grp <- factor(grp, levels = c("LeanPre", "LeanPost", "ObesePre", "ObesePost", "T2DPre", "T2DPost"))
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
cf <- duplicateCorrelation(M, design, block = subj)
fit <- lmFit(M, design, block = subj, correlation = cf$consensus)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(T2DPre - LeanPre, ObesePre - LeanPre, levels = design)))
tt <- topTable(fit2, coef = "T2DPre - LeanPre", number = Inf, sort.by = "none")
LA <- data.table(set = "larsen_2023_adipose_T2D_vs_lean", gene = sub(";.*$", "", lid), site = NA_character_, logFC = tt$logFC, p = tt$P.Value)
# Ubaida-Mohien 2019: muscle proteins associated with age (published table lists only p < 0.05); age slope = direction.
ub <- fread(UBAIDA)
UB <- data.table(set = "ubaida_mohien_2019_age", gene = ub$GenePrimary, site = NA_character_, logFC = ub$AgeBeta, p = ub$Pvalue)
# UK Biobank plasma (Olink): Sun 2023 Supplementary Table 5 (age beta, SE, log10 p; proteins associated with age, sex
# or BMI) and Gadd 2024 Supplementary Table 4 (Cox hazard ratio per protein for incident type 2 diabetes, age-adjusted).
st5 <- as.data.table(read_excel(file.path(EXT, "sun_2023_nature", "MOESM3.xlsx"), sheet = "ST5", skip = 5, col_names = FALSE, .name_repair = "minimal"))[, 1:5]
setnames(st5, c("id", "name", "beta", "se", "log10p")); st5 <- st5[!is.na(id)][, lapply(.SD, function(v) if (all(grepl("^[-0-9.eE+]+$", na.omit(v)))) as.numeric(v) else v)]
SU <- st5[, .(set = "sun_2023_ukb_age", gene = sub(":.*$", "", id), site = NA_character_, logFC = as.numeric(beta), p = 10^(-as.numeric(log10p)), zz = as.numeric(beta) / as.numeric(se))]
g4 <- as.data.table(read_excel(file.path(EXT, "gadd_2024_nataging", "MOESM3.xlsx"), sheet = "Supplementary Table 4", skip = 8, col_names = FALSE, .name_repair = "minimal"))[, 1:6]
setnames(g4, c("predictor", "outcome", "HR", "LCI", "UCI", "P")); g4 <- g4[outcome == "Type 2 diabetes"][, (3:6) := lapply(.SD, as.numeric), .SDcols = 3:6]
GA <- g4[, .(set = "gadd_2024_ukb_incident_T2D", gene = sub("\\..*$", "", predictor), site = NA_character_, logFC = log(HR), p = P,
             zz = log(HR) / ((log(UCI) - log(LCI)) / (2 * qnorm(0.975))))]   # z from the confidence interval (p underflows to 0)
# All sets, their chunk / tissue / level / whether all or only significant rows are published.
DS <- rbind(OLD, KJ, NE, LA, UB, SU, GA, fill = TRUE)
DS[, `:=`(logFC = suppressWarnings(as.numeric(logFC)), p = suppressWarnings(as.numeric(p)))]   # some supplementary cells are stored as text
DS <- DS[!is.na(gene) & gene != "" & is.finite(logFC) & is.finite(p)]
DS <- DS[order(p)][!duplicated(paste(set, gene, site))][, z := zof(logFC, p)][!is.na(zz), z := zz][, zz := NULL]   # UK Biobank: z from beta / SE
META <- data.table(set = c("ohman_2021", "chae_2018", "kjaergaard_2025_prot_discovery", "kjaergaard_2025_prot_validation", "kjaergaard_2025_phos_discovery",
                           "kjaergaard_2025_phos_validation", "needham_2024_prot_IR_vs_IS", "larsen_2023_adipose_T2D_vs_lean", "ubaida_mohien_2019_age",
                           "gadd_2024_ukb_incident_T2D", "sun_2023_ukb_age"),
                   chunk = c("old T2D", "old T2D", "new T2D", "new T2D", "new T2D", "new T2D", "new T2D", "new T2D", "ageing", "new T2D", "ageing"),
                   label = c("Öhman 2021 T2D muscle protein", "Chae 2018 T2D muscle protein", "Kjærgaard 2025 T2D muscle protein (discovery)", "Kjærgaard 2025 T2D muscle protein (validation)",
                             "Kjærgaard 2025 T2D muscle phosphosite (discovery)", "Kjærgaard 2025 T2D muscle phosphosite (validation)", "Needham 2024 insulin-resistant muscle protein",
                             "Larsen 2023 T2D adipose protein", "Ubaida-Mohien 2019 ageing muscle protein",
                             "UK Biobank (Gadd 2024) incident T2D plasma protein", "UK Biobank (Sun 2023) ageing plasma protein"),
                   tissue = c("muscle", "muscle", "muscle", "muscle", "muscle", "muscle", "muscle", "adipose", "muscle", "blood_prot", "blood_prot"),
                   level = c("protein", "protein", "protein", "protein", "site", "site", "protein", "protein", "protein", "protein", "protein"),
                   sig_only = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
                   alpha = c(rep(ALPHA, 9), 3.1e-6, 1.7e-5),         # "altered" threshold (UK Biobank: the paper's Bonferroni threshold)
                   primary = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
                   validates = NA_character_)
# Pooled sets (POST HOC — added after the Kjærgaard discovery result did not replicate in their validation cohort):
# per protein, Stouffer's z across cohorts that measured it in all of them (same contrast, independent people).
stouffer <- function(sets, name) { x <- DS[set %in% sets, .(z = sum(z) / sqrt(.N), logFC = mean(logFC), k = .N), by = gene][k == length(sets)]
  x[, `:=`(set = name, site = NA_character_, p = 2 * pnorm(-abs(z)))][, .(set, gene, site, logFC, p, z)] }
DS <- rbind(DS, stouffer(c("kjaergaard_2025_prot_discovery", "kjaergaard_2025_prot_validation"), "kjaergaard_2025_prot_pooled"),
                stouffer(c("ohman_2021", "kjaergaard_2025_prot_discovery", "kjaergaard_2025_prot_validation"), "t2d_muscle_3cohort_pooled"))
META <- rbind(META, data.table(set = c("kjaergaard_2025_prot_pooled", "t2d_muscle_3cohort_pooled"), chunk = c("new T2D", "old + new T2D"),
                               label = c("Kjærgaard 2025 T2D muscle protein, both cohorts pooled (post hoc)", "Öhman + Kjærgaard T2D muscle protein, 3 cohorts pooled (post hoc)"),
                               tissue = "muscle", level = "protein", sig_only = FALSE, alpha = ALPHA, primary = c(TRUE, FALSE), validates = NA_character_))
META[set == "kjaergaard_2025_prot_discovery", primary := FALSE]   # the pooled set replaces it as new T2D's primary protein set
DS <- DS[set %in% META$set]
fwrite(DS, file.path(OUT, "20_disease_scores.csv"))
# The sets themselves (chunk, label, readout tissue, level, "altered" threshold, primary), with their row counts.
fwrite(merge(META, DS[, .(n_rows = .N), by = set], by = "set", all.x = TRUE)[order(chunk, set)], file.path(OUT, "20_disease_sets.csv"))

# ---- 2. exercise responses ------------------------------------------------------------------------------------
EEg <- fread(file.path(OUT, "01_nodes_EE.csv")); REg <- fread(file.path(OUT, "01_nodes_RE.csv")); GENES471 <- EEg$gene_symbol
tmean <- function(Mx, tis) { cols <- grep(paste0("^", tis, "_"), names(Mx), value = TRUE); v <- rowMeans(as.matrix(Mx[, ..cols]), na.rm = TRUE); v[is.nan(v)] <- NA; setNames(v, Mx$gene_symbol) }
RESP <- list(muscle = list(EE = tmean(EEg, "muscle"), RE = tmean(REg, "muscle")), adipose = list(EE = tmean(EEg, "adipose"), RE = tmean(REg, "adipose")),
             blood_prot = list(EE = tmean(EEg, "blood_prot"), RE = tmean(REg, "blood_prot")))   # blood protein cells = Olink, as UK Biobank
PH <- fread(file.path(OUT, "17_phospho_site_stats.csv"))[tissue == "muscle" & arm %in% c("EE", "RE")]
SRESP <- dcast(PH[, .(r = mean(logFC, na.rm = TRUE)), by = .(protein, site, arm)], protein + site ~ arm, value.var = "r")   # site response per arm
E <- fread(file.path(OUT, "14_joint_edges.csv")); N <- fread(file.path(OUT, "14_joint_nodes.csv"))
gJ <- graph_from_data_frame(E[, .(node_a, node_b, edge_type, w_EE, w_RE, w_diff)], directed = FALSE, vertices = N[, .(node, node_type)])
PROT <- N[node_type == "protein", node]

# ---- helpers (the same tests as step 19) ----------------------------------------------------------------------
p2 <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (1 + length(null)) }
p1 <- function(obs, null) { null <- null[is.finite(null)]; (1 + sum(null >= obs - 1e-12)) / (1 + length(null)) }
p2c <- function(obs, null) { null <- null[is.finite(null)]; m <- mean(null); (1 + sum(abs(null - m) >= abs(obs - m) - 1e-12)) / (1 + length(null)) }
reversal_test <- function(z, a, b) {
  rz <- rank(z); ra <- rank(a); rb <- rank(b); oA <- -cor(rz, ra); oB <- -cor(rz, rb)
  P <- replicate(N_PERM, sample(rz)); nA <- -drop(cor(ra, P)); nB <- -drop(cor(rb, P))
  data.table(n = length(z), rev_EE = oA, p_EE = p2(oA, nA), rev_RE = oB, p_RE = p2(oB, nB), diff = oA - oB, p_diff = p2(oA - oB, nA - nB)) }
deg_sampler <- function(univ, nodes) { d <- degree(gJ)[univ]; b <- setNames(cut(rank(d, ties.method = "first"), 5, labels = FALSE), univ)
  tb <- table(b[nodes]); function() unlist(lapply(names(tb), function(k) { pool <- univ[b == as.integer(k)]; pool[sample.int(length(pool), tb[[k]])] })) }
nbJ <- lapply(as_adj_list(gJ, mode = "all"), as.integer)
grow <- function(k) { repeat { s <- sample.int(length(nbJ), 1); set <- s; fr <- setdiff(nbJ[[s]], set)
  while (length(set) < k && length(fr)) { x <- fr[sample.int(length(fr), 1)]; set <- c(set, x); fr <- setdiff(union(fr, nbJ[[x]]), set) }
  if (length(set) == k) return(set) } }
mean_edge <- function(v, attr = "w_diff") { w <- edge_attr(induced_subgraph(gJ, v), attr); if (length(w)) mean(w) else NA_real_ }
A <- as_adjacency_matrix(gJ, sparse = TRUE); A@x[] <- 1; VN <- rownames(A)
with_connectors <- function(seeds) { ind <- as.numeric(VN %in% seeds); cnt <- as.numeric(A %*% ind); c(seeds, VN[cnt >= 2 & ind == 0]) }
lcc_nodes <- function(v) { cm <- components(induced_subgraph(gJ, v)); names(cm$membership)[cm$membership == which.max(cm$csize)] }

# ---- 3. tests per set -----------------------------------------------------------------------------------------
TESTS <- list(); SUBG <- list()
for (i in seq_len(nrow(META))) { m <- META[i]; d <- DS[set == m$set]; set.seed(SEED)
  alt <- if (m$sig_only) d else d[p < m$alpha]
  if (m$level == "protein") {
    y <- alt[gene %in% GENES471, .(node = gene, z)][, `:=`(a = RESP[[m$tissue]]$EE[node], b = RESP[[m$tissue]]$RE[node])][is.finite(a) & is.finite(b)]
  } else {
    y <- merge(alt[, .(protein = gene, site, z)], SRESP, by = c("protein", "site"))[, .(node = paste(protein, site), z, a = EE, b = RE)][is.finite(a) & is.finite(b)]
  }
  if (nrow(y) >= 5) TESTS[[m$set]] <- cbind(data.table(set = m$set, test = "A reversal", n_measured_ours = nrow(y)), reversal_test(y$z, y$a, y$b))
  if (m$level == "protein" && m$primary) {                        # B: connector subgraph (muscle network; adipose sets use the same network)
    univ <- intersect(d$gene, PROT); DN <- intersect(alt$gene, PROT)
    if (length(DN) >= 5) { samp <- deg_sampler(univ, DN); S2 <- lcc_nodes(with_connectors(DN))
      nl <- replicate(N_PERM, length(lcc_nodes(with_connectors(samp())))); wd <- mean_edge(S2); nw <- replicate(N_PERM, mean_edge(grow(length(S2))))
      SUBG[[m$set]] <- data.table(set = m$set, test = "B connector subgraph", n_altered_network = length(DN), size = length(S2), p_size = p1(length(S2), nl),
                                  w_EE = mean_edge(S2, "w_EE"), w_RE = mean_edge(S2, "w_RE"), w_diff = wd, p_w_diff = p2c(wd, nw), nodes = paste(S2, collapse = ";")) } }
  message(sprintf("%-38s altered %4d; in ours %4d", m$set, nrow(alt), nrow(y)))
}
TS <- merge(META, rbindlist(TESTS), by = "set", all.x = TRUE); SG <- rbindlist(SUBG)
fwrite(TS, file.path(OUT, "20_tests.csv")); fwrite(SG, file.path(OUT, "20_subgraph_tests.csv"))
print(TS[, .(chunk, set, n_measured_ours, rev_EE = round(rev_EE, 3), p_EE, rev_RE = round(rev_RE, 3), p_RE, diff = round(diff, 3), p_diff)])
print(SG[, .(set, n_altered_network, size, p_size, w_EE = round(w_EE, 3), w_RE = round(w_RE, 3), p_w_diff)])

# ---- 4. the story per chunk -------------------------------------------------------------------------------------
cand <- TS[primary == TRUE & chunk %in% c("old T2D", "new T2D", "ageing") & !is.na(p_diff), .(chunk, set, label, kind = "A reversal (EE - RE)", p = p_diff, stat = diff)]
cand <- rbind(cand, merge(SG[, .(set, kind = "B connector subgraph", p = p_size, stat = size)], META[, .(set, chunk, label)], by = "set")[chunk %in% c("old T2D", "new T2D", "ageing")], fill = TRUE)
STORY <- cand[order(p)][, .SD[1], by = chunk][order(match(chunk, c("old T2D", "new T2D", "ageing")))]
fwrite(STORY, file.path(OUT, "20_stories.csv")); print(STORY)

# ---- 5. figures ---------------------------------------------------------------------------------------------------
COL_LOW <- "#5E3C99"; COL_HIGH <- "#E66100"; ARMCOL <- c(endurance = "#D7301F", resistance = "#2B8CBE", `endurance - resistance` = "#222222")
ARMLAB <- c(EE = "endurance vs control", RE = "resistance vs control")
rE <- RESP$muscle$EE; rR <- RESP$muscle$RE
# 20a — every set and test in one forest (chunks as panels)
fr <- TS[!is.na(p_diff)][, lab := sprintf("%s (n = %d)", label, n_measured_ours)]
fr <- rbind(fr[, .(chunk, lab, arm = "endurance", v = rev_EE, p = p_EE)], fr[, .(chunk, lab, arm = "resistance", v = rev_RE, p = p_RE)], fr[, .(chunk, lab, arm = "endurance - resistance", v = diff, p = p_diff)])
fr[, arm := factor(arm, levels = names(ARMCOL))][, chunk := factor(chunk, levels = c("old T2D", "new T2D", "old + new T2D", "ageing"))]
fr[, lab := factor(lab, levels = rev(unique(fr[order(chunk)]$lab)))]
f20a <- ggplot(fr, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 2.8, position = position_dodge(0.6)) +
  geom_text(aes(label = sprintf("p %.2g", p)), position = position_dodge(0.6), vjust = -0.85, size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = ARMCOL, name = NULL) + facet_grid(chunk ~ ., scales = "free_y", space = "free_y") +
  labs(x = "reversal = -Spearman(disease z, exercise response) over the altered proteins / sites (> 0: opposite to the disease)", y = NULL,
       title = "Fig. 20a | Endurance vs resistance against every disease set (tissue-matched; p: disease z permuted, 10,000 draws)",
       subtitle = "Old T2D = Amar et al. 2024 sets; new T2D = Kjærgaard 2025 (muscle protein and phosphosite), Needham 2024 (insulin-resistant muscle), Larsen 2023 (adipose);\nageing = Ubaida-Mohien 2019. Pooled sets (Stouffer's z per protein across cohorts) are post hoc.") +
  theme_motrpac(8) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom", strip.text.y = element_text(angle = 0))
ggsave(file.path(FIG, "20a_all_disease_tests.png"), f20a, width = 16, height = 10, dpi = 300, bg = "white"); message("-> ", file.path(FIG, "20a_all_disease_tests.png"))

# Story figures (20b new T2D, 20c ageing): the scatter behind the chosen test, and the altered network proteins with their
# links (connector subgraph) in both arms with the minimal PTM tags of step 19.
SITE <- PH[, .(ee = any(adj_p[arm == "EE"] < ALPHA, na.rm = TRUE), re = any(adj_p[arm == "RE"] < ALPHA, na.rm = TRUE)), by = .(protein, site)]
SITE[, cls := fifelse(ee & re, "both", fifelse(ee, "EE", fifelse(re, "RE", "none")))]
GLY <- fread(file.path(OUT, "17_glygen_protein_annotation.csv"))[, .(protein, glyco_N_sites, glyco_O_sites, glyco_OGlcNAc_sites, crosstalk_residues)]
TAGCOL <- c(EE = "#E41A1C", RE = "#377EB8", both = "#984EA3")
tag_data <- function(lay, rpx) { u <- max(diff(range(lay$x)), diff(range(lay$y)), 1e-6) * 0.0032
  rbindlist(lapply(lay$node, function(v) { st <- SITE[protein == v]; g <- GLY[protein == v]
    it <- rbindlist(lapply(names(TAGCOL), function(k) { n <- sum(st$cls == k); if (n) data.table(kind = "P", col = TAGCOL[[k]], n = n) }))
    if (nrow(g)) it <- rbind(it, if (g$glyco_N_sites > 0) data.table(kind = "sq", col = "#0072BC", n = g$glyco_N_sites), if (g$glyco_O_sites > 0) data.table(kind = "sq", col = "#FFD400", n = g$glyco_O_sites),
                             if (g$glyco_OGlcNAc_sites > 0) data.table(kind = "sqo", col = "#0072BC", n = g$glyco_OGlcNAc_sites), if (g$crosstalk_residues > 0) data.table(kind = "star", col = "#FFD700", n = g$crosstalk_residues))
    if (!nrow(it)) return(NULL); xy <- lay[node == v]; r <- rpx * u; a <- (70 - (seq_len(nrow(it)) - 1) * 32) * pi / 180
    it[, `:=`(node = v, x0 = xy$x + r * cos(a), y0 = xy$y + r * sin(a), x1 = xy$x + (r + 7 * u) * cos(a), y1 = xy$y + (r + 7 * u) * sin(a),
              hx = xy$x + (r + 10.5 * u) * cos(a), hy = xy$y + (r + 10.5 * u) * sin(a), nx = xy$x + (r + 16.5 * u) * cos(a), ny = xy$y + (r + 16.5 * u) * sin(a))] })) }
tag_layers <- function(td) { if (is.null(td) || !nrow(td)) return(list()); k <- 0.8; P <- td[kind == "P"]; Q <- td[kind %in% c("sq", "sqo")]; S <- td[kind == "star"]
  list(geom_segment(data = td, aes(x = x0, y = y0, xend = x1, yend = y1), colour = "#555555", linewidth = 0.3),
       geom_point(data = P, aes(hx, hy, colour = col), shape = 16, size = 2.6 * k), geom_point(data = P, aes(hx, hy), shape = 1, size = 2.6 * k, colour = "#222222", stroke = 0.3),
       geom_text(data = P, aes(hx, hy, label = "P"), colour = "white", fontface = "bold", size = 1.45 * k),
       geom_point(data = Q, aes(hx, hy, colour = col), shape = 15, size = 2.3 * k), geom_point(data = Q, aes(hx, hy), shape = 0, size = 2.3 * k, colour = "#222222", stroke = 0.3),
       geom_point(data = Q[kind == "sqo"], aes(hx, hy), shape = 16, size = 0.7 * k, colour = "white"), geom_text(data = S, aes(hx, hy, label = "★"), colour = "#FFD700", size = 3.2 * k),
       geom_text(data = td, aes(nx, ny, label = n), colour = "#222222", size = 1.9 * k)) }
key_panel <- function() { k <- data.table(y = c(7.6, 6.8, 6, 3.9, 3.1, 2.3), lab = c("after endurance only", "after resistance only", "after both arms", "N-linked (GlcNAc)", "O-linked (GalNAc)", "O-GlcNAc"),
                                         col = c(TAGCOL, "#0072BC", "#FFD400", "#0072BC"), kind = c("P", "P", "P", "sq", "sq", "sqo"))
  ggplot(k, aes(0, y)) + geom_point(data = k[kind == "P"], aes(colour = col), shape = 16, size = 3.6) + geom_point(data = k[kind == "P"], shape = 1, size = 3.6, colour = "#222222", stroke = 0.3) +
    geom_text(data = k[kind == "P"], label = "P", colour = "white", fontface = "bold", size = 2) +
    geom_point(data = k[kind != "P"], aes(colour = col), shape = 15, size = 3.3) + geom_point(data = k[kind != "P"], shape = 0, size = 3.3, colour = "#222222", stroke = 0.3) +
    geom_point(data = k[kind == "sqo"], shape = 16, size = 1, colour = "white") + annotate("text", 0, 1.2, label = "★", colour = "#FFD700", size = 4.5) +
    geom_text(aes(0.3, y, label = lab), hjust = 0, size = 2.8) + annotate("text", 0.3, 1.2, label = "phospho = O-glyco residue", hjust = 0, size = 2.8) +
    annotate("text", -0.2, 8.2, label = "PTM tags (number = count)\nMoTrPAC muscle phosphosites\nresponding (adj. p < 0.05, any time)", hjust = 0, vjust = 0, size = 2.7, fontface = "bold") +
    annotate("text", -0.2, 4.5, label = "glycosylation sites, SNFG symbols\n(UniProt via mnet)", hjust = 0, vjust = 0, size = 2.7, fontface = "bold") +
    scale_colour_identity() + coord_cartesian(xlim = c(-0.3, 2.6), ylim = c(0.4, 9.8), clip = "off") + theme_void() }
net_panel <- function(nodes, Zs, alt, arm, lim, wlim, title) {
  lay <- { set.seed(SEED); sg <- induced_subgraph(gJ, nodes); L <- layout_components(sg, layout = layout_with_kk); data.table(node = V(sg)$name, x = L[, 1], y = L[, 2]) }
  es <- E[node_a %in% nodes & node_b %in% nodes][, w := if (arm == "EE") w_EE else w_RE]
  es <- merge(merge(es, lay[, .(node_a = node, xa = x, ya = y)], by = "node_a"), lay[, .(node_b = node, xb = x, yb = y)], by = "node_b")
  nd <- copy(lay)[, r := if (arm == "EE") rE[node] else rR[node]][, outline := fifelse(!node %in% alt, "#222222", fifelse(Zs[node] < 0, COL_LOW, COL_HIGH))][, shape := fifelse(node %in% PROT, 21L, 24L)]
  ggplot() + geom_segment(data = es, aes(x = xa, y = ya, xend = xb, yend = yb, linewidth = abs(w), linetype = w < 0), colour = "grey45", alpha = 0.8) +
    tag_layers(tag_data(lay, 7)) + geom_point(data = nd, aes(x, y, fill = r, colour = outline, shape = shape), size = 3.6, stroke = 1.3) +
    geom_text_repel(data = nd, aes(x, y, label = node, fontface = ifelse(node %in% alt, "bold", "plain")), size = 2.4, seed = SEED, max.time = 60, max.iter = 1e4, box.padding = 0.3, min.segment.length = 0.2, segment.size = 0.2) +
    scale_colour_identity() + scale_shape_identity() + scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = "mean normalised\nresponse (set's tissue)") +
    scale_linewidth(range = c(0.2, 2.6), limits = c(0, wlim), name = "|edge weight|") + coord_equal(clip = "off") + labs(title = title) + theme_motrpac_void(8) +
    theme(plot.title = element_text(face = "bold", size = 10), legend.position = "bottom") }
rev_scatter <- function(tab, zv, title, labs_nodes, xlab, ylab = "mean normalised muscle response") {
  d <- rbind(data.table(node = names(zv), z = zv, r = rE[names(zv)], arm = ARMLAB[["EE"]]), data.table(node = names(zv), z = zv, r = rR[names(zv)], arm = ARMLAB[["RE"]]))[is.finite(r)]
  ann <- data.table(arm = ARMLAB[c("EE", "RE")], lab = c(sprintf("reversal %.2f, p %.3g", tab$rev_EE, tab$p_EE), sprintf("reversal %.2f, p %.3g", tab$rev_RE, tab$p_RE)))
  ggplot(d, aes(z, r)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(aes(colour = z < 0), size = 1.8, alpha = 0.85) + scale_colour_manual(values = c(`TRUE` = COL_LOW, `FALSE` = COL_HIGH), guide = "none") +
    geom_text_repel(data = d[node %in% labs_nodes], aes(label = node), size = 2.3, seed = SEED, max.time = 60, max.iter = 1e4, min.segment.length = 0.1) +
    geom_label(data = ann, aes(x = -Inf, y = -Inf, label = lab), hjust = -0.05, vjust = -0.3, size = 3, label.size = 0, fill = "white", alpha = 0.85, inherit.aes = FALSE) +
    facet_wrap(~arm) + labs(x = xlab, y = ylab, title = title) + theme_motrpac(8) + theme(plot.title = element_text(face = "bold", size = 10)) }
CELL <- fread(file.path(OUT, "17_node_cell_stats.csv"))
story_fig <- function(chunk_name, set_name, file, fig_title, xlab) {
  m <- META[set == set_name]; d <- DS[set == set_name]; alt <- if (m$sig_only) d else d[p < m$alpha]
  rE <<- RESP[[m$tissue]]$EE; rR <<- RESP[[m$tissue]]$RE                      # the readout of this set's tissue
  tlab <- c(muscle = "muscle", adipose = "adipose", blood_prot = "blood protein")[[m$tissue]]
  Zs <- setNames(d$z, d$gene); a471 <- intersect(alt$gene, GENES471); tab <- TS[set == set_name]
  ctis <- if (m$tissue == "blood_prot") "blood" else m$tissue; come <- if (m$tissue == "blood_prot") "prot" else c("rna", "prot")
  sigp <- unique(CELL[node %in% a471 & tissue == ctis & ome %in% come & arm %in% c("EE", "RE") & adj_p < ALPHA, node])
  psc <- rev_scatter(tab, Zs[a471], sprintf("%s: the %d altered proteins among the 471 (labels: a %s cell with adj. p < 0.05)", m$label, length(a471), tlab), sigp, xlab,
                     sprintf("mean normalised %s response", tlab))
  frc <- fr[chunk == chunk_name]
  pfc <- ggplot(frc, aes(v, lab, colour = arm)) + geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 2.8, position = position_dodge(0.6)) +
    geom_text(aes(label = sprintf("p %.2g", p)), position = position_dodge(0.6), vjust = -0.85, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = ARMCOL, name = NULL) + labs(x = "reversal (> 0: opposite to the disease)", y = NULL, title = sprintf("All %s sets", chunk_name)) +
    theme_motrpac(8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
  DN <- intersect(alt$gene, PROT); S2 <- lcc_nodes(with_connectors(DN)); sg <- SG[set == set_name]
  # Rule: a connector subgraph that is not beyond chance and too large to read (> 40 nodes) is replaced by the altered
  # proteins' own connected pieces (edges among altered proteins only, pieces with >= 2 proteins).
  pieces <- nrow(sg) && sg$p_size >= ALPHA && length(S2) > 40
  if (pieces) { cm <- components(induced_subgraph(gJ, DN)); S2 <- names(cm$membership)[cm$membership %in% which(cm$csize >= 2)] }
  big <- length(S2) > 60                                          # still too large to read: keep the 50 most strongly
  if (big) { top <- head(DN[order(-abs(Zs[DN]))], 50)              # altered proteins and their own connected pieces (>= 2)
    cm <- components(induced_subgraph(gJ, top)); S2 <- names(cm$membership)[cm$membership %in% which(cm$csize >= 2)] }
  es <- E[node_a %in% S2 & node_b %in% S2]; wl <- max(abs(c(es$w_EE, es$w_RE)), 1e-6); lim <- max(quantile(abs(c(rE[S2], rR[S2])), 0.95, na.rm = TRUE), 1e-6)
  sub <- if (nrow(sg)) sprintf("Connector subgraph: %d altered network proteins + nodes linked to >= 2 of them; largest piece %d nodes (size vs degree-matched random seeds p %.2g; w_EE %.3f vs w_RE %.3f, difference p %.2g).%s",
                               sg$n_altered_network, sg$size, sg$p_size, sg$w_EE, sg$w_RE, sg$p_w_diff,
                               if (pieces) sprintf(" Not beyond chance and too large to read, so the network shows the altered proteins' own connected pieces%s (%d nodes).", if (big) " (the 50 most strongly altered proteins)" else "", length(S2)) else " The network shows that subgraph.") else ""
  # the network row is ONE tagged panel (both arms + the PTM key), so the key does not take a panel letter
  row <- wrap_elements(full = (net_panel(S2, Zs, alt$gene, "EE", lim, wl, sprintf("Endurance: node fill = %s response, edge width = w_EE", tlab)) |
                               net_panel(S2, Zs, alt$gene, "RE", lim, wl, sprintf("Resistance: node fill = %s response, edge width = w_RE", tlab)) | key_panel()) +
    plot_layout(widths = c(1, 1, 0.26), guides = "collect") & theme(legend.position = "bottom"))
  f <- ((psc | pfc) + plot_layout(widths = c(1.4, 1))) / row + plot_layout(heights = c(1, 1.3)) +
    plot_annotation(tag_levels = "a", title = fig_title, subtitle = paste(strwrap(paste0("Outline / bold: purple = lower, orange = higher in the disease (or with age); black = connector. ", sub), 230), collapse = "\n"), theme = theme(plot.title = element_text(face = "bold", size = 10, family = "Helvetica"), plot.subtitle = element_text(size = 7.5, family = "Helvetica", colour = "grey25")))
  ggsave(file.path(FIG, file), f, width = 19, height = 14, dpi = 300, bg = "white"); message("-> ", file.path(FIG, file)) }
# The chosen story of each chunk (rule above), plus the supporting muscle sets of the story arc.
xl <- function(set) if (grepl("age", set)) "age z (signed; < 0 = lower with age)" else "T2D z (signed; < 0 = lower in T2D / at lower risk)"
story_fig("new T2D", STORY[chunk == "new T2D", set], "20b_new_t2d_story.png",
          sprintf("Fig. 20b | New T2D, chosen story: %s, endurance vs resistance", META[set == STORY[chunk == "new T2D", set], label]), xl(STORY[chunk == "new T2D", set]))
story_fig("ageing", STORY[chunk == "ageing", set], "20c_ageing_story.png",
          sprintf("Fig. 20c | Ageing, chosen story: %s, endurance vs resistance", META[set == STORY[chunk == "ageing", set], label]), xl(STORY[chunk == "ageing", set]))
story_fig("new T2D", "kjaergaard_2025_prot_pooled", "20d_new_t2d_muscle_pooled.png",
          "Fig. 20d | New T2D in muscle (Kjærgaard 2025, both cohorts pooled; post hoc): endurance vs resistance", xl("t2d"))
story_fig("ageing", "ubaida_mohien_2019_age", "20e_ageing_muscle.png",
          "Fig. 20e | Muscle ageing (Ubaida-Mohien 2019): age-associated muscle proteins, endurance vs resistance", xl("age"))

# ---- 6. report --------------------------------------------------------------------------------------------------
f3 <- function(x) formatC(x, digits = 3, format = "fg"); fp <- function(x) formatC(x, digits = 2, format = "g")
row_md <- function(d) paste0("| ", apply(d, 1, paste, collapse = " | "), " |")
md <- c("# Step 20 — the best endurance-vs-resistance story per disease chunk", "",
  "Same tests for every set (tissue-matched; 10,000 draws; seed 20260926). Reversal = -Spearman(disease z, exercise response) over the altered proteins / sites.",
  "Pooled sets (Stouffer's z per protein across cohorts) are POST HOC, added after Kjærgaard's discovery result was not reproduced by their validation cohort.", "",
  "## Protein / site level (A)", "", "| chunk | set | tissue | n | reversal EE | p | reversal RE | p | EE - RE | p |", "|---|---|---|---|---|---|---|---|---|---|",
  row_md(TS[!is.na(p_diff), .(chunk, label, tissue, n_measured_ours, f3(rev_EE), fp(p_EE), f3(rev_RE), fp(p_RE), f3(diff), fp(p_diff))]), "",
  "## Connector subgraph (B)", "", "| set | altered network proteins | largest piece | p (size) | w_EE | w_RE | p (w_EE - w_RE) |", "|---|---|---|---|---|---|---|",
  row_md(merge(SG, META[, .(set, label)], by = "set")[, .(label, n_altered_network, size, fp(p_size), f3(w_EE), f3(w_RE), fp(p_w_diff))]), "",
  "## Chosen story per chunk (smallest p among primary sets)", "", "| chunk | set | test | statistic | p |", "|---|---|---|---|---|",
  row_md(STORY[, .(chunk, label, kind, f3(stat), fp(p))]), "",
  "## Figures (caption skeletons)", "", "- **20a** — every set and test (forest, by chunk).",
  "- **20b / 20c** — the chosen story of new T2D / ageing: scatter (EE vs RE), all sets of the chunk, network in both arms with PTM tags.",
  "- **20d / 20e** — supporting muscle sets: Kjærgaard pooled (new T2D), Ubaida-Mohien (ageing).")
writeLines(md, file.path(OUT, "reports", "20_disease_chunks.md")); message("-> ", file.path(OUT, "reports", "20_disease_chunks.md"))
