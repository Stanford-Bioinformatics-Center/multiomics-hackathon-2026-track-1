#!/usr/bin/env Rscript
# =====================================================================================================
# 16d_disease_scores.R — STEP 16d: THE DISEASE SCORES (AMAR ET AL. 2024 DISEASE SETS), BUILT IN THIS REPOSITORY
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Which proteins go up or down in each published disease proteome, on one footing? Steps 17i, 18d, 18c, 18bc, 19
#   and 20 read the answer, $HACK_OUT/16d_disease_scores.csv.gz. Until 2026-09-27 they read a copy built outside this
#   repository (the Venus project, week 6: Vidal_code/week_6/00.1_inputs.R, section 2); this step is that code, moved
#   here unchanged, so the hackathon pipeline builds its own inputs from public sources. It reproduces the old file
#   exactly (checked when it was moved: every row and value identical).
#
# WHAT THIS SCRIPT DOES (plain language)
#   For each disease set (the 9 of Amar et al. 2024, Cell Metab 36:1411, doi 10.1016/j.cmet.2023.12.021, plus Yuan
#   2020 as an overlay), from the authors' published tables (downloaded by step 16f, pinned commit):
#     - per gene: the disease-vs-healthy log fold change and the two-sided p, as each paper reports it (or recomputed
#       from the per-sample values where the published p is unusable: Niu 2022, Welch t-tests);
#     - one row per gene (the first in the publishers' own order, which does not look at the statistic);
#     - direction read MEDIAN-CENTRED: logFC minus the set's median (published group means can carry a global offset,
#       e.g. Öhman's median protein is 20% lower in T2D from total-intensity normalisation; most proteins are assumed
#       not to change). The raw logFC is kept (logFC_raw); p is untouched. Overlays (significant proteins only) are
#       not centred, since their median is not a baseline;
#     - the one-sided p for "down in disease" and "up in disease": p/2 if the effect goes that way, 1 - p/2 if not.
#   Gene symbols: where a table lists several, the first one STRING v12 knows (Niu lists e.g. ";MYCBP"); Coats 2018
#   genes from the UniProt map the Amar et al. repository ships (its description text loses 58 genes).
#
# HOW TO RUN
#   Rscript network/16d_disease_scores.R          (run_all.sh step 16d, after 16f)
#
# INPUTS  ($AMAR_EXT, from step 16f): disease_datasets/*.csv / *.tsv (Amar et al. repository @ COMMIT_SHA),
#         9606.protein.info.v12.0.txt.gz (STRING v12)
# OUTPUTS ($HACK_OUT):
#   16d_disease_scores.csv.gz   set, disease, species, tissue, role, gene, logFC (median-centred), p, note,
#                               logFC_raw, p_down, p_up — one row per set x gene
#   16d_disease_inventory.csv   per set: genes, share p < 0.05, up / down, median offset, signs flipped by centring
#
# TECH STACK: R, data.table.   RUN TIME: a few seconds.
# =====================================================================================================
suppressPackageStartupMessages(library(data.table))
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
AMAR_EXT <- Sys.getenv("AMAR_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/amar_2024"))
DISEASE_DIR <- file.path(AMAR_EXT, "disease_datasets")
if (!dir.exists(DISEASE_DIR)) stop("missing ", DISEASE_DIR, " (run: bash network/16f_fetch_amar_disease_sets.sh)")

say <- function(...) cat(sprintf(...), "\n", sep = "")
canon <- function(x) toupper(trimws(as.character(x)))
one_sided <- function(p, eff, want) {             # want = +1 (up) or -1 (down)
  p <- pmin(pmax(p, .Machine$double.xmin), 1)
  fifelse(sign(eff) == want, p / 2, 1 - p / 2) }

# ---- from here to the end: Vidal_code/week_6/00.1_inputs.R section 2, unchanged (paths aside) -------------------
rd <- function(f, ...) fread(file.path(DISEASE_DIR, f), ...)
mk <- function(set, disease, species, tissue, role, gene, logFC, p, note) {
  d <- data.table(set, disease, species, tissue, role, gene = canon(gene), logFC, p, note)
  d <- d[!is.na(gene) & gene != "" & is.finite(logFC) & is.finite(p)]
  d[, p := pmin(pmax(p, .Machine$double.xmin), 1)]
  ## one row per gene: first occurrence in file order (the publishers' own ordering), which does not look at the
  ## statistic.
  d <- d[!duplicated(gene)]
  ## MEDIAN-CENTRED DIRECTION (see the header). Raw logFC is kept; p is untouched. Overlays list significant proteins
  ## only, so their median is not a baseline.
  d[, logFC_raw := logFC]
  if (role != "overlay") d[, logFC := logFC - median(logFC)]
  d[, `:=`(p_down = one_sided(p, logFC, -1), p_up = one_sided(p, logFC, 1))]
  d }
## first listed symbol that STRING knows (Niu lists e.g. ";MYCBP", "HDGFRP3;HDGFL3"); falls back to the first token
STRING_SYM <- canon(fread(file.path(AMAR_EXT, "9606.protein.info.v12.0.txt.gz"), showProgress = FALSE)[[2]])
first_sym <- function(x) vapply(strsplit(trimws(as.character(x)), "[ ;]+"), function(t) {
  t <- t[nzchar(t)]; if (!length(t)) return(NA_character_)
  k <- t[canon(t) %chin% STRING_SYM]; if (length(k)) k[1] else t[1] }, "")

DZ <- list()
## Ohman 2021 -- human T2D, vastus lateralis. 4-group ANOVA p over T2D/IGT/IFG/NGT; direction from the T2D vs NGT
## group means (the only per-group values published).
o <- rd("ohman_2021_iscience.csv")
DZ$ohman <- mk("ohman_2021", "T2D", "human", "SKM-VL", "primary", o[["Gene name"]],
               log2(o$SWATH_T2D / o$SWATH_NGT), as.numeric(o[["ANOVA_p-value"]]),
               "4-group ANOVA p; sign of log2(T2D/NGT) means")

## Coats 2018 -- human HCM, heart. Full distribution with per-sample intensities; direction and logFC from the group
## means. Gene symbol from the UniProt entry via the map the Amar et al. repository ships (the description's "GN" tag
## is missing / truncated / outdated for 58 of 527 genes); description parsing is the fallback only.
co <- rd("coats_2018_circgpm_stat.csv", check.names = TRUE)
cs <- grep("^Control[.][0-9]+$", names(co), value = TRUE); hs <- grep("^HCM[.][0-9]+$", names(co), value = TRUE)
cmap <- fread(file.path(DISEASE_DIR, "uniprot2gene_name_coats_bg.tsv"))
co_gene <- cmap$To[match(co$UniProt.Accession, cmap$From)]
co_gene[is.na(co_gene)] <- sub("^.* GN ([^ ]+) .*$", "\\1", co$Description[is.na(co_gene)])
co_lfc <- log2(rowMeans(co[, ..hs], na.rm = TRUE) / rowMeans(co[, ..cs], na.rm = TRUE))
DZ$coats <- mk("coats_2018", "HCM", "human", "HEART", "primary", co_gene, co_lfc,
               as.numeric(co[["Anova..p."]]), "published ANOVA p (2 groups); log2 ratio of group means")

## Niu 2022 -- human NASH and cirrhosis, liver. Recomputed from per-sample log2 intensities (45 samples, no missing
## values): Welch t-test of each disease group vs healthy controls (the published p is rounded to 2 dp and its F-test
## pools three groups).
ni <- rd("niu_2022_molsystbiol_set.csv", check.names = FALSE)
grp <- names(ni); X <- as.matrix(ni[, 3:(ncol(ni) - 3)]); g <- grp[3:(ncol(ni) - 3)]
stopifnot(setequal(unique(g), c("Cirrhosis", "Control", "NASH")))
welch <- function(a, b) {
  res <- t(apply(X, 1, function(r) { t <- t.test(r[g == a], r[g == b]); c(t$p.value, mean(r[g == a]) - mean(r[g == b])) }))
  list(p = res[, 1], lfc = res[, 2]) }
for (dis in c("NASH", "Cirrhosis")) {
  w <- welch(dis, "Control")
  DZ[[paste0("niu_", tolower(dis))]] <- mk(paste0("niu_2022_", tolower(dis)), dis, "human", "LIVER",
    "primary", first_sym(ni[["PG.Genes"]]), w$lfc, w$p,
    sprintf("recomputed Welch t, %s (n=%d) vs control (n=%d), log2 intensities", dis, sum(g == dis), sum(g == "Control")))
}

## Stocks 2022 -- mouse ob/ob vs lean, liver (secondary). Full t-test table.
st <- rd("stocks_2022_molcellproteomics_stat.csv", check.names = TRUE)
DZ$stocks <- mk("stocks_2022", "ob/ob", "mouse", "LIVER", "secondary", first_sym(st$Gene.names),
                as.numeric(st$Student.s.T.test.Difference.ob_lean),
                10^(-as.numeric(st$X.Log.Student.s.T.test.p.value.ob_lean)), "published t-test; difference ob - lean")

## Havlenova 2021 -- rat HF right ventricle (secondary). Only ADJUSTED p published.
hv <- rd("havlenova_2021_scirep.csv", check.names = TRUE)
DZ$havlenova <- mk("havlenova_2021", "HF (RV)", "rat", "HEART", "secondary", hv$Gene.Symbol,
                   log2(as.numeric(hv$RV.Fold.change)), as.numeric(hv$RV.Adj..P.Value),
                   "ADJUSTED p only (no raw p published); BUM fit is approximate")

## Park 2019 -- mouse MI 8 w vs sham, heart (secondary). The description field contains unquoted commas, so parse by
## position: the first five fields and the last three are fixed.
pk <- readLines(file.path(DISEASE_DIR, "park_2019_celldeathdis.csv"), encoding = "UTF-8")[-1]
pk <- rbindlist(lapply(strsplit(pk, ","), function(f)
  data.table(gene = f[2], fc = as.numeric(f[3]), p = as.numeric(f[5]))))
DZ$park <- mk("park_2019", "MI 8w", "mouse", "HEART", "secondary", first_sym(pk$gene),
              log2(pk$fc), pk$p, "published p; FC MI/sham; parsed by position")

## Overlays only: significant proteins published, no full distribution.
ch <- rd("chae_2018_emm_stat.csv", check.names = TRUE)
DZ$chae <- mk("chae_2018", "T2D", "human", "SKM", "overlay", ch$Symbol,
              as.numeric(ch$log2.Fold.change..T2DM...NGT.), as.numeric(ch$FDR),
              "SIGNIFICANT PROTEINS ONLY (FDR<0.1); overlay, never a network side")
yu <- rd("yuan_2020_jproteomics.csv", check.names = TRUE)
DZ$yuan <- mk("yuan_2020", "NAFLD", "human", "LIVER", "overlay", yu$GENE.name,
              log2(as.numeric(yu$Fold.changed..NAFLD.MHO.)), as.numeric(yu$p.value),
              "SIGNIFICANT PROTEINS ONLY (p<0.05); overlay, never a network side")

DS <- rbindlist(DZ)
INV <- DS[, .(genes = .N, frac_p05 = round(mean(p < 0.05), 3), up = sum(logFC > 0), down = sum(logFC < 0),
              median_offset = round(median(logFC_raw), 3),
              sig_flipped_by_centring = sum(p < 0.05 & sign(logFC) != sign(logFC_raw)),
              p_max = signif(max(p), 3), note = note[1]), by = .(set, disease, species, tissue, role)]
say("\n===== DISEASE INPUTS ====="); print(INV[, !"note"])
## sanity: a real full distribution has p reaching ~1; an overlay does not
stopifnot(INV[role != "overlay", all(p_max > 0.9)])
fwrite(DS, file.path(OUT, "16d_disease_scores.csv.gz"))
fwrite(INV, file.path(OUT, "16d_disease_inventory.csv"))
say("written: 16d_disease_scores.csv.gz (%d rows, %d sets), 16d_disease_inventory.csv", nrow(DS), uniqueN(DS$set))
