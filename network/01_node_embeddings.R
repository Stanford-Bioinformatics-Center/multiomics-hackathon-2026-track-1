#!/usr/bin/env Rscript
# 01_node_embeddings.R — node table for the exercise network.
#
# Universe: the 471 Entrez genes measured in ALL of
#   transcriptomics (adipose, blood, muscle) + proteomics (adipose MS, blood OLINK, muscle MS),
#   both arms, every post-exercise timepoint each tissue x ome has.
#
# Embedding (18 dims), order = tissue (adipose -> blood -> muscle), then ome (rna -> prot),
# then timepoint (0.5 -> 4 -> 24 h). Value = logFC of the delta-delta contrast
# (arm change from baseline minus control change from baseline), scaled per ome block.
# One file per arm: EE-CON and RE-CON.
#
# Scaling: omes have very different logFC magnitudes (RNA ~0.1, OLINK ~0.5+), so each
# tissue x ome block is divided by its root-mean-square logFC, computed over the 471 genes,
# all timepoints the block has and BOTH arms. One factor per block means timepoint and arm
# differences inside a block are preserved, zero stays "no change", and signs are untouched.
# Factors are written to 01_scale_factors.csv; raw logFC to 01_nodes_<arm>_raw_logFC.csv.
#
# Structural NAs: adipose protein exists only at 4 h (0.5 h and 24 h are NA by design).
#
# Many-to-one: when several features map to one gene within a tissue x ome, the feature with
# the highest AveExpr is kept (arm- and contrast-independent, so no selection on the response).
# The chosen feature ids are written to a separate provenance file.
#
# Output root: $HACK_OUT (default ~/Desktop/output/hackathon-2026-track1/network).

suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TISSUES <- c("adipose", "blood", "muscle")
OMES    <- c(rna = "TRNSCRPT", prot = NA)                       # prot object is tissue-specific
PROT_OBJ <- c(adipose = "ADIPOSE_PROT_PR_DA", blood = "BLOOD_PROT_OL_DA", muscle = "MUSCLE_PROT_PR_DA")
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
ARMS <- c("EE-CON", "RE-CON")

map <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[
  assay %in% c("transcript-rna-seq", "prot-pr", "prot-ol"),
  .(assay = as.character(assay), feature_id = as.character(feature_id),
    entrez_gene = as.character(entrez_gene), gene_symbol = as.character(gene_symbol))
])[!is.na(entrez_gene)]

# ---- load every tissue x ome, delta-delta contrasts at the 3 post timepoints -------------
cells <- CJ(tissue = TISSUES, ome = names(OMES), sorted = FALSE)
da <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  tis <- cells$tissue[i]; ome <- cells$ome[i]
  obj <- if (ome == "rna") paste0(toupper(tis), "_TRNSCRPT_DA") else PROT_OBJ[[tis]]
  x <- as.data.table(get(obj))[contrast_category %in% ARMS & Timepoint %in% TPS,
         .(assay, feature_id, arm = as.character(contrast_category),
           tp = names(TPS)[match(as.character(Timepoint), TPS)], logFC, AveExpr)]
  x[, `:=`(tissue = tis, ome = ome)]
  merge(x, map, by = c("assay", "feature_id"))
}))

# ---- universe: genes present in all 6 tissue x ome cells --------------------------------
genes_by_cell <- da[, .(genes = list(unique(entrez_gene))), by = .(tissue, ome)]
universe <- Reduce(intersect, genes_by_cell$genes)
stopifnot(length(universe) == 471)
da <- da[entrez_gene %in% universe]

# ---- one feature per gene per tissue x ome: highest AveExpr -----------------------------
feat_expr <- unique(da[, .(tissue, ome, entrez_gene, feature_id, AveExpr)])[
  , .(AveExpr = mean(AveExpr)), by = .(tissue, ome, entrez_gene, feature_id)]
setorder(feat_expr, tissue, ome, entrez_gene, -AveExpr, feature_id)
feat_expr[, n_candidates := .N, by = .(tissue, ome, entrez_gene)]
chosen <- feat_expr[, .SD[1], by = .(tissue, ome, entrez_gene)]
da <- da[chosen[, .(tissue, ome, entrez_gene, feature_id)], on = .(tissue, ome, entrez_gene, feature_id)]
stopifnot(!anyDuplicated(da[, .(tissue, ome, entrez_gene, arm, tp)]))

# ---- per-block RMS scaling --------------------------------------------------------------
scale_f <- da[, .(rms_logFC = sqrt(mean(logFC^2)), n_values = .N,
                  timepoints = paste(names(TPS)[names(TPS) %in% tp], collapse = ",")),
              by = .(tissue, ome)]
da <- scale_f[, .(tissue, ome, rms_logFC)][da, on = .(tissue, ome)]
da[, scaled := logFC / rms_logFC]
fwrite(scale_f, file.path(OUT, "01_scale_factors.csv"))
print(scale_f)

# ---- 18-dim embedding, one file per arm ------------------------------------------------
dims <- CJ(tissue = TISSUES, ome = names(OMES), tp = names(TPS), sorted = FALSE)[
  , dim := paste(tissue, ome, tp, sep = "_")]
symbols <- unique(map[entrez_gene %in% universe, .(entrez_gene, gene_symbol)])[
  , .(gene_symbol = gene_symbol[1]), by = entrez_gene]

to_wide <- function(a, value) {
  long <- da[arm == a][dims, on = .(tissue, ome, tp)]
  wide <- dcast(long[!is.na(entrez_gene)], entrez_gene ~ dim, value.var = value)
  # re-add structurally absent dims (adipose prot 0.5h/24h) as NA columns, then fix order
  for (d in setdiff(dims$dim, names(wide))) set(wide, j = d, value = NA_real_)
  wide <- symbols[wide, on = "entrez_gene"]
  setcolorder(wide, c("entrez_gene", "gene_symbol", dims$dim))
  setorder(wide, gene_symbol)
  stopifnot(nrow(wide) == 471, ncol(wide) == 2 + 18)
  wide
}

for (a in ARMS) {
  tag <- sub("-CON", "", a)
  f <- file.path(OUT, sprintf("01_nodes_%s.csv", tag))
  fwrite(to_wide(a, "scaled"), f)
  fwrite(to_wide(a, "logFC"), file.path(OUT, sprintf("01_nodes_%s_raw_logFC.csv", tag)))
  message(sprintf("%s: 471 nodes x 18 dims -> %s", a, f))
}

prov <- dcast(chosen, entrez_gene ~ paste(tissue, ome, sep = "_"),
              value.var = c("feature_id", "n_candidates"))
prov <- symbols[prov, on = "entrez_gene"]
fwrite(prov, file.path(OUT, "01_nodes_feature_provenance.csv"))
message("provenance -> ", file.path(OUT, "01_nodes_feature_provenance.csv"))
