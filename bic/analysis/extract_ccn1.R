# Read published human MoTrPAC model results; no participant-level model fitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)
input <- normalizePath(args[1], mustWork=TRUE)
output <- args[2]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
read_object <- function(name) {
  env <- new.env(parent=emptyenv())
  load(file.path(input, paste0(name, ".rda")), env)
  as.data.frame(env[[name]])
}
mapping <- read_object("HUMAN_FEATURE_TO_GENE")
map <- unique(mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol %in% c("CCN1", "CYR61"),
                      c("assay", "feature_id", "gene_symbol", "entrez_gene", "ensembl_gene", "uniprot")])
objects <- c("BLOOD_PROT_OL_DA", "BLOOD_TRNSCRPT_DA", "MUSCLE_TRNSCRPT_DA", "ADIPOSE_TRNSCRPT_DA",
             "MUSCLE_PROT_PR_DA", "MUSCLE_PROT_PH_DA", "ADIPOSE_PROT_PR_DA", "ADIPOSE_PROT_PH_DA")
results <- list(); coverage <- list(); backgrounds <- list(); analyzed_assays <- character()
columns <- c("source_object", "tissue", "assay", "platform", "feature_id", "gene_symbol",
             "contrast", "contrast_short", "contrast_type", "contrast_category", "Timepoint",
             "full_model", "logFC", "CI.L_calculated", "CI.R_calculated", "p_value", "adj_p_value", "t", "z.std")
for (object in objects) {
  dat <- read_object(object)
  if (!"platform" %in% names(dat)) dat$platform <- as.character(dat$assay)
  analyzed_assays <- union(analyzed_assays, as.character(unique(dat$assay)))
  lookup <- unique(map[map$assay %in% unique(dat$assay), c("assay", "feature_id", "gene_symbol")])
  stopifnot(!anyDuplicated(lookup[, c("assay", "feature_id")]))
  hit <- merge(dat, lookup, by=c("assay", "feature_id"))
  hit$source_object <- rep(object, nrow(hit))
  results[[object]] <- hit[, columns]
  coverage[[object]] <- data.frame(source_object=object, tissue=as.character(unique(dat$tissue)),
    assay=as.character(unique(dat$assay)), total_analyzed_features=length(unique(dat$feature_id)),
    mapped_CCN1_features=length(unique(lookup$feature_id)),
    analyzed_CCN1_features=length(unique(hit$feature_id)), CCN1_contrast_rows=nrow(hit))
  # Full feature universe, retained before filtering to CCN1, for every focal assay.
  if (nrow(hit) > 0) {
    bg <- dat[dat$contrast_type %in% c("exercise_with_controls", "Endur_vs_Resist"), ]
    bg$source_object <- object
    backgrounds[[object]] <- bg[, c("source_object", "tissue", "assay", "platform", "feature_id",
      "contrast", "contrast_type", "contrast_category", "Timepoint", "logFC", "p_value", "adj_p_value")]
  }
  cat(object, ":", nrow(hit), "CCN1 contrast rows\n")
}
combined <- do.call(rbind, results[vapply(results, nrow, integer(1)) > 0])
rownames(combined) <- NULL
stopifnot(nrow(combined) > 0)
write.csv(combined, file.path(output, "ccn1_source_contrasts.csv"), row.names=FALSE)
write.csv(do.call(rbind, coverage), file.path(output, "ccn1_assay_coverage.csv"), row.names=FALSE)
write.csv(map[map$assay %in% analyzed_assays, ], file.path(output, "ccn1_feature_mapping.csv"), row.names=FALSE)
write.csv(map[!map$assay %in% analyzed_assays, ], file.path(output, "ccn1_other_mapping_not_analyzed.csv"), row.names=FALSE)
con <- gzfile(file.path(output, "assay_background.csv.gz"), "wt")
write.csv(do.call(rbind, backgrounds), con, row.names=FALSE)
close(con)
cat("Saved", nrow(combined), "source contrasts and full backgrounds for", length(backgrounds), "assays.\n")
