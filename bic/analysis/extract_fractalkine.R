# Extract the human CX3CL1 results and full plasma Olink comparison universe.
# Base R only; no network access, participant-level data, or model fitting.
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
gene <- "CX3CL1"
mapping <- read_object("HUMAN_FEATURE_TO_GENE")
map <- unique(mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol == gene,
                      c("assay", "feature_id", "gene_symbol")])
objects <- c("BLOOD_PROT_OL_DA", "BLOOD_TRNSCRPT_DA",
             "MUSCLE_TRNSCRPT_DA", "ADIPOSE_TRNSCRPT_DA",
             "MUSCLE_PROT_PR_DA", "MUSCLE_PROT_PH_DA",
             "ADIPOSE_PROT_PR_DA", "ADIPOSE_PROT_PH_DA")
results <- list()
coverage <- list()
for (object in objects) {
  dat <- read_object(object)
  if (!"platform" %in% names(dat)) dat$platform <- as.character(dat$assay)
  lookup <- map[map$assay %in% unique(dat$assay), ]
  hit <- merge(dat, lookup, by=c("assay", "feature_id"))
  hit$source_object <- rep(object, nrow(hit))
  columns <- c("source_object", "tissue", "assay", "platform", "feature_id", "gene_symbol",
               "contrast", "contrast_short", "contrast_type", "contrast_category", "Timepoint",
               "full_model", "logFC", "CI.L_calculated", "CI.R_calculated", "p_value",
               "adj_p_value", "t", "z.std")
  stopifnot(all(columns %in% names(hit)))
  results[[object]] <- hit[, columns]
  coverage[[object]] <- data.frame(
    source_object=object, tissue=as.character(unique(dat$tissue)),
    assay=as.character(unique(dat$assay)),
    total_analyzed_features=length(unique(dat$feature_id)),
    mapped_CX3CL1_features=length(unique(lookup$feature_id)),
    analyzed_CX3CL1_features=length(unique(hit$feature_id)),
    CX3CL1_contrast_rows=nrow(hit), stringsAsFactors=FALSE)
  if (object == "BLOOD_PROT_OL_DA") {
    background <- dat[dat$contrast_type %in% c("exercise_with_controls", "Endur_vs_Resist"), ]
    background$source_object <- object
    mp <- unique(mapping[mapping$assay == "prot-ol" & !is.na(mapping$gene_symbol), c("feature_id", "gene_symbol")])
    symbols <- aggregate(gene_symbol ~ feature_id, mp, function(x) paste(sort(unique(x)), collapse=";"))
    background$gene_symbol <- symbols$gene_symbol[match(background$feature_id, symbols$feature_id)]
    write.csv(background, file.path(output, "plasma_olink_background.csv"), row.names=FALSE)
  }
  cat(object, ":", nrow(hit), "CX3CL1 contrast rows\n")
}
combined <- do.call(rbind, results[vapply(results, nrow, integer(1)) > 0])
rownames(combined) <- NULL
stopifnot(nrow(combined) > 0, all(combined$gene_symbol == gene))
write.csv(combined, file.path(output, "cx3cl1_source_contrasts.csv"), row.names=FALSE)
write.csv(do.call(rbind, coverage), file.path(output, "cx3cl1_assay_coverage.csv"), row.names=FALSE)
write.csv(map, file.path(output, "cx3cl1_feature_mapping.csv"), row.names=FALSE)
cat("Extracted", nrow(combined), "CX3CL1 rows from human c2.0 cached source objects.\n")
