# Audit the extracted candidate rows against the pinned source R objects.
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
root <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
all <- read.csv(file.path(root, "table3/all_candidate_results.csv"), stringsAsFactors=FALSE)
plasma <- read.csv(file.path(root, "table3/plasma_results.csv"), stringsAsFactors=FALSE)
ids <- unique(plasma$candidate_id[plasma$identity_resolved])
d <- all[all$candidate_id %in% ids, ]
keycols <- c("tissue", "assay", "platform", "feature_id", "contrast_type", "contrast_category", "Timepoint")
numcols <- c("logFC", "CI.L_calculated", "CI.R_calculated", "p_value", "adj_p_value")
key <- function(x) do.call(paste, c(x[keycols], sep="\t"))
checked <- 0L
for (obj in unique(d$source_object)) {
  e <- new.env()
  load(file.path(root, paste0(obj, ".rda")), e)
  source <- as.data.frame(e[[obj]])
  if (!"platform" %in% names(source)) source$platform <- as.character(source$assay)
  source_key <- key(source)
  stopifnot(!anyDuplicated(source_key))
  extracted <- d[d$source_object == obj, ]
  m <- match(key(extracted), source_key)
  stopifnot(!anyNA(m))
  for (column in numcols) {
    stopifnot(isTRUE(all.equal(extracted[[column]], source[[column]][m], tolerance=1e-12, check.attributes=FALSE)))
  }
  checked <- checked + nrow(extracted)
  cat(obj, nrow(extracted), "rows verified\n")
}
cat("VERIFIED", checked, "rows; all effects, CIs, raw p-values, and source adjusted p-values match pinned R objects.\n")
