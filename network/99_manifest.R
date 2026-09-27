#!/usr/bin/env Rscript
# =====================================================================================================
# 99_manifest.R — FINGERPRINT EVERY OUTPUT, AND CHECK A RERUN REPRODUCES THEM EXACTLY
# =====================================================================================================
# PURPOSE: the reproducibility check behind "run the path and compare expected outputs". Lists every output
#   the pipeline wrote (tables in $HACK_OUT, figures / pages / PDF in $HACK_FIG) with its MD5 checksum.
#   The first run saves a REFERENCE; every later run is compared with it file by file.
# HOW TO RUN: Rscript network/99_manifest.R            (compare with the reference, or create it if absent)
#             MANIFEST_RESET=1 Rscript network/99_manifest.R   (make the current outputs the new reference)
#   network/run_all.sh calls it last.
# OUTPUTS: $HACK_OUT/99_manifest.csv (current), 99_manifest_reference.csv, 99_manifest_diff.csv (differences)
# EXCLUDED (not pipeline results, or inherently variable): logs, 00_environment, the GlyGen API cache and
#   downloaded files, archived outputs, the neo4j export (owned by teammates), LaTeX auxiliary files.
# =====================================================================================================
suppressMessages(library(data.table))
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# every file under the two output folders, minus the excluded ones
f_out <- list.files(OUT, recursive = TRUE, full.names = TRUE)
f_fig <- list.files(FIG, recursive = TRUE, full.names = TRUE, pattern = "\\.(png|pdf|html|tex)$")
skip <- "(/logs/|/00_environment/|/glygen_cache/|/glygen_files/|/archive|/neo4j_import/|99_manifest|\\.(log|aux)$|/\\.DS_Store$)"
files <- c(f_out, f_fig)
files <- files[!grepl(skip, files)]
# resources (the team's mnet folder) live under FIG but are inputs, not outputs
files <- files[!grepl("/resources/", files)]
M <- data.table(where = fifelse(startsWith(files, OUT), "HACK_OUT", "HACK_FIG"),
                file = sub(paste0("^(", OUT, "|", FIG, ")/"), "", files), bytes = file.size(files), md5 = unname(tools::md5sum(files)))[order(where, file)]
fwrite(M, file.path(OUT, "99_manifest.csv"))
ref <- file.path(OUT, "99_manifest_reference.csv")
if (!file.exists(ref) || Sys.getenv("MANIFEST_RESET") == "1") {
  fwrite(M, ref); message(sprintf("reference saved: %d output files fingerprinted (%s)", nrow(M), ref))
} else {
  R <- fread(ref)
  D <- merge(R[, .(where, file, md5_reference = md5)], M[, .(where, file, md5_now = md5)], by = c("where", "file"), all = TRUE)
  D[, status := fcase(is.na(md5_reference), "new", is.na(md5_now), "missing now", md5_reference == md5_now, "identical", default = "DIFFERENT")]
  fwrite(D[status != "identical"], file.path(OUT, "99_manifest_diff.csv"))
  print(D[, .N, by = status])
  if (any(D$status != "identical")) print(D[status != "identical", .(where, file, status)], nrows = 100)
  message(if (all(D$status == "identical")) sprintf("REPRODUCED: all %d output files are byte-identical to the reference", nrow(D))
          else sprintf("%d of %d files differ from the reference (see 99_manifest_diff.csv)", sum(D$status != "identical"), nrow(D)))
}
