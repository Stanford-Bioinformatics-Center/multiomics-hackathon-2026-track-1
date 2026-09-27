#!/usr/bin/env Rscript
# =====================================================================================================
# 00_environment.R — STEP 0: RECORD THE ENVIRONMENT AND THE EXACT EXTERNAL INPUTS (REPRODUCIBILITY)
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   "Can a future learner set up the environment, obtain the data, run the path and compare expected
#   outputs?" (hackathon documentation standard). This script records what the results depend on: the
#   software versions, and a fingerprint (MD5 checksum, size, date) of every external input file, so a rerun
#   elsewhere can be checked against the same inputs.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Writes the R session (R version, platform, every attached / loaded package with its version), the
#      versions of Python, pandoc and TinyTeX, and the versions of the key R packages the pipeline uses.
#   2. Fingerprints every external input: the MoTrPAC packages (version), the team's mnet resource files, the
#      legacy STRING file, the cached Rhea files (+ release), the cached GlyGen responses (+ release, read live
#      from GlyGen's dataset service), and the Downloads parquet if present.
#   3. Writes network/ENVIRONMENT.md (committed: versions and checksums only, no data) and the same tables to
#      $HACK_OUT/00_environment/.
#
# HOW TO RUN:  Rscript network/00_environment.R   (seconds; run by network/run_all.sh first)
# TECH STACK:  R 4.4; data.table, tools (md5sum).
# OUTPUTS:     network/ENVIRONMENT.md; $HACK_OUT/00_environment/{session_info.txt, packages.csv, inputs.csv}
# =====================================================================================================

# Load packages quietly.
suppressMessages(library(data.table))

# Folders (outside the repo) and the repo's network folder (this script's folder).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
MNET_DIR <- Sys.getenv("MNET_DIR", unset = path.expand("~/Desktop/output/hackathon/resources/mo_annotation"))
EXT <- Sys.getenv("HACK_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/rhea"))
STRING_PARQUET <- Sys.getenv("STRING_PARQUET", unset = path.expand("~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet"))
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) "network")
ENV_DIR <- file.path(OUT, "00_environment"); dir.create(ENV_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. software versions ---------------------------------------------------------------------------------
# Helper: first line of a command's version output ("" if the tool is missing).
tool_version <- function(cmd, args) tryCatch(system2(cmd, args, stdout = TRUE, stderr = TRUE)[1], error = function(e) "not found", warning = function(w) "not found")
# pdflatex from TinyTeX (installed under the user's library, often not on the shell PATH)
pdflatex_bin <- function() { p <- Sys.which("pdflatex"); if (nzchar(p)) return(p)
  r <- tryCatch(tinytex::tinytex_root(), error = function(e) ""); c(list.files(file.path(r, "bin"), pattern = "^pdflatex$", recursive = TRUE, full.names = TRUE), "pdflatex")[1] }
KEY_PKGS <- c("data.table", "igraph", "ggplot2", "ggrepel", "ggforce", "ggnewscale", "scales", "nanoparquet", "jsonlite", "visNetwork",
              "htmlwidgets", "htmltools", "rmarkdown", "tinytex", "limma", "readxl", "checkmate", "testthat", "pkgload", "patchwork", "exnet", "TMSig", "fgsea", "MotrpacHumanPreSuspensionAnalysis",
              "MotrpacHumanPreSuspensionData", "MotrpacHumanPreSuspension", "Matrix", "repfdr", "curl")
pk <- data.table(package = KEY_PKGS, version = sapply(KEY_PKGS, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")))
tools <- data.table(tool = c("R", "Python", "pandoc", "TinyTeX (pdflatex)", "platform"),
                    version = c(R.version.string, tool_version("python3", "--version"),
                                tryCatch(as.character(rmarkdown::pandoc_version()), error = function(e) "not found"),
                                tool_version(pdflatex_bin(), "--version"), R.version$platform))
fwrite(pk, file.path(ENV_DIR, "packages.csv"))
writeLines(capture.output(sessionInfo()), file.path(ENV_DIR, "session_info.txt"))

# ---- 2. external inputs: fingerprints ----------------------------------------------------------------------
fp <- function(files, source) {
  files <- files[file.exists(files)]
  if (!length(files)) return(data.table(source = source, file = NA_character_, bytes = NA_real_, modified = NA_character_, md5 = "MISSING"))
  data.table(source = source, file = basename(files), bytes = file.size(files),
             modified = format(file.mtime(files), "%Y-%m-%d %H:%M"), md5 = unname(tools::md5sum(files)))
}
# GlyGen: the release the service reports now (cached per-protein responses are what the pipeline uses)
glygen_release <- tryCatch(jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory("https://dsapi.glygen.org/dataset/init",
  handle = curl::handle_setopt(curl::new_handle(), postfields = "{}", httpheader = "Content-Type: application/json"))$content))$record$dataversion,
  error = function(e) "unavailable (offline)")
gg_cache <- list.files(file.path(OUT, "inventory", "glygen_cache"), full.names = TRUE)
rhea_rel <- tryCatch(paste(sub(".*=", "", grep("release", readLines(file.path(EXT, "rhea-release.properties")), value = TRUE)), collapse = " / "), error = function(e) "not cached")
INP <- rbind(
  fp(list.files(MNET_DIR, full.names = TRUE, pattern = "\\.(csv|parquet|md)$"), "team mnet resource (MNET_DIR)"),
  fp(STRING_PARQUET, "legacy curated STRING file (STRING_PARQUET; used only with EDGE_SOURCE=legacy)"),
  fp(Sys.getenv("DISEASE_SCORES", unset = path.expand("~/Desktop/output/week_6/_shared/disease_scores.csv.gz")),
     "Amar et al. 2024 disease proteomics sets, processed in Venus week 6 (DISEASE_SCORES; step 18d)"),
  fp(list.files(Sys.getenv("DISEASE_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/disease")), recursive = TRUE, full.names = TRUE,
                pattern = "\\.(xlsx|txt|R|md)$|COMMIT_SHA"), "published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f)"),
  fp(Sys.getenv("UBAIDA", unset = path.expand("~/Desktop/github/motrpac/motrpac-rat-training-mitochondria/disease_datasets/ubaida_mohien_ 2019_elife_stat.csv")),
     "Ubaida-Mohien 2019 muscle ageing proteome, Amar et al. repository copy (UBAIDA; step 20)"),
  fp(list.files(EXT, full.names = TRUE), paste0("Rhea cache (HACK_EXT), release ", rhea_rel)),
  data.table(source = paste0("GlyGen API cache (inventory/glygen_cache), service release now: ", glygen_release), file = sprintf("%d per-protein files", length(gg_cache)),
             bytes = sum(file.size(gg_cache)), modified = if (length(gg_cache)) format(max(file.mtime(gg_cache)), "%Y-%m-%d %H:%M") else NA_character_,
             md5 = if (length(gg_cache)) "see below" else "MISSING"),
  data.table(source = "MoTrPAC packages (installed)", file = c("MotrpacHumanPreSuspensionAnalysis", "MotrpacHumanPreSuspensionData"),
             bytes = NA_real_, modified = NA_character_, md5 = paste0("version ", pk[package %in% c("MotrpacHumanPreSuspensionAnalysis", "MotrpacHumanPreSuspensionData"), version])),
  fill = TRUE)
# a single fingerprint of the GlyGen cache: the MD5 of the sorted list of per-file MD5s
if (length(gg_cache)) { tf <- tempfile(); writeLines(unname(tools::md5sum(sort(gg_cache))), tf); INP[grepl("^GlyGen", source), md5 := unname(tools::md5sum(tf))] }
fwrite(INP, file.path(ENV_DIR, "inputs.csv"))

# ---- 3. ENVIRONMENT.md (committed; no data) ------------------------------------------------------------------
md <- c("# Environment and input fingerprints", "",
        sprintf("Recorded by `network/00_environment.R` on %s. Re-run it (or `network/run_all.sh`) and compare: the same", format(Sys.Date())),
        "versions and input checksums should give identical outputs (checked with `network/99_manifest.R`).", "",
        "## Software", "", "| Tool | Version |", "|---|---|", sprintf("| %s | %s |", tools$tool, gsub("\\|", "/", tools$version)), "",
        "## R packages used by the pipeline", "", "| Package | Version |", "|---|---|", sprintf("| %s | %s |", pk$package, pk$version), "",
        "## External inputs (MD5 fingerprints)", "", "| Source | File | Bytes | Modified | MD5 |", "|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s |", INP$source, INP$file, format(INP$bytes, big.mark = ",", scientific = FALSE, trim = TRUE), INP$modified, INP$md5), "",
        "Full `sessionInfo()` is written to `$HACK_OUT/00_environment/session_info.txt` (not committed).")
writeLines(md, file.path(HERE, "ENVIRONMENT.md"))
message("-> ", file.path(HERE, "ENVIRONMENT.md"), " and ", ENV_DIR)
print(tools); print(INP[, .(source = substr(source, 1, 60), file, md5 = substr(md5, 1, 12))])
