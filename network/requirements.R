# =====================================================================================================
# network/requirements.R — INSTALL EVERY R PACKAGE THE MAIN PIPELINE (run_all.sh) NEEDS
# =====================================================================================================
# PURPOSE: one command from a fresh R (>= 4.4) to all packages used by the steps in network/run_all.sh,
#   then a check against the versions of the reference run (network/ENVIRONMENT.md).
# HOW TO RUN (from the repository root):
#   Rscript network/requirements.R            # current CRAN / Bioconductor versions (binaries, fast)
#   Rscript network/requirements.R --exact    # the reference run's exact CRAN versions (may build from source)
#   Rscript network/requirements.R --check    # install nothing; only report missing / different versions
# LIST MADE FROM: every library() / require() / pkg:: call in the scripts run_all.sh runs, the files they
#   source (R/figure_style.R, engine/), and the engine's DESCRIPTION. Not included: network/interaction_db
#   (its own requirements.txt) and network/video (Node; package.json) — run_all.sh does not use them.
# ALSO NEEDED (not R packages; see network/requirements.txt): Python >= 3.9 (standard library only), curl,
#   pandoc (steps 17i; ships with RStudio / Positron / Quarto), TinyTeX (step 13; installed below if missing).
# =====================================================================================================

args  <- commandArgs(trailingOnly = TRUE)
EXACT <- "--exact" %in% args
CHECK <- "--check" %in% args
if (getRversion() < "4.4") stop("R >= 4.4 is required (MotrpacHumanPreSuspensionAnalysis); this is R ", getRversion())
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Package, source, reference version (network/ENVIRONMENT.md, 2026-09-26; repfdr from step 18c's header), steps.
REQ <- read.csv(text = "
package,source,version,used_by
data.table,CRAN,1.18.2.1,all steps
igraph,CRAN,2.2.2,network components (02-21)
Matrix,CRAN,1.7.5,sparse matrices (igraph / limma)
nanoparquet,CRAN,0.4.3,reading .parquet edges (02 / 08)
ggplot2,CRAN,3.5.2,figures (10-21)
ggrepel,CRAN,0.9.8,figure labels (10-21)
ggforce,CRAN,0.5.0,class outlines (15-16)
ggnewscale,CRAN,0.5.2,PTM tag legends (16)
patchwork,CRAN,1.3.2,multi-panel figures (18-21 / M)
scales,CRAN,1.4.0,figure scales
visNetwork,CRAN,2.1.4,interactive pages (17i)
htmlwidgets,CRAN,1.6.4,interactive pages (17i)
htmltools,CRAN,0.5.9,interactive pages (17i)
jsonlite,CRAN,2.0.0,interactive pages / web lookups
curl,CRAN,7.0.0,GlyGen release lookup (00 / inventory)
rmarkdown,CRAN,2.30,pandoc detection (00 / 17i)
tinytex,CRAN,0.58,LaTeX table to PDF (13)
readxl,CRAN,1.4.5,published disease tables (20)
repfdr,CRAN,1.2.3,repfdr states (18c)
checkmate,CRAN,2.3.4,exnet engine
testthat,CRAN,3.3.2,exnet engine tests (14t)
pkgload,CRAN,1.5.0,loading the engine from source (14t / M / 21)
limma,Bioconductor,3.62.2,module tests (17s / 18d)
TMSig,Bioconductor,1.0.0,CAMERA-PR / ORA via MoTrPAC (17s / 18d / 18c)
fgsea,Bioconductor,1.32.4,dependency of the MoTrPAC package
MotrpacHumanPreSuspensionAnalysis,GitHub:MoTrPAC/MotrpacHumanPreSuspensionAnalysis,0.2.4,the data (all steps)
", strip.white = TRUE, stringsAsFactors = FALSE)
# Installed with the MoTrPAC package (checked, not installed separately).
DEPS_CHECKED <- c(MotrpacHumanPreSuspensionData = "0.0.1.4")

have <- function(p) requireNamespace(p, quietly = TRUE)
ver  <- function(p) if (have(p)) as.character(utils::packageVersion(p)) else NA_character_

if (!CHECK) {
  if (!have("pak")) install.packages("pak")
  if (!have("BiocManager")) install.packages("BiocManager")
  cran <- REQ[REQ$source == "CRAN", ]
  bioc <- REQ[REQ$source == "Bioconductor", ]
  gh   <- REQ[startsWith(REQ$source, "GitHub:"), ]
  # CRAN: exact reference versions with --exact, otherwise the current binaries of whatever is missing.
  stale <- if (EXACT) mapply(function(p, v) !identical(ver(p), v), cran$package, cran$version) else !vapply(cran$package, have, logical(1))
  todo  <- cran[stale, ]
  if (nrow(todo)) pak::pkg_install(if (EXACT) paste0(todo$package, "@", todo$version) else todo$package, ask = FALSE)
  # Bioconductor: the release that matches this R (3.20 for R 4.4, as in the reference run).
  todo <- bioc$package[!vapply(bioc$package, have, logical(1))]
  if (length(todo)) BiocManager::install(todo, update = FALSE, ask = FALSE)
  # GitHub: the MoTrPAC package (brings MotrpacHumanPreSuspensionData with it).
  todo <- gh[!vapply(gh$package, have, logical(1)), ]
  if (nrow(todo)) pak::pkg_install(sub("^GitHub:", "", todo$source), ask = FALSE)
  # TinyTeX for step 13, unless a LaTeX with pdflatex is already on the PATH.
  if (!nzchar(Sys.which("pdflatex")) && !tinytex::is_tinytex()) tinytex::install_tinytex()
}

# Report: every package against the reference run, then the non-R tools.
all_req <- rbind(REQ[, c("package", "version")],
                 data.frame(package = names(DEPS_CHECKED), version = unname(DEPS_CHECKED)))
all_req$installed <- vapply(all_req$package, ver, character(1))
all_req$status <- ifelse(is.na(all_req$installed), "MISSING",
                  ifelse(all_req$installed == all_req$version, "ok", "different version"))
cat("\nR", as.character(getRversion()), "(reference run: 4.4.3)\n")
print(all_req[, c("package", "version", "installed", "status")], row.names = FALSE, right = FALSE)
pandoc <- tryCatch(as.character(rmarkdown::pandoc_version()), error = function(e) NA)
cat("\npandoc:  ", if (length(pandoc) && !is.na(pandoc)) pandoc else "NOT FOUND (needed by step 17i; install RStudio / Quarto or pandoc)",
    "\npdflatex:", if (nzchar(Sys.which("pdflatex")) || (have("tinytex") && tinytex::is_tinytex())) "found" else "NOT FOUND (needed by step 13)",
    "\npython3: ", if (nzchar(Sys.which("python3"))) "found" else "NOT FOUND (needed by steps 01c / 01d / inv1)",
    "\ncurl:    ", if (nzchar(Sys.which("curl"))) "found" else "NOT FOUND (needed by step 20f)", "\n")
n_missing <- sum(all_req$status == "MISSING")
if (n_missing) { cat("\n", n_missing, "package(s) missing.\n"); quit(status = 1) }
cat("\nAll R packages present. Different versions can change outputs; 99_manifest.R will show which.\n")
