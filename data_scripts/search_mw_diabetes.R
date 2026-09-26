#!/usr/bin/env Rscript
# =====================================================================================================
# search_mw_diabetes.R — RANK DIABETES METABOLOMICS STUDIES BY HOW MANY OF OUR METABOLITES THEY MEASURED
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Metabolomics Workbench is the Common Fund Metabolomics program's repository, one of the programs in the
#   CFDE (cfde.cloud). Its public REST API (no login) lets us:
#     1. find studies whose title contains a diabetes-related word (KEYWORDS below), keeping only human
#        studies (SPECIES below; override with HACK_SPECIES);
#     2. download the list of named metabolites each study measured, with their RefMet names;
#     3. match those to our metabolites (metabolites_all.csv from export_all_features.R, also RefMet names)
#        and count how many of the network's 450, and of all 1,460, each study measured.
#   Studies measuring many of the network's metabolites are the best candidates for checking the network in
#   diabetes. Only names and counts are written; no measurements are downloaded.
#
# MATCHING
#   Case-insensitive exact match of names. The study's RefMet name is used when it has one, otherwise the name
#   the submitter gave. Untargeted studies that report no named metabolites score 0 (n_named == 0), which does
#   not mean the study is useless, only that it cannot be scored this way.
#
# CAVEATS
#   Title search is by substring, so "glucose" also finds glucose-tracer and cell-culture studies; the
#   keywords column says which word matched, and species / sample counts help sort them out.
#
# TECH STACK
#   R 4.4; data.table; jsonlite (reads the API's JSON).
#
# INPUTS:  $HACK_RES/metabolites_all.csv (run export_all_features.R first)
# OUTPUTS: $HACK_DIAB/mw_diabetes_studies_ranked.csv   one row per study, best first
#          $HACK_DIAB/mw_diabetes_study_metabolites.csv  one row per study x matched metabolite
#          $HACK_DIAB/mw_cache/                         raw API answers, so reruns do not refetch
#          ($HACK_DIAB defaults to ~/Desktop/output/hackathon-2026-track1/diabetes_data; never inside the repo)
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(jsonlite) })

# Results folder of the network pipeline (override with HACK_OUT); used only to find the resource folder.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where export_all_features.R wrote its tables (override with HACK_RES).
RES <- Sys.getenv("HACK_RES", unset = file.path(OUT, "resource"))
# Where this script writes (override with HACK_DIAB); outside the repo so no data is committed.
DIAB <- Sys.getenv("HACK_DIAB", unset = path.expand("~/Desktop/output/hackathon-2026-track1/diabetes_data"))
# Cache of raw API answers (delete it to refetch).
CACHE <- file.path(DIAB, "mw_cache")
# Create both if needed.
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

# Metabolomics Workbench REST API.
API <- "https://www.metabolomicsworkbench.org/rest"
# Title words searched (substring match, so "diabet" also finds "diabetic").
KEYWORDS <- c("diabet", "insulin", "glucose", "glycemi", "T2D", "T1D", "islet")
# Species kept (override with HACK_SPECIES, ";"-separated; empty keeps all species).
SPECIES <- strsplit(Sys.getenv("HACK_SPECIES", unset = "Homo sapiens"), ";", fixed = TRUE)[[1]]
# Allow slow API answers (large studies can take over a minute).
options(timeout = 300)
# Requests that failed after retries.
FAILED <- character(0)

# Helper: fetch one API URL as parsed JSON, from the cache if already fetched.
get_json <- function(path) {
  # Cache file named after the request path.
  f <- file.path(CACHE, paste0(gsub("[^A-Za-z0-9]+", "_", path), ".json"))
  # Fetch only if not cached: up to 3 tries, into a temporary file renamed only on success so a failed
  # download is never cached; a short pause between requests to be polite to the server.
  if (!file.exists(f)) {
    tmp <- paste0(f, ".part")
    for (i in 1:3) {
      ok <- tryCatch(download.file(paste0(API, "/", path), tmp, quiet = TRUE) == 0, error = function(e) FALSE, warning = function(w) FALSE)
      Sys.sleep(0.2)
      if (ok) break
    }
    # Still failing: warn and treat as no answer (the study is reported with fetch_failed = TRUE).
    if (!ok) { unlink(tmp); warning("could not fetch ", path); FAILED <<- c(FAILED, path); return(list()) }
    file.rename(tmp, f)
  }
  # Read the text; the API answers "[]" or an empty body when nothing matches.
  txt <- paste(readLines(f, warn = FALSE), collapse = "")
  if (!nzchar(txt) || txt == "[]") return(list())
  # One match comes back as a single record, several as {"1": record, "2": record, ...}: always return a list of records.
  x <- fromJSON(txt, simplifyVector = FALSE)
  if (!is.null(x$study_id)) list(x) else unname(x)
}

# Helper: list of records -> data.table (NULL fields become NA).
as_dt <- function(recs) rbindlist(lapply(recs, function(r) lapply(r, function(v) if (is.null(v)) NA_character_ else as.character(v))), fill = TRUE)

# Helper: name normalised for matching.
norm <- function(x) tolower(trimws(x))

# ---- our metabolites ---------------------------------------------------------------------------------
# The table from export_all_features.R.
f_met <- file.path(RES, "metabolites_all.csv")
if (!file.exists(f_met)) stop("missing ", f_met, "; run data_scripts/export_all_features.R first")
ours <- fread(f_met)[, .(metabolite, in_network, key = norm(metabolite))]
# Safety check: same counts as the export script, and no two names collide after normalising.
stopifnot(nrow(ours) == 1460, sum(ours$in_network) == 450, !anyDuplicated(ours$key))

# ---- candidate studies -------------------------------------------------------------------------------
# One search per keyword, keeping which keyword found each study.
studies <- rbindlist(lapply(KEYWORDS, function(k) {
  x <- as_dt(get_json(paste0("study/study_title/", k, "/summary")))
  if (nrow(x)) x[, keyword := k]
  x
}), fill = TRUE)
# One row per study, with all the keywords that found it.
studies <- studies[, .(keywords = paste(keyword, collapse = ";")),
                   by = .(study_id, study_title, species, analysis_type, number_of_samples, institute, release_date)]
# Safety check: one row per study.
stopifnot(!anyDuplicated(studies$study_id))
message(sprintf("%d studies matched the title keywords", nrow(studies)))
# Keep only the chosen species (human by default).
if (length(SPECIES)) studies <- studies[species %in% SPECIES]
message(sprintf("%d of them are %s", nrow(studies), paste(SPECIES, collapse = " or ")))

# ---- metabolites measured by each study ----------------------------------------------------------------
# Named metabolites of every study (one row per analysis x metabolite).
mw <- rbindlist(lapply(studies$study_id, function(s) as_dt(get_json(paste0("study/study_id/", s, "/metabolites")))), fill = TRUE)
# Name used for matching: the RefMet name if given, otherwise the submitted name.
mw[, key := norm(fifelse(is.na(refmet_name) | refmet_name == "", metabolite_name, refmet_name))]
# Distinct metabolites per study (a metabolite measured in several analyses of one study counts once).
mw <- unique(mw[!is.na(key) & key != "", .(study_id, key)])
# Number of distinct named metabolites per study.
n_named <- mw[, .(n_named = .N), by = study_id]
# The ones that are ours.
hits <- ours[mw, on = "key", nomatch = NULL][, .(study_id, metabolite, in_network)]

# ---- score and rank ----------------------------------------------------------------------------------
# Counts per study: of our 1,460, and of the network's 450.
score <- hits[, .(n_matched_all = .N, n_matched_network = sum(in_network)), by = study_id]
# Attach to the studies (studies with no named metabolites or no matches get 0).
ranked <- score[n_named[studies, on = "study_id"], on = "study_id"]
for (col in c("n_named", "n_matched_all", "n_matched_network")) set(ranked, which(is.na(ranked[[col]])), col, 0L)
# Studies whose metabolite list could not be fetched (their 0s are unknown, not true zeros; rerun to retry).
ranked[, fetch_failed := paste0("study/study_id/", study_id, "/metabolites") %in% FAILED]
# Share of the network's 450 the study measured, and a link to the study page.
ranked[, pct_network := round(100 * n_matched_network / 450, 1)]
ranked[, url := paste0("https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=", study_id)]
# Best first: most network metabolites, then most metabolites overall.
setorder(ranked, -n_matched_network, -n_matched_all, study_id)
# Column order.
setcolorder(ranked, c("study_id", "study_title", "species", "analysis_type", "number_of_samples", "keywords",
                      "n_named", "n_matched_all", "n_matched_network", "pct_network"))
# Save.
fwrite(ranked, file.path(DIAB, "mw_diabetes_studies_ranked.csv"))

# The matched metabolites of each study, network ones first, in the ranked study order.
hits[, study_id := factor(study_id, levels = ranked$study_id)]
setorder(hits, study_id, -in_network, metabolite)
fwrite(hits, file.path(DIAB, "mw_diabetes_study_metabolites.csv"))

# ---- report ------------------------------------------------------------------------------------------
message(sprintf("%d studies report named metabolites; %d match at least one network metabolite; %d could not be fetched",
                sum(ranked$n_named > 0), sum(ranked$n_matched_network > 0), sum(ranked$fetch_failed)))
# Top 15 in the log.
print(ranked[1:min(15, .N), .(study_id, species, analysis_type, n = number_of_samples, n_named,
                              n_matched_network, title = substr(study_title, 1, 60))])
message("written to ", DIAB)
