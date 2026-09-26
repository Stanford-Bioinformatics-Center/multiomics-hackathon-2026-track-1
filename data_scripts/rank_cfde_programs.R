#!/usr/bin/env Rscript
# =====================================================================================================
# rank_cfde_programs.R — WHICH CFDE PROGRAM'S DIABETES DATA MEASURED THE MOST OF OUR NETWORK, PER OME
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   search_mw_diabetes.R ranked individual Metabolomics Workbench studies. This script moves up one level, to
#   the Common Fund programs of the CFDE (cfde.cloud), and asks for each program's open, human, diabetes-
#   relevant data: how many of the network's 471 genes/proteins and 450 metabolites were measured, per ome
#   (transcriptomics, proteomics, metabolomics), and are the omes measured in the same people?
#   Only data that a CFDE program itself holds is counted. A survey of all 19 programs (their sites and papers,
#   September 2026) found open human data usable for diabetes in three programs only:
#     HMP     the iHMP type 2 diabetes / prediabetes study (HMP project "T2D"; Zhou et al. 2019, Nature): blood
#             RNA (PBMC), plasma proteome and metabolome in the same ~105 people, labelled insulin resistant /
#             sensitive. CFDE indexes it through HMP's C2M2 package, but only as raw files (FASTQ, mzML, RAW,
#             mzXML at drs://drs.hmpdacc.org). The gene/protein/metabolite tables used here are the study team's
#             processed version of those same samples, from the iPOP server at Stanford (checked September 2026:
#             887 of the 908 proteome samples in HMP's C2M2 package, 103 of 104 people, are in the protein table).
#             The iPOP EXERCISE sub-study on the same server is NOT in HMP's C2M2 package and is not used.
#     GTEx    RNA-seq of adipose, skeletal muscle and whole blood (the T2D donor flag is dbGaP-protected).
#     Metabolomics Workbench (NMDR)  the human studies ranked by search_mw_diabetes.R (best study, and the
#             union of all of them).
#   The other programs (no open human omics with diabetes labels) are listed in cfde_program_status.csv with
#   the reason.
#
# MATCHING
#   Genes/proteins: gene symbol (iHMP) or Ensembl gene ID without version (GTEx) against genes_all.csv.
#     Protein column names lose duplicate suffixes (".1", "_2").
#   Metabolites: each name is converted to its RefMet name with the Metabolomics Workbench RefMet API (our
#     names are already RefMet names), then matched case-insensitively. Lipids given with chains ("PC 16:0_18:1",
#     "TG 54:2_18:1") also match our sum-composition name ("PC 34:1", "TG 54:2").
#   Denominators: network = 471 genes / 450 metabolites. "all" = our genes measured as RNA (transcriptomics),
#     as protein (proteomics), or all 1,460 metabolites.
#
# TECH STACK
#   R 4.4; data.table; jsonlite; xml2 (reads the one .xlsx without extra packages).
#
# INPUTS:  $HACK_RES/genes_all.csv, metabolites_all.csv            (export_all_features.R)
#          $HACK_DIAB/mw_diabetes_studies_ranked.csv, mw_diabetes_study_metabolites.csv  (search_mw_diabetes.R)
#          iHMP (iPOP server) and GTEx files, downloaded once into $HACK_DIAB/cfde_cache/
# OUTPUTS: $HACK_DIAB/cfde_program_ranking.csv    one row per program x dataset x ome, with coverage counts
#          $HACK_DIAB/cfde_program_features.csv   one row per dataset x ome x matched feature of ours
#          $HACK_DIAB/cfde_program_status.csv     all 19 programs: what they have and why (not) counted
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(jsonlite); library(xml2) })

# Folders (same defaults and overrides as the other data_scripts).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
RES <- Sys.getenv("HACK_RES", unset = file.path(OUT, "resource"))
DIAB <- Sys.getenv("HACK_DIAB", unset = path.expand("~/Desktop/output/hackathon-2026-track1/diabetes_data"))
# Downloads and API answers (delete to refetch).
CACHE <- file.path(DIAB, "cfde_cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
# Allow slow servers.
options(timeout = 900)

# iPOP public file server and GTEx v8 gene median TPM table.
IPOP <- "http://ipop_public.gbsc.os.scg.stanford.edu"
GTEX <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
# GTEx tissues matching ours.
GTEX_TISSUES <- c("Adipose - Subcutaneous", "Muscle - Skeletal", "Whole Blood")

# Helper: download a URL once into the cache and return the local path.
fetch <- function(url) {
  f <- file.path(CACHE, basename(url))
  if (!file.exists(f)) {
    tmp <- paste0(f, ".part")
    download.file(url, tmp, quiet = TRUE, mode = "wb")
    file.rename(tmp, f)
  }
  f
}

# Helper: participant ID from an iPOP sample ID ("ZOZOW1T-1013" -> "ZOZOW1T").
subject_of <- function(x) sub("-.*$", "", x)

# ---- our features ------------------------------------------------------------------------------------
genes <- fread(file.path(RES, "genes_all.csv"), colClasses = list(character = c("entrez_gene", "uniprot", "ensembl_gene")))
mets <- fread(file.path(RES, "metabolites_all.csv"))
stopifnot(sum(genes$in_network) == 471, sum(mets$in_network) == 450)
# Which of our genes were measured as RNA, and as protein (the "all" denominators).
genes[, as_rna := rna_adipose | rna_blood | rna_muscle][, as_prot := prot_adipose | prot_blood | prot_muscle]
# Lookup keys: upper-case symbol, and each Ensembl ID (one row per ID).
sym_key <- genes[!is.na(gene_symbol) & gene_symbol != "", .(key = toupper(gene_symbol), entrez_gene)]
ens_key <- genes[!is.na(ensembl_gene) & ensembl_gene != "", .(key = unlist(strsplit(ensembl_gene, ";"))), by = entrez_gene]
# Metabolite keys: lower-case RefMet name; names like "PE P-39:6 or PE O-39:7" also match on either half.
met_key <- mets[, .(key = unique(tolower(c(metabolite, trimws(strsplit(metabolite, " or ", fixed = TRUE)[[1]]))))),
                by = .(metabolite, refmet_id, in_network)]

# Helper: our genes matched by a vector of symbols or Ensembl IDs (version suffix removed).
match_genes <- function(ids, by = c("symbol", "ensembl")) {
  by <- match.arg(by)
  ids <- unique(if (by == "symbol") toupper(ids) else sub("\\.[0-9]+$", "", ids))
  lk <- if (by == "symbol") sym_key else ens_key
  unique(lk[key %in% ids, .(entrez_gene)])[genes, on = "entrez_gene", nomatch = NULL][
    , .(feature = gene_symbol, id = entrez_gene, in_network, as_rna, as_prot)]
}

# ---- RefMet name lookup (cached) ---------------------------------------------------------------------
# All names looked up so far: name -> RefMet name ("" if RefMet does not know it).
REFMET_FILE <- file.path(CACHE, "refmet_lookup.csv")
refmet_map <- if (file.exists(REFMET_FILE)) fread(REFMET_FILE, colClasses = "character", na.strings = NULL) else
  data.table(name = character(), refmet_name = character())

# Helper: RefMet names of a vector of metabolite names (API called only for names not yet looked up).
refmet <- function(x) {
  todo <- setdiff(unique(x), refmet_map$name)
  if (length(todo)) message(sprintf("RefMet lookup of %d names...", length(todo)))
  # Short timeout per request (the global one is long for the big downloads), up to 3 tries.
  h <- curl::new_handle(timeout = 20, connecttimeout = 10)
  for (i in seq_along(todo)) {
    nm <- todo[i]
    # The API reads "/" as a path separator; RefMet accepts "_" between chains.
    q <- URLencode(gsub("/", "_", nm), reserved = TRUE)
    ans <- NULL
    for (try in 1:3) {
      r <- tryCatch(curl::curl_fetch_memory(paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/", q), handle = h),
                    error = function(e) NULL)
      if (!is.null(r) && r$status_code == 200) {
        # An unknown name gives an empty answer or an HTML error page: recorded as "" (not known to RefMet).
        ans <- tryCatch(fromJSON(rawToChar(r$content))$refmet_name, error = function(e) "")
        if (is.null(ans)) ans <- ""
        break
      }
      Sys.sleep(2)
    }
    # Names that still failed (network) are not recorded, so the next run retries them.
    if (!is.null(ans)) refmet_map <<- rbind(refmet_map, data.table(name = nm, refmet_name = ans))
    # Save every 50 names, so an interrupted run keeps its progress.
    if (i %% 50 == 0 || i == length(todo)) { fwrite(refmet_map, REFMET_FILE); message(sprintf("  %d / %d", i, length(todo))) }
    Sys.sleep(0.1)
  }
  refmet_map[.(x), on = "name", refmet_name]
}

# Helper: sum composition of a lipid given with chains ("PC 16:0_18:1" -> "pc 34:1", "TG 54:2_18:1" -> "tg 54:2").
sum_key <- function(x) {
  x <- tolower(x)
  out <- rep(NA_character_, length(x))
  # Two chains: add carbons and double bonds.
  m <- regmatches(x, regexec("^([a-z]+) ([op]-)?([0-9]+):([0-9]+)_([0-9]+):([0-9]+)$", x))
  two <- lengths(m) == 7
  out[two] <- vapply(m[two], function(v) sprintf("%s %s%d:%d", v[2], v[3], as.integer(v[4]) + as.integer(v[6]),
                                                  as.integer(v[5]) + as.integer(v[7])), "")
  # TG/DG with total and one named chain: keep the total.
  tg <- !two & grepl("^(tg|dg) [0-9]+:[0-9]+_[0-9]+:[0-9]+$", x)
  out[tg] <- sub("_.*$", "", x[tg])
  out
}

# Helper: our metabolites matched by a vector of names from a dataset.
match_mets <- function(names) {
  names <- unique(trimws(names[!is.na(names) & names != ""]))
  rm <- refmet(names)
  # Candidate keys: the raw name, its RefMet name, and the sum composition of either.
  keys <- unique(tolower(c(names, rm[!is.na(rm) & rm != ""], na.omit(sum_key(c(names, rm[!is.na(rm) & rm != ""]))))))
  unique(met_key[key %in% keys, .(feature = metabolite, id = refmet_id, in_network)])
}

# Helper: first sheet of an .xlsx as a data.table (header = first row), using xml2 only.
read_xlsx1 <- function(f) {
  d <- tempfile(); unzip(f, exdir = d)
  ss <- xml_text(xml_find_all(read_xml(file.path(d, "xl/sharedStrings.xml")), "//d1:si", xml_ns(read_xml(file.path(d, "xl/sharedStrings.xml")))))
  sh <- read_xml(file.path(d, "xl/worksheets/sheet1.xml")); ns <- xml_ns(sh)
  rows <- lapply(xml_find_all(sh, "//d1:row", ns), function(r) {
    cs <- xml_find_all(r, "d1:c", ns)
    v <- vapply(cs, function(c) { t <- xml_text(xml_find_first(c, "d1:v", ns)); if (!is.na(t) && identical(xml_attr(c, "t"), "s")) ss[as.integer(t) + 1] else t }, "")
    # Place values by column letter so empty cells keep their position.
    names(v) <- gsub("[0-9]", "", xml_attr(cs, "r")); v
  })
  cols <- names(rows[[1]])
  dt <- rbindlist(lapply(rows[-1], function(v) as.list(setNames(v[cols], cols))))
  setnames(dt, unname(rows[[1]]))
  dt
}

# ---- collect the datasets ----------------------------------------------------------------------------
hits <- list()   # one data.table of matched features per dataset x ome
info <- list()   # one row of description per dataset x ome
add <- function(program, dataset, ome, tissue, n_subjects, n_samples, found, n_measured) {
  k <- paste(dataset, ome)
  hits[[k]] <<- cbind(data.table(program, dataset, ome), found)
  info[[k]] <<- data.table(program, dataset, ome, tissue, n_subjects, n_samples, n_features_in_dataset = n_measured)
}

# HMP: iHMP T2D / prediabetes study (processed tables of the samples CFDE indexes; see header).
IHMP <- "iHMP T2D prediabetes (HMP project T2D)"
f <- fetch(paste0(IPOP, "/HMP/RNAseq_abundance/RNAseq_abundance.txt"))
x <- fread(f, select = 1L); h <- names(fread(f, nrows = 0))[-1]
add("HMP", IHMP, "transcriptomics", "blood (PBMC)", uniqueN(subject_of(x[[1]])), nrow(x), match_genes(h, "symbol"), length(h))
f <- fetch(paste0(IPOP, "/HMP/proteome_abundance/proteome_abundance.txt"))
x <- fread(f, select = 1L); h <- names(fread(f, nrows = 0))[-1]
add("HMP", IHMP, "proteomics", "plasma", uniqueN(subject_of(x[[1]])), nrow(x),
    match_genes(sub("([._][0-9]+)$", "", h), "symbol"), length(h))
f <- fetch(paste0(IPOP, "/HMP/metabolome_abundance/metabolome_abundance.txt"))
x <- fread(f, select = 1L)
ann <- read_xlsx1(fetch(paste0(IPOP, "/HMP/metabolome_abundance/iPOP_Metablolite_Annotation.xlsx")))
# Named features only; "(1)", "(2)" mark isomers of one name.
nm <- unique(sub("\\([0-9]+\\)$", "", ann$Metabolite[!is.na(ann$Metabolite)]))
add("HMP", IHMP, "metabolomics", "plasma", uniqueN(subject_of(x[[1]])), nrow(x), match_mets(nm), length(nm))

# GTEx v8: genes in the RNA-seq table, and those expressed (median TPM > 1) in any of our three tissues.
g <- fread(fetch(GTEX), skip = 2, select = c("Name", GTEX_TISSUES))
g[, expressed := do.call(pmax, .SD) > 1, .SDcols = GTEX_TISSUES]
add("GTEx", "GTEx v8", "transcriptomics", "adipose, skeletal muscle, whole blood (+ all tissues)", 948, 17382,
    match_genes(g$Name, "ensembl"), nrow(g))
add("GTEx", "GTEx v8 (median TPM > 1 in adipose, muscle or blood)", "transcriptomics", "adipose, skeletal muscle, whole blood",
    948, 17382, match_genes(g[expressed == TRUE, Name], "ensembl"), g[, sum(expressed)])

# Metabolomics Workbench: the human studies ranked by search_mw_diabetes.R.
mw <- fread(file.path(DIAB, "mw_diabetes_studies_ranked.csv"))[species == "Homo sapiens"]
mwm <- fread(file.path(DIAB, "mw_diabetes_study_metabolites.csv"))[study_id %in% mw$study_id]
best <- mw[1]
add("Metabolomics Workbench", sprintf("%s (best single study: %s)", best$study_id, substr(best$study_title, 1, 60)), "metabolomics",
    "serum", NA, best$number_of_samples, mets[metabolite %in% mwm[study_id == best$study_id, metabolite], .(feature = metabolite, id = refmet_id, in_network)],
    best$n_named)
add("Metabolomics Workbench", sprintf("union of %d human diabetes-related studies", nrow(mw)), "metabolomics",
    "various (mostly plasma/serum)", NA, sum(as.numeric(mw$number_of_samples), na.rm = TRUE),
    mets[metabolite %in% mwm$metabolite, .(feature = metabolite, id = refmet_id, in_network)], NA)

# ---- score -------------------------------------------------------------------------------------------
feat <- rbindlist(hits, fill = TRUE)
feat <- unique(feat, by = c("dataset", "ome", "id"))
# Denominators per ome.
den <- data.table(ome = c("transcriptomics", "proteomics", "metabolomics"),
                  network_total = c(471L, 471L, 450L),
                  all_total = c(genes[, sum(as_rna)], genes[, sum(as_prot)], nrow(mets)))
score <- feat[, .(n_network = sum(in_network),
                  n_all = if (ome[1] == "transcriptomics") sum(as_rna, na.rm = TRUE)
                          else if (ome[1] == "proteomics") sum(as_prot, na.rm = TRUE) else .N), by = .(dataset, ome)]
ranking <- den[score[rbindlist(info), on = .(dataset, ome)], on = "ome"]
ranking[, `:=`(pct_network = round(100 * n_network / network_total, 1), pct_all = round(100 * n_all / all_total, 1))]
# Number of omes the same people were measured in (datasets sharing a cohort).
ranking[, omes_same_people := fifelse(program == "HMP", uniqueN(ome), 1L), by = .(program, sub(" \\(.*", "", dataset))]
# Where to find each dataset: its landing page, and where its data files can be downloaded.
mw_page <- "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID="
ranking[, `:=`(url = fcase(
  dataset == IHMP, "https://data.cfde.cloud/matrix/HMP (C2M2 package; raw files at portal.hmpdacc.org, project T2D)",
  startsWith(dataset, "GTEx"), "https://gtexportal.org/home/downloads/adult-gtex",
  startsWith(dataset, best$study_id), paste0(mw_page, best$study_id),
  startsWith(dataset, "union"), "https://www.metabolomicsworkbench.org"),
  # Processed tables used for the counts (the iPOP server has no folder listings, so each row gives its exact file(s)).
  data_url = fcase(
  dataset == IHMP & ome == "transcriptomics", paste0(IPOP, "/HMP/RNAseq_abundance/RNAseq_abundance.txt"),
  dataset == IHMP & ome == "proteomics", paste0(IPOP, "/HMP/proteome_abundance/proteome_abundance.txt"),
  dataset == IHMP & ome == "metabolomics",
    paste0(IPOP, "/HMP/metabolome_abundance/metabolome_abundance.txt ; ", IPOP, "/HMP/metabolome_abundance/iPOP_Metablolite_Annotation.xlsx"),
  startsWith(dataset, "GTEx"), "https://gtexportal.org/home/downloads/adult-gtex (T2D donor flag: dbGaP phs000424)",
  startsWith(dataset, best$study_id), paste0(mw_page, best$study_id, " (\"Download named metabolite data\")"),
  startsWith(dataset, "union"), "one link per study: url column of mw_diabetes_studies_ranked.csv"))]
setorder(ranking, ome, -n_network)
setcolorder(ranking, c("ome", "program", "dataset", "tissue", "n_subjects", "n_samples", "n_features_in_dataset",
                       "n_network", "network_total", "pct_network", "n_all", "all_total", "pct_all", "omes_same_people",
                       "url", "data_url"))
fwrite(ranking, file.path(DIAB, "cfde_program_ranking.csv"))
setorder(feat, ome, dataset, -in_network, feature)
fwrite(feat[, .(ome, program, dataset, feature, id, in_network)], file.path(DIAB, "cfde_program_features.csv"))

# ---- status of all 19 programs (from the survey; see header) ------------------------------------------
status <- fread(text = "program|counted|human omics|diabetes-relevant data|access
HMP|yes|iHMP T2D study: PBMC RNA-seq; plasma proteome; plasma metabolome (same people). In CFDE as raw files only; processed tables from the study team|iHMP T2D / prediabetes: ~105 people, insulin resistant vs sensitive (SSPG), A1C, glucose|open (CC0)
GTEx|yes|bulk RNA-seq in ~50 tissues incl. adipose, muscle, blood, pancreas|T2D donor flag exists but is dbGaP-protected|expression open; phenotypes dbGaP phs000424
Metabolomics Workbench|yes|metabolomics / lipidomics|many diabetes/insulin studies (search_mw_diabetes.R), usually one ome each|open
Bridge2AI|no|none released yet (blood banked for future omics)|AI-READI: 2,280 people stratified by T2D severity; CGM, A1c|registered / controlled
NPH|no|metabolomics, proteomics, microbiome planned|diet-response study with CGM (8,000-10,000 people)|All of Us Researcher Workbench
A2CPS|no|blood omics planned, not in current release|diabetes only as a comorbidity item|NDA, no omics yet
SenNet|no|single-cell / spatial|pancreas and adipose donors, no diabetes annotation|open
HuBMAP|no|imaging, ATAC, some single-cell|pancreas donors are non-diabetic|open
exRNA|no|plasma small RNA|no diabetes study in the open atlas|open
LINCS|no|cell-line perturbation signatures|only via drug/insulin signatures, no people|open
IMPC|no|mouse only|glucose tolerance in knockout mice (gene-level check only)|open
Kids First|no|WGS, some RNA|paediatric cancer / birth defects|controlled
MoTrPAC|no|(our own source data)|—|—
SPARC|no|nerve, mostly animal|—|open
4DN|no|chromatin structure|—|open
SCGE|no|gene-editing tools|—|open
SMaHT|no|somatic mosaicism DNA|—|controlled
IDG|no|knowledge base (drug targets)|—|open
GlyGen|no|knowledge base (glycoproteins)|—|open", sep = "|")
stopifnot(nrow(status) == 19)
fwrite(status, file.path(DIAB, "cfde_program_status.csv"))

# ---- report ------------------------------------------------------------------------------------------
print(ranking[, .(ome, program, dataset = substr(dataset, 1, 45), tissue = substr(tissue, 1, 22), n_subjects,
                  network = sprintf("%d/%d (%s%%)", n_network, network_total, pct_network),
                  all = sprintf("%d/%d", n_all, all_total), omes_same_people)], nrows = 50)
message("written to ", DIAB)
