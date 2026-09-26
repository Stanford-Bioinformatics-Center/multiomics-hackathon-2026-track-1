# Reproduce an exploratory screen of public MoTrPAC aggregate results.
# Base R only. No participant-level data or model refitting.
# The official load_differential_analysis() returns these stored objects unchanged
# unless optional annotation merging is requested; that merge is explicit below.
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sha <- "535b4044e7417413de471104c619120337602b77"
base <- paste0("https://raw.githubusercontent.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/", sha, "/data/")
objects <- c("HUMAN_FEATURE_TO_GENE", "BLOOD_PROT_OL_DA",
             "MUSCLE_TRNSCRPT_DA", "ADIPOSE_TRNSCRPT_DA",
             "ADIPOSE_PROT_PH_DA", "BLOOD_METAB_DA", "MUSCLE_METAB_DA")
for (obj in objects) {
  file <- file.path(outdir, paste0(obj, ".rda"))
  # Always download the pinned version, never silently reuse another release.
  download.file(paste0(base, obj, ".rda"), file, mode = "wb", quiet = TRUE)
}
read_object <- function(obj) {
  env <- new.env(parent = emptyenv())
  loaded <- load(file.path(outdir, paste0(obj, ".rda")), envir = env)
  stopifnot(obj %in% loaded)
  as.data.frame(env[[obj]])
}
mapping <- read_object("HUMAN_FEATURE_TO_GENE")
atlas_map <- data.frame(
  atlas_name = c("FSTL1", "ELA", "NRG-1", "Fractalkine", "HSP72", "HSP72",
    "TGF-beta2", "Omentin", "GDF15", "Follistatin", "Decorin", "Myonectin",
    "Lubricin", "VEGF", "FGF2", "SeP", "BDNF", "Cathepsin B", "GPLD1",
    "Clusterin", "PF4", "VIP", "GDNF", "Apelin", "Adiponectin", "IL-6",
    "Irisin", "IGF-1", "Klotho"),
  gene_symbol = c("FSTL1", "APELA", "NRG1", "CX3CL1", "HSPA1A", "HSPA1B",
    "TGFB2", "ITLN1", "GDF15", "FST", "DCN", "ERFE", "PRG4", "VEGFA",
    "FGF2", "SELENOP", "BDNF", "CTSB", "GPLD1", "CLU", "PF4", "VIP",
    "GDNF", "APLN", "ADIPOQ", "IL6", "FNDC5", "IGF1", "KL"))
write.csv(atlas_map, file.path(outdir, "atlas_protein_gene_mapping.csv"), row.names = FALSE)
genes <- unique(c(atlas_map$gene_symbol, "CCN1", "MME", "PLAT", "CXCL12",
  "CX3CR1", "CXCR4", "ACKR3", "IL6R", "IL6ST", "SOCS3", "STAT3",
  "BCKDHA", "BCKDHB", "BCKDK", "PPM1K", "BCAT2", "GLUL", "AHNAK2",
  "SPTBN1", "SLC9A3R1", "PPARGC1A", "AADAT", "KYAT1", "KYAT3", "GOT2"))
keep_contrast <- c("exercise_with_controls", "Endur_vs_Resist")
cols <- c("tissue", "assay", "feature_id", "gene_symbol", "platform",
  "contrast_type", "contrast_category", "Timepoint", "logFC",
  "CI.L_calculated", "CI.R_calculated", "p_value", "adj_p_value")
pieces <- list()
for (obj in c("BLOOD_PROT_OL_DA", "MUSCLE_TRNSCRPT_DA",
              "ADIPOSE_TRNSCRPT_DA", "ADIPOSE_PROT_PH_DA")) {
  dat <- read_object(obj)
  dat <- dat[dat$contrast_type %in% keep_contrast, ]
  mp <- mapping[mapping$assay %in% unique(dat$assay) & mapping$gene_symbol %in% genes,
                c("assay", "feature_id", "gene_symbol")]
  mp <- unique(mp)
  dat <- merge(dat, mp, by = c("assay", "feature_id"))
  dat$platform <- if ("platform" %in% names(dat)) dat$platform else NA_character_
  pieces[[obj]] <- dat[, cols]
}
metabolites <- c("Kynurenic acid", "Kynurenine", "Leucine", "Isoleucine",
                 "Valine", "Ketoleucine", "3-Hydroxybutyric acid", "Lactic acid")
for (obj in c("BLOOD_METAB_DA", "MUSCLE_METAB_DA")) {
  dat <- read_object(obj)
  dat <- dat[dat$contrast_type %in% keep_contrast & dat$feature_id %in% metabolites, ]
  dat$gene_symbol <- NA_character_
  pieces[[obj]] <- dat[, cols]
}
result <- do.call(rbind, pieces)
rownames(result) <- NULL
stopifnot(all(result$adj_p_value >= 0 & result$adj_p_value <= 1))
write.csv(result, file.path(outdir, "candidate_evidence.csv"), row.names = FALSE)
plasma <- pieces[["BLOOD_PROT_OL_DA"]]
plasma <- plasma[plasma$contrast_type == "exercise_with_controls", ]
coverage <- merge(atlas_map, unique(plasma[, c("gene_symbol", "feature_id")]), all.x = TRUE)
coverage$status <- ifelse(is.na(coverage$feature_id), "absent_from_analyzed_table", "present")
write.csv(coverage, file.path(outdir, "atlas_plasma_coverage.csv"), row.names = FALSE)
writeLines(c(paste("MoTrPAC commit:", sha),
  paste("Retrieved UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Source FDR values retained: BH within tissue, assay, platform and contrast.",
  "No new omnibus FDR, enrichment, mediation, or causal analysis performed.",
  "Atlas aliases are manual; precursors and mature secreted products are not interchangeable.",
  "A missing assay is not evidence that an exerkine is absent biologically."),
  file.path(outdir, "analysis_provenance.txt"))
cat("Exported", nrow(result), "candidate rows.\n")
