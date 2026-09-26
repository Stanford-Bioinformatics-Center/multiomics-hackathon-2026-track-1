#!/usr/bin/env Rscript
# =====================================================================================================
# network/inventory/export_phospho_features.R — ONE CSV OF ALL MOTRPAC PHOSPHOPROTEOMIC FEATURE IDS
# =====================================================================================================
# PURPOSE: a shareable list of every MoTrPAC human phosphosite feature (muscle and adipose) with its
#   identifiers and site annotation, flagged for our 471 network proteins, for planning the phospho layer.
# HOW TO RUN: Rscript network/inventory/export_phospho_features.R   (after glygen_protein_inventory.py,
#   for the "known in GlyGen" flag; without it that column is left empty)
# DATA: MotrpacHumanPreSuspensionAnalysis v0.2.4 (MUSCLE/ADIPOSE_PROT_PH_DA, HUMAN_FEATURE_TO_GENE);
#   MotrpacHumanPreSuspensionData (MUSCLE/ADIPOSE_PROT_PH_QC feature metadata).
# OUTPUT: $HACK_OUT/inventory/motrpac_phospho_feature_ids.csv, one row per feature_id:
#   feature_id           MoTrPAC ID: UniProt accession[-isoform]_<residue><position><lowercase residue>...
#   uniprot, isoform     accession as in the ID, and whether it is an isoform accession (position is on it)
#   sites, n_sites       the site string (e.g. S758s) and number of phosphosites in the feature
#   residues, positions  e.g. "S;T" and "663;665"
#   gene_symbol, entrez_gene, ensembl_gene   from HUMAN_FEATURE_TO_GENE
#   in_muscle, in_adipose                    measured in that tissue's differential analysis
#   confident_site, flanking_sequence, redundant_ids   from the QC feature metadata (muscle, else adipose)
#   in_471               the gene is one of our 471 network proteins
#   known_in_glygen      single-site feature whose canonical accession, position and residue GlyGen lists
# =====================================================================================================

# Load packages quietly.
suppressMessages(library(data.table))
# Folders (outside the repo).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
INV <- file.path(OUT, "inventory"); dir.create(INV, recursive = TRUE, showWarnings = FALSE)

# Features measured per tissue (differential analysis tables) and their QC annotation.
tis <- c(muscle = "MUSCLE", adipose = "ADIPOSE")
meas <- lapply(tis, function(t) { data(list = paste0(t, "_PROT_PH_DA"), package = "MotrpacHumanPreSuspensionAnalysis")
  unique(as.character(get(paste0(t, "_PROT_PH_DA"))$feature_id)) })
ann <- rbindlist(lapply(names(tis), function(n) { data(list = paste0(tis[[n]], "_PROT_PH_QC"), package = "MotrpacHumanPreSuspensionData")
  as.data.table(get(paste0(tis[[n]], "_PROT_PH_QC"))$feature_metadata)[, .(feature_id = id, confident_site, flanking_sequence, redundant_ids, tissue = n)] }))
# One annotation per feature (muscle first, then adipose).
ann <- ann[order(feature_id, tissue != "muscle")][!duplicated(feature_id)][, tissue := NULL]

# All measured features with identifiers.
X <- data.table(feature_id = sort(unique(unlist(meas))))
X[, `:=`(in_muscle = feature_id %in% meas$muscle, in_adipose = feature_id %in% meas$adipose)]
data("HUMAN_FEATURE_TO_GENE", package = "MotrpacHumanPreSuspensionAnalysis")
f2g <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay == "prot-ph", .(feature_id = as.character(feature_id), gene_symbol = as.character(gene_symbol),
                                                                       entrez_gene = as.character(entrez_gene), ensembl_gene = as.character(ensembl_gene))], by = "feature_id")
X <- f2g[X, on = "feature_id"]
# Parse the ID: accession, isoform flag, site string, residues and positions.
X[, `:=`(uniprot = sub("_.*$", "", feature_id), sites = sub("^[^_]*_", "", feature_id))]
X[, isoform := grepl("-[0-9]+$", uniprot)]
st <- regmatches(X$sites, gregexpr("[STY][0-9]+", X$sites))
X[, `:=`(n_sites = lengths(st), residues = sapply(st, function(v) paste(substr(v, 1, 1), collapse = ";")),
         positions = sapply(st, function(v) paste(substring(v, 2), collapse = ";")))]
X <- ann[X, on = "feature_id"]
# Flags: our 471 proteins; site already listed by GlyGen.
genes <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))
X[, in_471 := entrez_gene %in% genes$entrez_gene]
gp <- file.path(INV, "glygen_phosphosites.csv")
if (file.exists(gp)) {
  g <- fread(gp); gk <- paste(sub("-.*$", "", g$glygen_ac), g$position, substr(g$residue, 1, 1))
  X[, known_in_glygen := fifelse(n_sites == 1 & in_471, paste(uniprot, positions, residues) %in% gk, NA)]
} else X[, known_in_glygen := NA]
# Column order and save.
setcolorder(X, c("feature_id", "uniprot", "isoform", "sites", "n_sites", "residues", "positions", "gene_symbol", "entrez_gene", "ensembl_gene",
                 "in_muscle", "in_adipose", "confident_site", "flanking_sequence", "redundant_ids", "in_471", "known_in_glygen"))
stopifnot(!anyDuplicated(X$feature_id), all(X$n_sites >= 1))
fwrite(X, file.path(INV, "motrpac_phospho_feature_ids.csv"))
message(sprintf("-> %s: %d features (muscle %d, adipose %d, both %d; on our 471 proteins %d)", file.path(INV, "motrpac_phospho_feature_ids.csv"),
                nrow(X), sum(X$in_muscle), sum(X$in_adipose), sum(X$in_muscle & X$in_adipose), sum(X$in_471)))
