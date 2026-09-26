#!/usr/bin/env Rscript
# =====================================================================================================
# export_all_features.R — STEP 9b: SHAREABLE LISTS OF EVERY MEASURED GENE/PROTEIN AND METABOLITE
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   export_feature_lists.R lists only the network's features: the 471 genes measured in all six tissue x
#   layer blocks and the 450 metabolites measured in all three tissues. This script lists EVERY feature
#   the MoTrPAC package measured, with columns saying where each one was measured and whether it is in
#   the network, so teammates can look up any feature, not just the network's:
#     genes_all.csv            one row per Entrez gene measured as RNA or protein in any tissue (17,832)
#     features_no_entrez.csv   transcripts/proteins with no Entrez gene ID (not matchable across layers;
#                              listed by their Ensembl or UniProt ID instead)
#     metabolites_all.csv      one row per metabolite measured in any tissue (1,460)
#   Names, identifiers and "measured where" flags only; no measurements. Like step 9, the TABLES ARE NOT
#   COMMITTED: they are written to $HACK_RES, outside the repo.
#
# HOW IT DIFFERS FROM export_feature_lists.R
#   It reads only the package (the "*_DA" result tables and HUMAN_FEATURE_TO_GENE), so no pipeline step
#   has to be run first. It does not include the web-looked-up metabolite IDs (ChEBI, PubChem, InChIKey,
#   classes; step 1c), which exist only for the 450. Phosphoproteomics (PROT_PH) is not included.
#
# WHERE THE IDENTIFIERS COME FROM (package table HUMAN_FEATURE_TO_GENE)
#   Genes: symbol; UniProt accessions (MS and OLINK protein entries, isoform suffix removed, as in step 2);
#     Ensembl gene IDs (RNA-seq entries). Measured-where flags from the six tables used in step 1.
#   Metabolites: RefMet ID and KEGG compound ID. The table's feature_id matches only 651 of the 1,460 names;
#     its refmet_name column matches far more, so both columns are used as lookup keys (checked: no name
#     gets two different RefMet IDs). The platform that measured it in each tissue comes from the *_METAB_DA
#     tables (one platform per metabolite per tissue, the consortium's lowest-CV rule; see step 1b).
#
# CONSISTENCY WITH THE PIPELINE
#   in_network uses the same rules as steps 1 and 1b (Entrez gene in all six blocks; exact RefMet name in all
#   three tissues), and the script stops if those counts are not exactly 471 and 450.
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis (v0.2.4).
#
# INPUTS:  the MoTrPAC package only
# OUTPUTS: $HACK_RES/genes_all.csv, features_no_entrez.csv, metabolites_all.csv
#          ($HACK_RES defaults to $HACK_OUT/resource, i.e. ~/Desktop/output/hackathon-2026-track1/network/resource;
#          never inside the repo)
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Results folder of the pipeline (override with HACK_OUT); used only to place the resource folder.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where the shareable tables are written (override with HACK_RES); outside the repo so no data is committed.
RES <- Sys.getenv("HACK_RES", unset = file.path(OUT, "resource"))
# Create it if needed.
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
# The package's feature-to-ID lookup table, all columns as text.
map <- as.data.table(HUMAN_FEATURE_TO_GENE)[, lapply(.SD, as.character)]

# The three tissues.
TISSUES <- c("adipose", "blood", "muscle")
# The six gene tables, named by block (same tables as step 1: MS protein for adipose/muscle, OLINK for blood).
BLOCKS <- c(rna_adipose = "ADIPOSE_TRNSCRPT_DA", rna_blood = "BLOOD_TRNSCRPT_DA", rna_muscle = "MUSCLE_TRNSCRPT_DA",
            prot_adipose = "ADIPOSE_PROT_PR_DA", prot_blood = "BLOOD_PROT_OL_DA", prot_muscle = "MUSCLE_PROT_PR_DA")

# Helper: join several IDs into one ";"-separated text, sorted, blanks dropped.
collapse_ids <- function(x) paste(sort(unique(x[!is.na(x) & x != ""])), collapse = ";")

# ---- genes -------------------------------------------------------------------------------------------
# Every measured feature in each of the six blocks (one row per block x feature).
feat <- rbindlist(lapply(names(BLOCKS), function(b)
  unique(as.data.table(get(BLOCKS[[b]]))[, .(block = b, assay = as.character(assay), feature_id = as.character(feature_id))])))
# Attach the gene IDs of each feature.
feat <- map[, .(assay, feature_id, entrez_gene, gene_symbol, ensembl_gene, uniprot)][feat, on = .(assay, feature_id)]
# Safety check: every measured feature is in the lookup table.
stopifnot(feat[is.na(gene_symbol) & is.na(ensembl_gene) & is.na(uniprot), .N] == 0)

# One row per Entrez gene with a TRUE/FALSE column per block.
genes <- dcast(unique(feat[!is.na(entrez_gene), .(entrez_gene, block)]), entrez_gene ~ block,
               fun.aggregate = length, value.var = "block")
# Turn counts into TRUE/FALSE, in the fixed block order.
for (b in names(BLOCKS)) genes[, (b) := get(b) > 0]
# Number of blocks the gene is measured in, and whether it is one of the network's 471 (all six).
genes[, n_blocks := rowSums(.SD), .SDcols = names(BLOCKS)][, in_network := n_blocks == length(BLOCKS)]
# Gene symbol (the package gives one per Entrez ID; collapse just in case).
sym <- map[entrez_gene %in% genes$entrez_gene & !is.na(gene_symbol), .(gene_symbol = collapse_ids(gene_symbol)), by = entrez_gene]
# UniProt accessions per gene (MS and OLINK entries, isoform suffix removed), as in export_feature_lists.R.
uni <- map[assay %in% c("prot-pr", "prot-ol") & entrez_gene %in% genes$entrez_gene,
           .(uniprot = collapse_ids(sub("-[0-9]+$", "", uniprot))), by = entrez_gene]
# Ensembl gene IDs per gene (RNA-seq entries).
ens <- map[assay == "transcript-rna-seq" & entrez_gene %in% genes$entrez_gene,
           .(ensembl_gene = collapse_ids(ensembl_gene)), by = entrez_gene]
# Put the columns together, one row per gene.
genes <- ens[uni[sym[genes, on = "entrez_gene"], on = "entrez_gene"], on = "entrez_gene"]
# Column order.
setcolorder(genes, c("gene_symbol", "entrez_gene", "uniprot", "ensembl_gene", names(BLOCKS), "n_blocks", "in_network"))
# Rows in alphabetical order of gene symbol (the ~300 genes the package gives no symbol go last).
genes <- genes[order(is.na(gene_symbol), gene_symbol, entrez_gene)]
# Safety check: one row per gene, and the network's gene count matches step 1.
stopifnot(!anyDuplicated(genes$entrez_gene), sum(genes$in_network) == 471)
# Save.
fwrite(genes, file.path(RES, "genes_all.csv"))

# Features without an Entrez gene ID: one row per feature, with the blocks it was measured in.
noent <- feat[is.na(entrez_gene), .(blocks = paste(names(BLOCKS)[names(BLOCKS) %in% block], collapse = ";")),
              by = .(assay, feature_id, gene_symbol, ensembl_gene, uniprot = sub("-[0-9]+$", "", uniprot))]
# Rows ordered by assay, then feature ID.
setorder(noent, assay, feature_id)
# Save.
fwrite(noent, file.path(RES, "features_no_entrez.csv"))

# ---- metabolites -------------------------------------------------------------------------------------
# Every metabolite in each tissue, with the platform that measured it there.
met <- rbindlist(lapply(TISSUES, function(tis) {
  x <- unique(as.data.table(get(paste0(toupper(tis), "_METAB_DA")))[
    , .(metabolite = as.character(feature_id), platform = as.character(platform))])
  x[, tissue := tis]
}))
# Safety check: one platform per metabolite per tissue (see step 1b).
stopifnot(!anyDuplicated(met[, .(tissue, metabolite)]))
# One row per metabolite, one platform column per tissue (blank = not measured there).
met <- dcast(met, metabolite ~ tissue, value.var = "platform")
# Rename the columns to platform_adipose, platform_blood, platform_muscle.
setnames(met, TISSUES, paste0("platform_", TISSUES))
# Number of tissues it is measured in, and whether it is one of the network's 450 (all three).
met[, n_tissues := rowSums(!is.na(.SD)), .SDcols = paste0("platform_", TISSUES)][, in_network := n_tissues == 3L]
# RefMet and KEGG IDs from the package, looked up by either refmet_name or feature_id (see header).
mm <- map[assay == "metab"]
ids <- unique(rbind(mm[!is.na(refmet_name), .(metabolite = refmet_name, refmet_id, kegg_id)],
                    mm[, .(metabolite = feature_id, refmet_id, kegg_id)]))[metabolite %in% met$metabolite]
ids <- ids[, .(refmet_id = collapse_ids(refmet_id), kegg_id = collapse_ids(kegg_id)), by = metabolite]
# Safety check: no metabolite gets two RefMet IDs.
stopifnot(!any(grepl(";", ids$refmet_id)))
# Attach the IDs (metabolites not in the lookup table get blanks).
met <- ids[met, on = "metabolite"]
# Column order.
setcolorder(met, c("metabolite", "refmet_id", "kegg_id", paste0("platform_", TISSUES), "n_tissues", "in_network"))
# Rows in alphabetical order of metabolite name (ignoring upper/lower case).
met <- met[order(tolower(metabolite))]
# Safety check: the network's metabolite count matches step 1b.
stopifnot(!anyDuplicated(met$metabolite), sum(met$in_network) == 450)
# Save.
fwrite(met, file.path(RES, "metabolites_all.csv"))

# ---- report ------------------------------------------------------------------------------------------
# What was written and how complete each identifier is.
message(sprintf("genes_all.csv: %d genes (in network %d; UniProt %d, Ensembl %d); by number of blocks: %s",
                nrow(genes), sum(genes$in_network), sum(genes$uniprot != "", na.rm = TRUE),
                sum(genes$ensembl_gene != "", na.rm = TRUE),
                paste(sprintf("%d:%d", 1:6, tabulate(genes$n_blocks, 6)), collapse = " ")))
# The features without an Entrez ID.
message(sprintf("features_no_entrez.csv: %d features", nrow(noent)))
# The metabolite table.
message(sprintf("metabolites_all.csv: %d metabolites (in network %d; RefMet ID %d, KEGG %d); by number of tissues: %s",
                nrow(met), sum(met$in_network), sum(met$refmet_id != "", na.rm = TRUE),
                sum(met$kegg_id != "", na.rm = TRUE),
                paste(sprintf("%d:%d", 1:3, tabulate(met$n_tissues, 3)), collapse = " ")))
message("written to ", RES)