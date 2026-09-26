#!/usr/bin/env Rscript
# =====================================================================================================
# 09_export_resources.R — STEP 9: SHAREABLE LISTS OF OUR 471 PROTEINS AND 450 METABOLITES
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Writes two small reference tables into the repository folder network/resource/, so teammates can
#   look up our network's features in other datasets without running the pipeline:
#     resource/proteins_471.csv     the 471 genes/proteins (measured as RNA and protein in all three tissues)
#     resource/metabolites_450.csv  the 450 metabolites (measured in adipose, blood and muscle)
#   Each table carries the identifiers most other resources use, so a feature can be matched by whichever
#   ID the other dataset has. Unlike every other step, these files ARE committed to the repo on purpose:
#   they are small reference lists, not results. They contain names and IDs only, no measurements.
#
# WHERE THE IDENTIFIERS COME FROM
#   Proteins: Entrez gene ID and symbol (step 1); UniProt accessions (package table HUMAN_FEATURE_TO_GENE,
#     MS and OLINK protein entries, isoform suffix removed, as in step 2); Ensembl gene IDs (same table,
#     RNA-seq entries); whether the protein is in the STRING file and its number of network partners (step 2).
#   Metabolites: RefMet name, RefMet ID, classes, PubChem CID, InChIKey and ChEBI IDs (step 1c; looked up
#     2026-09-26); KEGG compound ID (package table HUMAN_FEATURE_TO_GENE, where available); the assay
#     platform that measured it in each tissue (step 1b).
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis.
#
# INPUTS:  $HACK_OUT/01_nodes_EE.csv, 02_nodes_string.csv, 01c_metabolite_ids.csv, 01b_metab_nodes_provenance.csv
# OUTPUTS: network/resource/proteins_471.csv and network/resource/metabolites_450.csv (inside the repo)
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Results folder where the pipeline's outputs live (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# This script's own folder, found from how it was launched (falls back to the working directory).
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
                 error = function(e) getwd())
# The resource folder inside the repo, next to this script.
RES <- file.path(HERE, "resource")
# Create it if needed.
dir.create(RES, showWarnings = FALSE)
# The package's feature-to-ID lookup table, all columns as text.
map <- as.data.table(HUMAN_FEATURE_TO_GENE)[, lapply(.SD, as.character)]

# ---- proteins -----------------------------------------------------------------------------------
# The 471 genes (ID and symbol) from step 1.
prot <- fread(file.path(OUT, "01_nodes_EE.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol)]
# UniProt accessions per gene (MS and OLINK entries, isoform suffix removed; several joined by ";").
uni <- unique(map[assay %in% c("prot-pr", "prot-ol") & entrez_gene %in% prot$entrez_gene & !is.na(uniprot),
                  .(entrez_gene, uniprot = sub("-[0-9]+$", "", uniprot))])[
  , .(uniprot = paste(sort(unique(uniprot)), collapse = ";")), by = entrez_gene]
# Ensembl gene IDs per gene (RNA-seq entries; several joined by ";").
ens <- unique(map[assay == "transcript-rna-seq" & entrez_gene %in% prot$entrez_gene & !is.na(ensembl_gene),
                  .(entrez_gene, ensembl_gene)])[
  , .(ensembl_gene = paste(sort(unique(ensembl_gene)), collapse = ";")), by = entrez_gene]
# STRING presence and number of partners in our network (step 2).
str2 <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[
  , .(entrez_gene, in_string, network_degree = degree)]
# Put the columns together, one row per gene.
prot <- str2[ens[uni[prot, on = "entrez_gene"], on = "entrez_gene"], on = "entrez_gene"]
# Column order and alphabetical row order.
setcolorder(prot, c("gene_symbol", "entrez_gene", "uniprot", "ensembl_gene", "in_string", "network_degree"))
# Rows in alphabetical order of gene symbol.
setorder(prot, gene_symbol)
# Safety check: 471 genes, each with a UniProt and an Ensembl ID.
stopifnot(nrow(prot) == 471, !anyNA(prot$uniprot), !anyNA(prot$ensembl_gene))
# Save.
fwrite(prot, file.path(RES, "proteins_471.csv"))

# ---- metabolites --------------------------------------------------------------------------------
# Identifiers and classes from step 1c.
met <- fread(file.path(OUT, "01c_metabolite_ids.csv"), colClasses = list(character = "pubchem_cid"))[
  , .(metabolite, refmet_id, super_class, main_class, pubchem_cid, inchi_key, chebi_id, chebi_all)]
# KEGG compound IDs from the package where available (several joined by ";").
kegg <- unique(map[assay == "metab" & feature_id %in% met$metabolite & !is.na(kegg_id), .(metabolite = feature_id, kegg_id)])[
  , .(kegg_id = paste(sort(unique(kegg_id)), collapse = ";")), by = metabolite]
# The platform that measured each metabolite in each tissue (step 1b).
plat <- fread(file.path(OUT, "01b_metab_nodes_provenance.csv"))
# Put the columns together, one row per metabolite.
met <- plat[kegg[met, on = "metabolite"], on = "metabolite"]
# Column order and alphabetical row order (ignoring upper/lower case).
setcolorder(met, c("metabolite", "refmet_id", "super_class", "main_class", "chebi_id", "chebi_all",
                   "pubchem_cid", "inchi_key", "kegg_id", "platform_adipose", "platform_blood", "platform_muscle"))
# Rows in alphabetical order of metabolite name (ignoring upper/lower case).
met <- met[order(tolower(metabolite))]
# Safety check: 450 metabolites, each with its three platforms.
stopifnot(nrow(met) == 450, !anyNA(met$platform_muscle))
# Save.
fwrite(met, file.path(RES, "metabolites_450.csv"))

# Report what was written and how complete each identifier is.
message(sprintf("proteins_471.csv: %d rows (UniProt %d, Ensembl %d, in STRING %d)", nrow(prot),
                sum(prot$uniprot != ""), sum(prot$ensembl_gene != ""), sum(prot$in_string)))
# The same for the metabolite table.
message(sprintf("metabolites_450.csv: %d rows (RefMet ID %d, ChEBI %d, PubChem %d, InChIKey %d, KEGG %d)", nrow(met),
                sum(met$refmet_id != "" & !is.na(met$refmet_id)), sum(met$chebi_id != "" & !is.na(met$chebi_id)),
                sum(met$pubchem_cid != "" & !is.na(met$pubchem_cid)), sum(met$inchi_key != "" & !is.na(met$inchi_key)),
                sum(!is.na(met$kegg_id))))
