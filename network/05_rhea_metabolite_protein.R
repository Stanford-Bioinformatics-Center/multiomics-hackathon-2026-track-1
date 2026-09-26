#!/usr/bin/env Rscript
# =====================================================================================================
# 05_rhea_metabolite_protein.R — STEP 5: WHICH OF OUR METABOLITES ARE USED BY WHICH OF OUR PROTEINS (RHEA)
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Links each of our 450 metabolites to the proteins (among our 471 genes) that act on it as enzymes.
#   A link means: in the Rhea database, the metabolite is a substrate or product of a reaction that the
#   protein catalyses. This is enzyme-substrate evidence (the protein physically handles the molecule
#   during catalysis). It is not a general binding assay: transporters, carriers and regulatory binding
#   are only included where Rhea records them as a reaction. Step 6 uses these links to decide which
#   metabolites are connected to each other.
#
# THE DATABASE: Rhea (www.rhea-db.org), expert-curated biochemical reactions (SIB / EMBL-EBI)
#   Downloaded from ftp.expasy.org/databases/rhea/ (release number and date are recorded in the output).
#   Licence CC BY 4.0. Files used:
#     txt/rhea-kegg.reaction.gz     every reaction and its participants as ChEBI IDs ("EQUATION" lines)
#     tsv/rhea2uniprot_sprot.tsv    reactions -> reviewed (Swiss-Prot) UniProt enzymes, all organisms
#     tsv/chebi_pH7_3_mapping.tsv   ChEBI ID -> the form of the molecule that dominates at pH 7.3
#     rhea-release.properties       release number and date
#   Upstream curation: Rhea curators balance every reaction and describe molecules as they exist at
#   physiological pH (e.g. "L-lactate", not "L-lactic acid"); enzyme links come from UniProtKB/Swiss-Prot
#   curation.
#
# HOW OUR METABOLITES ARE MATCHED TO RHEA
#   1. Our ChEBI IDs come from step 1c (column chebi_all: every ChEBI ID found for the metabolite). Most
#      are the neutral molecule. Each is translated to its pH 7.3 form with Rhea's own mapping, and both
#      the original and the translated IDs are used, so "L-lactic acid" (CHEBI:422) finds reactions
#      written with "L-lactate" (CHEBI:16651).
#   2. Metabolites with no ChEBI ID (mostly lipid species without chain positions) cannot be matched.
#      Rhea's generic lipid entries (e.g. "a 1,2-diacyl-sn-glycero-3-phosphocholine") are NOT expanded to
#      our specific species here: no ontology-based guessing.
#
# HOW PROTEINS ARE MATCHED
#   Only our 471 genes count (the team's rule). Their UniProt accessions come from the package's lookup
#   table (HUMAN_FEATURE_TO_GENE, MS and OLINK protein entries, isoform suffix removed), exactly as in
#   step 2. Because the accessions are human, only human enzymes can match, although Rhea covers all
#   organisms.
#
# NO HUB REMOVAL HERE (the team's decision: build first, decide later). Very common molecules (water,
#   protons, ATP, NAD+, CoA) take part in thousands of reactions and will link to many proteins; they are
#   flagged by step 7, not removed.
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis (gene -> UniProt lookup). Internet only the first
#   time (files are cached in $HACK_EXT).
#
# INPUTS:  $HACK_OUT/01c_metabolite_ids.csv (ChEBI IDs), $HACK_OUT/01_nodes_EE.csv (the 471 genes),
#          Rhea files (downloaded to $HACK_EXT, default ~/Desktop/output/hackathon-2026-track1/external/rhea)
# OUTPUTS: $HACK_OUT/05_metabolite_protein_links.csv   one row per metabolite-protein link
#          $HACK_OUT/05_rhea_summary.csv               counts, coverage, Rhea release
# =====================================================================================================

# Load the packages quietly: data.table (tables) and the MoTrPAC package (gene -> protein IDs).
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Results folder (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Folder for downloaded external files (override with HACK_EXT); kept outside the code repo.
EXT <- Sys.getenv("HACK_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/rhea"))
# Create it if needed.
dir.create(EXT, recursive = TRUE, showWarnings = FALSE)
# Base address of Rhea's download site.
RHEA_FTP <- "https://ftp.expasy.org/databases/rhea/"
# The files we need (path on the site -> local file name).
FILES <- c("txt/rhea-kegg.reaction.gz", "tsv/rhea2uniprot_sprot.tsv", "tsv/chebi_pH7_3_mapping.tsv",
           "rhea-release.properties")

# ---- download (only files not already cached) ---------------------------------------------------
# Download each file once; later runs reuse the cached copy so results stay reproducible.
for (f in FILES) {
  # the local path for this file
  dest <- file.path(EXT, basename(f))
  # download only if it is not already there (mode "wb" keeps the gzip file intact)
  if (!file.exists(dest)) download.file(paste0(RHEA_FTP, f), dest, mode = "wb", quiet = TRUE)
}
# Read the release number and date (lines like "rhea.release.number=142").
rel <- readLines(file.path(EXT, "rhea-release.properties"))
# Keep just the values.
rhea_release <- sub(".*=", "", grep("release.number", rel, value = TRUE))
# The release date, the same way.
rhea_date <- sub(".*=", "", grep("release.date", rel, value = TRUE))
# Report which release is being used.
message(sprintf("Rhea release %s (%s)", rhea_release, rhea_date))

# ---- Rhea reactions -> participants (ChEBI) -------------------------------------------------------
# Read the reaction file line by line (it is small enough to hold in memory).
rx <- readLines(gzfile(file.path(EXT, "rhea-kegg.reaction.gz")))
# Keep the lines that start a reaction ("ENTRY RHEA:10000") or list its molecules ("EQUATION ...").
keep <- grepl("^(ENTRY|EQUATION)", rx)
# Drop all other lines.
rx <- rx[keep]
# For every EQUATION line, the reaction it belongs to is the most recent ENTRY line above it.
entry_idx <- cumsum(grepl("^ENTRY", rx))
# The reaction number from each ENTRY line ("RHEA:10000" -> "10000").
entries <- sub("^ENTRY\\s+RHEA:", "", rx[grepl("^ENTRY", rx)])
# The equation lines and their reaction IDs.
eq <- data.table(rhea_id = entries[entry_idx[grepl("^EQUATION", rx)]], line = rx[grepl("^EQUATION", rx)])
# Pull every ChEBI ID out of each equation (both sides; coefficients like "2 CHEBI:15378" are ignored).
parts <- eq[, .(chebi = regmatches(line, gregexpr("CHEBI:[0-9]+", line))[[1]]), by = rhea_id]
# One row per reaction and participant.
parts <- unique(parts)
# Report the size.
message(sprintf("Rhea: %d reactions with participants, %d distinct molecules",
                uniqueN(parts$rhea_id), uniqueN(parts$chebi)))

# ---- Rhea reactions -> enzymes (UniProt) ----------------------------------------------------------
# Reviewed enzyme links for every reaction ID (undirected and directed versions are both listed).
enz <- fread(file.path(EXT, "rhea2uniprot_sprot.tsv"), colClasses = "character")
# Keep reaction ID and protein accession.
enz <- unique(enz[, .(rhea_id = RHEA_ID, uniprot = ID)])

# ---- our 471 genes -> UniProt --------------------------------------------------------------------
# The 471 gene IDs from step 1.
genes <- fread(file.path(OUT, "01_nodes_EE.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol)]
# Their UniProt accessions from the package (MS and OLINK protein entries), isoform suffix removed.
acc <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay %in% c("prot-pr", "prot-ol")][
  , .(entrez_gene = as.character(entrez_gene), uniprot = sub("-[0-9]+$", "", as.character(uniprot)))][
  entrez_gene %in% genes$entrez_gene & !is.na(uniprot)])
# Safety check: every gene has at least one accession (same as step 2).
stopifnot(uniqueN(acc$entrez_gene) == 471)

# ---- our 450 metabolites -> ChEBI (original + pH 7.3 form) ----------------------------------------
# Step 1c's identifier table.
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
# One row per metabolite and ChEBI ID (a metabolite can have several, separated by ";").
mchebi <- ids[chebi_all != "" & !is.na(chebi_all), .(chebi = unlist(strsplit(chebi_all, ";"))), by = metabolite]
# Rhea's mapping to the dominant form at pH 7.3 (numbers only in the file, so add the "CHEBI:" prefix).
ph <- fread(file.path(EXT, "chebi_pH7_3_mapping.tsv"), colClasses = "character")[
  , .(chebi = paste0("CHEBI:", CHEBI), chebi_ph73 = paste0("CHEBI:", CHEBI_PH7_3))]
# Add the pH 7.3 form of each of our IDs (where the mapping has one).
mchebi <- merge(mchebi, ph, by = "chebi", all.x = TRUE)
# Use both the original and the pH 7.3 ID for matching.
mmatch <- unique(rbind(mchebi[, .(metabolite, chebi, via = "as_is")],
                       mchebi[!is.na(chebi_ph73) & chebi_ph73 != chebi, .(metabolite, chebi = chebi_ph73, via = "pH7.3_form")]))

# ---- links: metabolite takes part in a reaction catalysed by one of our proteins -----------------
# Metabolite -> reactions it takes part in.
m_rx <- merge(mmatch, parts, by = "chebi", allow.cartesian = TRUE)
# Reaction -> our proteins that catalyse it.
rx_p <- merge(enz, acc, by = "uniprot", allow.cartesian = TRUE)
# Metabolite -> protein, through a shared reaction.
link <- merge(m_rx, rx_p, by = "rhea_id", allow.cartesian = TRUE)
# One row per metabolite-gene pair: how many reactions support it, example reactions, and matching route.
link <- link[, .(n_reactions = uniqueN(rhea_id),
                 example_reactions = paste(head(sort(unique(paste0("RHEA:", rhea_id))), 5), collapse = ";"),
                 uniprot = paste(sort(unique(uniprot)), collapse = ";"),
                 matched_via = paste(sort(unique(via)), collapse = "+")),
             by = .(metabolite, entrez_gene)]
# Add gene symbols and the metabolite's class (from step 1c).
link <- genes[link, on = "entrez_gene"]
# Add each metabolite's RefMet classes.
link <- ids[, .(metabolite, super_class, main_class)][link, on = "metabolite"]
# Order: metabolite, then gene symbol.
setcolorder(link, c("metabolite", "super_class", "main_class", "entrez_gene", "gene_symbol", "uniprot",
                    "n_reactions", "example_reactions", "matched_via"))
# Sort rows by metabolite, then gene.
setorder(link, metabolite, gene_symbol)
# Save.
fwrite(link, file.path(OUT, "05_metabolite_protein_links.csv"))

# ---- summary -------------------------------------------------------------------------------------
summ <- data.table(
  metric = c("rhea_release", "rhea_release_date", "metabolites", "metabolites_with_chebi",
             "metabolites_found_in_rhea", "metabolites_linked_to_our_genes", "genes", "genes_linked_to_our_metabolites",
             "metabolite_gene_links"),
  value = c(rhea_release, rhea_date, nrow(ids), uniqueN(mchebi$metabolite),
            uniqueN(m_rx$metabolite), uniqueN(link$metabolite), 471, uniqueN(link$entrez_gene), nrow(link)))
# Save and show.
fwrite(summ, file.path(OUT, "05_rhea_summary.csv"))
# Show the summary on screen.
print(summ)
