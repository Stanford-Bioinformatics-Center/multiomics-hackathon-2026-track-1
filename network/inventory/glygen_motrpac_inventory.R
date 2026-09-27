#!/usr/bin/env Rscript
# =====================================================================================================
# network/inventory/glygen_motrpac_inventory.R — WHICH DATA EXIST FOR OUR 471 PROTEINS AND 450 METABOLITES?
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Before annotating the Neo4j graph with phosphoproteomics (MoTrPAC + GlyGen) and glycosylation / other
#   GlyGen data, count how many of our 471 proteins and 450 metabolites have each data type, so the team can
#   decide what to integrate.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Proteins:
#     - MoTrPAC phosphoproteomics (muscle: 0.5 / 4 / 24 h; adipose: 4 h only; both arms): phosphosites
#       measured per protein and tissue, how many have adj. p < 0.05 in an exercise-vs-control contrast, how
#       many are confidently localised, and how many match a site GlyGen already knows (same canonical
#       protein, position and residue; single-site features only);
#     - GlyGen: the per-protein counts from glygen_protein_inventory.py (every data section).
#   Metabolites:
#     - whether the metabolite is a GlyGen glycan (matched by ChEBI, PubChem compound ID, or the first
#       InChIKey block, i.e. ignoring stereochemistry). GlyGen's human Rhea file lists enzymes only (no small
#       molecules), so metabolite - enzyme links come from our own step 5 (Rhea directly), not from GlyGen.
#   Then a coverage table: for each data type, how many proteins / metabolites have at least one record.
#
# HOW TO RUN
#   After glygen_protein_inventory.py:   Rscript network/inventory/glygen_motrpac_inventory.R
#
# DATA AND PROVENANCE
#   MotrpacHumanPreSuspensionAnalysis v0.2.4 (MUSCLE/ADIPOSE_PROT_PH_DA, HUMAN_FEATURE_TO_GENE) and
#   MotrpacHumanPreSuspensionData (MUSCLE/ADIPOSE_PROT_PH_QC feature metadata); GlyGen release 2.11.1 (small
#   files: glycan_xref_chebi, glycan_xref_pubchem, glycan_sequences_inchi).
#
# TECH STACK:  R 4.4; data.table.
# INPUTS:  $HACK_OUT/02_nodes_string.csv, 01c_metabolite_ids.csv, inventory/glygen_protein_counts.csv,
#          inventory/glygen_phosphosites.csv
# OUTPUTS: $HACK_OUT/inventory/protein_inventory.csv      one row per protein (MoTrPAC phospho + GlyGen counts)
#          $HACK_OUT/inventory/metabolite_inventory.csv   one row per metabolite
#          $HACK_OUT/inventory/coverage_summary.csv       data type, entities with >= 1 record, median records
# =====================================================================================================

# Load packages quietly.
suppressMessages(library(data.table))
# Folders (outside the repo).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
INV <- file.path(OUT, "inventory")
# GlyGen download folder for the few small files used here.
GG <- "https://data.glygen.org/ln2data/releases/data/current/reviewed/"
# Significance threshold for "responds" (as elsewhere in the project).
ALPHA <- 0.05

# ---- our proteins and metabolites ----------------------------------------------------------------------
# Our 471 network proteins (step 2) and 450 metabolites (step 1c); check the counts.
genes <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol, uniprot)]
mets <- fread(file.path(OUT, "01c_metabolite_ids.csv"), colClasses = list(character = c("pubchem_cid")))
stopifnot(nrow(genes) == 471, nrow(mets) == 450)

# ---- MoTrPAC phosphoproteomics per protein ------------------------------------------------------------------
# Load the package's feature-to-gene table and keep the phosphosite rows (feature ID -> Entrez gene).
data("HUMAN_FEATURE_TO_GENE", package = "MotrpacHumanPreSuspensionAnalysis")
f2g <- as.data.table(HUMAN_FEATURE_TO_GENE)[assay == "prot-ph", .(feature_id = as.character(feature_id), entrez_gene = as.character(entrez_gene))]
# GlyGen's known phosphosites on our proteins (canonical accession, position, residue).
gps <- fread(file.path(INV, "glygen_phosphosites.csv"), colClasses = list(character = "entrez_gene"))
# match key: canonical accession (isoform suffix removed) + position + FIRST LETTER of GlyGen's three-letter
# residue (Ser -> S, Thr -> T, but Tyr also -> T, so tyrosine sites can never match; AA below is not used)
gps[, key := paste(sub("-.*$", "", glygen_ac), position, substr(residue, 1, 1))]
# One-letter to three-letter residue codes (S = serine, T = threonine, Y = tyrosine); defined but not used below.
AA <- c(S = "Ser", T = "Thr", Y = "Tyr")
# For each tissue: per phosphosite feature, whether it responds to either arm, plus its confidence and gene.
ph <- rbindlist(lapply(c(muscle = "MUSCLE", adipose = "ADIPOSE"), function(t) {
  # differential analysis: exercise-vs-control (delta-delta) contrasts for both arms
  data(list = paste0(t, "_PROT_PH_DA"), package = "MotrpacHumanPreSuspensionAnalysis")
  da <- as.data.table(get(paste0(t, "_PROT_PH_DA")))[contrast_type == "exercise_with_controls"]
  # arm from the contrast name: Endurance... = EE, otherwise RE
  da[, arm := fifelse(grepl("^Endur", contrast_short), "EE", "RE")]
  # per feature: significant (adj. p < 0.05) in any endurance contrast? any resistance contrast? number of contrasts
  s <- da[, .(sig_EE = any(adj_p_value < ALPHA & arm == "EE"), sig_RE = any(adj_p_value < ALPHA & arm == "RE"), n_contrasts = .N), by = feature_id]
  # site localisation confidence from the QC feature metadata
  data(list = paste0(t, "_PROT_PH_QC"), package = "MotrpacHumanPreSuspensionData")
  fm <- as.data.table(get(paste0(t, "_PROT_PH_QC"))$feature_metadata)[, .(feature_id = id, confident_site)]
  # attach the confidence flag, then keep only features that have a gene (inner join with the gene table)
  s <- fm[s, on = "feature_id"][f2g, on = "feature_id", nomatch = 0]
  # label the tissue ("muscle" or "adipose")
  s[, tissue := names(which(c(muscle = "MUSCLE", adipose = "ADIPOSE") == t))]
  s
}))
# single-site features: accession, position, residue -> is it a site GlyGen knows?
# Split each feature ID into its protein accession and its site string (e.g. "S758s").
ph[, `:=`(acc = sub("_.*$", "", feature_id), sites = sub("^[^_]*_", "", feature_id))]
# Single-site feature = exactly one residue + position + lower-case residue.
ph[, single := grepl("^[STY][0-9]+[sty]$", sites)]
# Build the same match key as for GlyGen (only for single-site features).
ph[single == TRUE, key := paste(acc, sub("^[STY]([0-9]+).*$", "\\1", sites), substr(sites, 1, 1))]
# Known in GlyGen = the key is in GlyGen's list.
ph[, known_in_glygen := !is.na(key) & key %in% gps$key]
# one row per protein: sites measured and responding, per tissue
# Keep our 471 genes; per gene and tissue count sites, responding sites (EE, RE, either), confident sites, single
# sites and GlyGen-known sites; then one column per count x tissue (0 where the tissue has none).
php <- dcast(ph[entrez_gene %in% genes$entrez_gene, .(sites = .N, sig_EE = sum(sig_EE), sig_RE = sum(sig_RE), sig_any = sum(sig_EE | sig_RE),
                                                        confident = sum(confident_site, na.rm = TRUE), single_site = sum(single), known_in_glygen = sum(known_in_glygen)),
                 by = .(entrez_gene, tissue)], entrez_gene ~ tissue, value.var = c("sites", "sig_EE", "sig_RE", "sig_any", "confident", "single_site", "known_in_glygen"), fill = 0)
# Prefix the count columns with "motrpac_ph_".
setnames(php, names(php)[-1], paste0("motrpac_ph_", names(php)[-1]))

# ---- GlyGen per protein -------------------------------------------------------------------------------------
# GlyGen per-protein counts (from glygen_protein_inventory.py).
gg <- fread(file.path(INV, "glygen_protein_counts.csv"), colClasses = list(character = "entrez_gene"))
# One row per our 471 genes: GlyGen counts and MoTrPAC counts attached (left joins).
P <- php[gg[, !c("gene_symbol", "uniprot")][genes, on = "entrez_gene"], on = "entrez_gene"]
# Genes with no MoTrPAC phosphosite get 0 instead of NA.
for (v in grep("^motrpac_ph_", names(P), value = TRUE)) set(P, which(is.na(P[[v]])), v, 0L)
# Put the ID columns first, then save.
setcolorder(P, c("entrez_gene", "gene_symbol", "uniprot", "glygen_ac", "mapped_via"))
fwrite(P, file.path(INV, "protein_inventory.csv"))

# ---- metabolites: GlyGen glycans and GlyGen (Rhea) reaction participants --------------------------------------
# GlyGen files are downloaded once into inventory/glygen_files and reused, so reruns are reproducible
gget <- function(f) { d <- file.path(INV, "glygen_files"); dir.create(d, showWarnings = FALSE); p <- file.path(d, f)
  if (!file.exists(p)) download.file(paste0(GG, f), p, quiet = TRUE); p }
# GlyGen glycan cross-references: to ChEBI, to PubChem (compound entries only), and InChIKeys.
xc <- fread(gget("glycan_xref_chebi.csv"), colClasses = "character")
xp <- fread(gget("glycan_xref_pubchem.csv"), colClasses = "character")[grepl("compound", xref_key)]
xi <- fread(gget("glycan_sequences_inchi.csv"), colClasses = "character", select = c("glytoucan_ac", "inchi_key"))
# our metabolites' ChEBI IDs (all candidates, numbers only), PubChem CID and InChIKey
mets[, chebi_nums := lapply(strsplit(fifelse(is.na(chebi_all) | chebi_all == "", fifelse(is.na(chebi_id), "", chebi_id), chebi_all), "[;|, ]+"),
                            function(v) unique(sub("^CHEBI:", "", v[v != ""])))]
# Helper: for each metabolite's key(s), the GlyGen glycans (GlyTouCan IDs) whose column matches.
hit <- function(keys, table, col) lapply(keys, function(k) unique(table$glytoucan_ac[table[[col]] %in% k]))
# Glycan matches by ChEBI, by PubChem CID, and by the first InChIKey block (14 characters = the molecular skeleton,
# ignoring stereochemistry); missing IDs become "".
mets[, `:=`(glycan_by_chebi = hit(chebi_nums, xc, "xref_id"),
            glycan_by_pubchem = hit(as.list(fifelse(is.na(pubchem_cid), "", pubchem_cid)), xp, "xref_id"),
            glycan_by_inchikey = hit(as.list(fifelse(is.na(inchi_key), "", substr(inchi_key, 1, 14))), xi[, .(glytoucan_ac, sk = substr(inchi_key, 1, 14))], "sk"))]
# All matched glycans per metabolite, joined with ";" ("" if none).
mets[, glygen_glycans := mapply(function(a, b, c) paste(sort(unique(c(a, b, c))), collapse = ";"), glycan_by_chebi, glycan_by_pubchem, glycan_by_inchikey)]
# Rhea enzyme links from our own step 5 (for comparison; GlyGen's Rhea file has no small molecules)
# Our step-5 metabolite-enzyme links.
lk <- fread(file.path(OUT, "05_metabolite_protein_links.csv"))
# Number of distinct enzymes (genes) linked to each metabolite.
mets[, rhea_enzymes_among_471 := sapply(metabolite, function(m) uniqueN(lk[metabolite == m, gene_symbol]))]
# One row per metabolite with the match flags; save.
M <- mets[, .(metabolite, super_class, chebi_id, pubchem_cid, inchi_key, has_chebi = lengths(chebi_nums) > 0,
              glygen_glycan = glygen_glycans != "", glygen_glycans,
              via_chebi = lengths(glycan_by_chebi) > 0, via_pubchem = lengths(glycan_by_pubchem) > 0, via_inchikey = lengths(glycan_by_inchikey) > 0,
              rhea_enzymes_among_471)]
fwrite(M, file.path(INV, "metabolite_inventory.csv"))

# ---- coverage summary ----------------------------------------------------------------------------------------
# Helper: one coverage row = how many entities have a value above 0, out of how many, and the median among those.
cov <- function(entity, source, type, v, what) data.table(entity = entity, source = source, data_type = type, what = what,
                                                         n_with_data = sum(v > 0, na.rm = TRUE), of = length(v),
                                                         median_if_present = if (any(v > 0, na.rm = TRUE)) as.numeric(median(v[v > 0], na.rm = TRUE)) else NA_real_)
# Shortcut for a protein row taken from one column of the protein table.
pr <- function(src, type, col, what) cov("protein", src, type, P[[col]], what)
# The coverage table: one row per data type (proteins first, then metabolites).
S <- rbind(
  cov("protein", "MoTrPAC", "phosphoproteomics (either tissue)", P$motrpac_ph_sites_muscle + P$motrpac_ph_sites_adipose, "phosphosites measured"),
  pr("MoTrPAC", "phosphoproteomics, muscle (0.5/4/24 h)", "motrpac_ph_sites_muscle", "phosphosites measured"),
  pr("MoTrPAC", "phosphoproteomics, adipose (4 h)", "motrpac_ph_sites_adipose", "phosphosites measured"),
  cov("protein", "MoTrPAC", "phosphosite responds (adj. p < 0.05, EE or RE vs control)", P$motrpac_ph_sig_any_muscle + P$motrpac_ph_sig_any_adipose, "responding sites"),
  cov("protein", "MoTrPAC", "measured site also known in GlyGen", P$motrpac_ph_known_in_glygen_muscle + P$motrpac_ph_known_in_glygen_adipose, "matching sites"),
  pr("GlyGen", "phosphorylation sites (UniProtKB, iPTMnet)", "phospho_sites_unique", "sites"),
  pr("GlyGen", "phosphorylation sites with a known kinase", "phospho_sites_with_kinase", "sites"),
  pr("GlyGen", "glycosylation: any site (position known)", "gly_sites_unique", "sites"),
  pr("GlyGen", "glycosylated, site unknown (protein-level evidence)", "gly_protein_level_no_site", "yes"),
  cov("protein", "GlyGen", "glycosylated at all (site known or protein-level)", P$gly_sites_unique + P$gly_protein_level_no_site, "sites or yes"),
  pr("GlyGen", "O-GlcNAc sites (position known)", "gly_oglcnac_sites", "sites"),
  pr("GlyGen", "glycosylation: reported sites", "gly_sites_reported", "sites"),
  pr("GlyGen", "glycosylation: reported sites with glycan structure", "gly_sites_reported_with_glycan", "sites"),
  pr("GlyGen", "glycosylation: predicted sites", "gly_sites_predicted", "sites"),
  pr("GlyGen", "glycosylation: literature-mined sites", "gly_sites_automatic_literature_mining", "sites"),
  pr("GlyGen", "glycosylation: N-linked sites", "gly_sites_N", "sites"),
  pr("GlyGen", "glycosylation: O-linked sites", "gly_sites_O", "sites"),
  pr("GlyGen", "glycan structures observed at a site (GlyTouCan)", "glycans_at_sites", "glycans"),
  pr("GlyGen", "glycation sites", "glycation_sites_unique", "sites"),
  pr("GlyGen", "PTM annotation (UniProtKB)", "n_ptm_annotation", "records"),
  pr("GlyGen", "site annotation (active / binding sites)", "n_site_annotation", "records"),
  pr("GlyGen", "mutations / SNVs", "n_snv", "records"),
  pr("GlyGen", "mutagenesis", "n_mutagenesis", "records"),
  pr("GlyGen", "disease associations", "n_disease", "records"),
  pr("GlyGen", "biomarkers", "n_biomarkers", "records"),
  pr("GlyGen", "expression, normal tissue (Bgee)", "n_expression_tissue", "tissues"),
  pr("GlyGen", "expression, cancer (BioXpress)", "n_expression_disease", "records"),
  pr("GlyGen", "pathways (Reactome, KEGG)", "n_pathway", "pathways"),
  pr("GlyGen", "reactions (Reactome, Rhea)", "n_reactions", "reactions"),
  pr("GlyGen", "reaction participants (other molecules in its reactions)", "n_reaction_participants", "participants"),
  pr("GlyGen", "enzyme annotation (EC)", "n_enzyme_annotation", "records"),
  pr("GlyGen", "GO annotation", "n_go_terms", "terms"),
  pr("GlyGen", "function text", "n_function", "records"),
  pr("GlyGen", "3D structures (PDB)", "n_structures", "structures"),
  pr("GlyGen", "isoforms", "n_isoforms", "isoforms"),
  pr("GlyGen", "orthologs", "n_orthologs", "orthologs"),
  pr("GlyGen", "publications", "n_publication", "publications"),
  pr("GlyGen", "cross-references to other databases", "n_crossref", "links"),
  cov("metabolite", "GlyGen", "is a GlyGen glycan (ChEBI / PubChem / InChIKey)", as.integer(M$glygen_glycan), "yes"),
  cov("metabolite", "Rhea (our step 5)", "enzyme among our 471 proteins (for comparison)", M$rhea_enzymes_among_471, "enzymes"),
  cov("metabolite", "(ours)", "has a ChEBI ID (needed for GlyGen matching)", as.integer(M$has_chebi), "yes"))
# Save the coverage table and show it.
fwrite(S, file.path(INV, "coverage_summary.csv"))
print(S[, .(entity, source, data_type, n_with_data, of, median_if_present)], nrows = 100)
