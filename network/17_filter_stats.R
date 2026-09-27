#!/usr/bin/env Rscript
# =====================================================================================================
# 17_filter_stats.R — STEP 17 (PART 1): STATISTICS BEHIND THE FILTERABLE INTERACTIVE NETWORKS
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   The interactive pages (17_interactive_networks.R) let the user filter the networks by ome, tissue, time
#   point and arm. For every node and every such "cell" the page needs the change and its significance, and
#   for groups of connected nodes ("modules") it needs a test of whether the group as a whole responds. This
#   script precomputes both, with the consortium's own statistics and gene-set method.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Node statistics: for each gene (the exact transcript / protein used in step 1) and metabolite (the
#      platform used in step 1b), the log fold change and adjusted p-value in every tissue x ome x time point
#      for three contrasts (delta-delta, from the MoTrPAC differential analysis): endurance vs control (EE),
#      resistance vs control (RE), and endurance vs resistance (EE-RE).
#   2. Modules: communities found from each network's STRUCTURE alone (Louvain on the unweighted edges, fixed
#      seed; the exercise data are not used), for the joint, gene and metabolite networks; modules with fewer
#      than 5 members are not tested (MoTrPAC's minimum set size).
#   3. Module tests: MoTrPAC's run_cameraPR() (limma CAMERA-PR via TMSig; the method behind the package's
#      published pathway results) with the modules as custom gene sets (gene symbols + metabolite names), on
#      every tissue x ome x contrast. Each module is split into its genes (tested in RNA and protein) and its
#      metabolites (tested in metabolomics); a part is tested if at least 5 members, and at least 70% of them
#      (MoTrPAC's default), are measured in that ome; competitive test against all features of that ome and
#      tissue; FDR across the modules of one network within each cell.
#   5. Module names: ORA (hypergeometric, MoTrPAC run_ORA) of each module's genes against our 471 genes and of its
#      metabolites against our 450 metabolites; a module is named after its most significant pathway (adj. p < 0.05;
#      Reactome, KEGG, WikiPathways, PID, BioCarta, GO BP, MitoCarta) and, if significant, its RefMet class.
#   4. Annotation layers: MoTrPAC phosphosites of our proteins (logFC / adj. p per tissue x arm x time, with
#      GlyGen knowledge per site: known site, kinase, O-GlcNAc crosstalk) and GlyGen per-protein annotation counts (needs network/inventory and step 16).
#
# HOW TO RUN
#   After steps 1, 1b, 3, 6, 10, 14, 15:   Rscript network/17_filter_stats.R   (a few minutes; loads all
#   differential-analysis tables). Then:   Rscript network/17_interactive_networks.R
#
# DATA AND PROVENANCE
#   MotrpacHumanPreSuspensionAnalysis v0.2.4: *_TRNSCRPT_DA, *_PROT_PR_DA, BLOOD_PROT_OL_DA, *_METAB_DA,
#   run_cameraPR(), CONTRAST_CONVERTER. Feature choices from steps 1 / 1b (provenance tables).
#
# TECH STACK:  R 4.4; data.table, igraph (Louvain), MotrpacHumanPreSuspensionAnalysis + TMSig (CAMERA-PR).
#
# INPUTS ($HACK_OUT)
#   01_nodes_feature_provenance.csv, 01b_metab_nodes_provenance.csv, 03_weighted_edges.csv,
#   06_metabolite_edges.csv, 14_joint_edges.csv
#   02_nodes_string.csv (step 2: gene symbols and UniProt accessions),
#   inventory/protein_inventory.csv and inventory/glygen_glycosites.csv (network/inventory)
# INPUTS (MNET_DIR, default ~/Desktop/output/hackathon/resources/mo_annotation)
#   phosphosites.csv, glycosites.csv, motrpac_feature_site_map.csv, proteins_ptm.csv (the team's mnet resource)
# OUTPUTS ($HACK_OUT)
#   17_node_cell_stats.csv   node, node_type, tissue, ome (rna / prot / metab), arm (EE / RE / ER), time, logFC, adj_p
#   17_modules.csv           network (joint / gene / metabolite), module, node, node_type
#   17_module_modules.gmt    the modules as gene sets (input to run_cameraPR)
#   17_module_camera.csv     network, module, tissue, ome, arm, time, n_members_tested, direction, z, p, fdr
#   17_phospho_site_stats.csv   protein, feature_id, site, tissue, arm, time, logFC, adj_p, known_in_glygen, kinases, crosstalk
#   17_module_ora.csv        module x gene set: ORA against our universe (471 genes / 450 metabolites), overlap
#                            members, p, BH adj. p (MoTrPAC run_ORA; pathway, GO, MitoCarta, CellMarker, RefMet sets)
#   17_module_names.csv      one readable name per module (most significant pathway / metabolite class, or hubs)
#   17_glygen_protein_annotation.csv   one row per protein: GlyGen counts (phosphosites, glycosylation, mutations,
#                                disease, pathways, ...) and crosstalk residues
#   17_crosstalk_sites.csv   protein, residue, site_id, features, responds (EE / RE / both / no), glycosylation
#                            type and source, for residues that are both a MoTrPAC phosphosite and an O-glycosylation site
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   Every drawn node has statistics for each of its measured cells; the script stops if a node's chosen feature
#   is missing from its differential-analysis table, or if the logFC read here differs from the step 1 / 1b
#   raw logFC for the EE and RE contrasts.
#
# KNOWN LIMITS
#   Modules are structural communities, one reasonable definition among several (data-driven "active modules"
#   are a possible extension); CAMERA-PR asks whether a module moves MORE than other features of the same ome
#   (competitive), not whether it moves at all. Multiple testing is corrected across modules within each cell,
#   not across cells.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph) })

# Folders (outside the repo).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Fixed seed for the community detection.
SEED <- 20260926
# Minimum members per module and ome for a test (MoTrPAC's min_size).
MIN_SIZE <- 5L
# Constants shared with steps 1 / 1b.
TISSUES <- c("adipose", "blood", "muscle")
TPS <- c("0.5h" = "post_15_30_45_min", "4h" = "post_3.5_4_hr", "24h" = "post_24_hr")
ARM_OF <- c("EE-CON" = "EE", "RE-CON" = "RE", "EE-RE" = "ER")
PROT_OBJ <- c(adipose = "ADIPOSE_PROT_PR_DA", blood = "BLOOD_PROT_OL_DA", muscle = "MUSCLE_PROT_PR_DA")
PKG <- "MotrpacHumanPreSuspensionAnalysis"

# ---- 1. node statistics ---------------------------------------------------------------------------------
# Helper: one tissue's differential-analysis rows for the three contrasts at the three post-exercise times.
da_rows <- function(obj) {
  data(list = obj, package = PKG)
  # load the table from the package
  x <- as.data.table(get(obj))[contrast_category %in% names(ARM_OF) & Timepoint %in% TPS]
  # keep the three contrasts (EE-CON, RE-CON, EE-RE) at the three post-exercise time points
  x[, .(feature_id = as.character(feature_id), platform = if ("platform" %in% names(x)) as.character(platform) else NA_character_,
        # standard columns: feature, platform (metabolomics only), arm code, time label, logFC, adjusted p-value
        arm = ARM_OF[as.character(contrast_category)], time = names(TPS)[match(as.character(Timepoint), TPS)], logFC, adj_p = adj_p_value)]
}
# Genes: the feature step 1 chose for each tissue x ome.
prov <- fread(file.path(OUT, "01_nodes_feature_provenance.csv"), colClasses = list(character = "entrez_gene"))
# for every tissue and each gene ome (RNA, protein):
G <- rbindlist(lapply(TISSUES, function(t) rbindlist(lapply(c("rna", "prot"), function(o) {
  # the results table for this tissue x ome (blood protein = OLINK)
  obj <- if (o == "rna") paste0(toupper(t), "_TRNSCRPT_DA") else PROT_OBJ[[t]]
  # the feature step 1 chose for each gene in this tissue x ome (genes without one are skipped)
  ch <- prov[, .(node = gene_symbol, feature_id = get(paste0("feature_id_", t, "_", o)))][!is.na(feature_id) & feature_id != ""]
  # that feature's rows from the results table
  x <- da_rows(obj)[ch, on = "feature_id", nomatch = 0]
  # safety check: every chosen feature is in its table
  stopifnot(all(ch$feature_id %in% x$feature_id))
  # label the rows with tissue and ome
  x[, `:=`(tissue = t, ome = o)]
}))))
# every gene row is a "protein" node (the network's name for gene / protein nodes)
G[, node_type := "protein"]
# Metabolites: the platform step 1b kept in each tissue.
mprov <- fread(file.path(OUT, "01b_metab_nodes_provenance.csv"))
M <- rbindlist(lapply(TISSUES, function(t) {
  # this tissue's metabolomics results table
  x <- da_rows(paste0(toupper(t), "_METAB_DA"))
  # the platform step 1b chose for each metabolite in this tissue (metabolites without one are skipped)
  ch <- mprov[, .(node = metabolite, platform = get(paste0("platform_", t)))][!is.na(platform) & platform != ""]
  # keep only the rows of the chosen platform (for metabolites the feature ID is the metabolite name)
  x <- x[ch, on = c(feature_id = "node", "platform"), nomatch = 0]
  # call the feature column "node"
  setnames(x, "feature_id", "node")
  # label the rows
  x[, `:=`(tissue = t, ome = "metab", node_type = "metabolite")]
}))
# All node statistics in one table.
S <- rbind(G[, .(node, node_type, tissue, ome, arm, time, logFC, adj_p)], M[, .(node, node_type, tissue, ome, arm, time, logFC, adj_p)])
# Safety check: the EE / RE logFCs equal the step 1 / 1b raw values used for the embeddings.
# Helper: read a step 1 / 1b raw-logFC table in long form (node, dimension name, value), empty cells dropped.
chk <- function(file, id) { r <- fread(file.path(OUT, file)); d <- setdiff(names(r), c("entrez_gene", "gene_symbol", "metabolite"))
  for (k in d) set(r, j = k, value = as.numeric(r[[k]]))   # always-empty columns are read as logical
  melt(r[, c(id, d), with = FALSE], id.vars = id, variable.name = "dim", value.name = "raw", na.rm = TRUE) }
# the step 1 / 1b ENDURANCE values of genes and metabolites, with one common node column
rawE <- rbind(chk("01_nodes_EE_raw_logFC.csv", "gene_symbol")[, node := gene_symbol][, gene_symbol := NULL],
              chk("01b_metab_nodes_EE_raw_logFC.csv", "metabolite")[, node := metabolite][, metabolite := NULL])
rawE[, dim := as.character(dim)]
# match this script's EE values to them by node and dimension name (e.g. "muscle_rna_4h")
cmp <- S[arm == "EE", .(node, dim = paste(tissue, ome, time, sep = "_"), logFC)][rawE, on = c("node", "dim"), nomatch = 0]
# stop unless some rows matched and all matched values agree (only the EE arm is compared here)
stopifnot(nrow(cmp) > 0, isTRUE(all.equal(cmp$logFC, cmp$raw)))
# save and report
fwrite(S, file.path(OUT, "17_node_cell_stats.csv"))
message(sprintf("node statistics: %d rows (%d nodes)", nrow(S), uniqueN(S$node)))

# ---- 2. modules: structural communities of each network ---------------------------------------------------
# Edge lists of the three networks (steps 3, 6, 14).
w3 <- fread(file.path(OUT, "03_weighted_edges.csv")); w6 <- fread(file.path(OUT, "06_metabolite_edges.csv")); je <- fread(file.path(OUT, "14_joint_edges.csv"))
# (only the two end nodes of each edge: the module search ignores the weights)
nets <- list(joint = je[, .(a = node_a, b = node_b)], gene = w3[, .(a = symbol_a, b = symbol_b)], metabolite = w6[, .(a = metabolite_a, b = metabolite_b)])
# names of all measured metabolites, to tell metabolite nodes from protein nodes
mets <- unique(M$node)
# Find the communities (modules) of each network.
MOD <- rbindlist(lapply(names(nets), function(n) {
  # build the unweighted graph
  g <- graph_from_data_frame(nets[[n]], directed = FALSE)
  # Louvain community detection, with a fixed seed so it gives the same answer every run
  set.seed(SEED); cl <- cluster_louvain(g)
  # each node's community
  m <- data.table(network = n, node = V(g)$name, community = membership(cl))
  # keep communities with at least MIN_SIZE members; name them by size (largest first)
  sz <- m[, .N, by = community][N >= MIN_SIZE][order(-N)]
  m <- m[community %in% sz$community][, module := sprintf("%s_M%02d", n, match(community, sz$community))][, community := NULL]
  # node type
  m[, node_type := fifelse(node %in% mets, "metabolite", "protein")]
}))
# save and report the number of modules per network
fwrite(MOD, file.path(OUT, "17_modules.csv"))
message("modules: ", paste(MOD[, .(k = uniqueN(module)), by = network][, sprintf("%s %d", network, k)], collapse = ", "))
# the modules as a GMT file: each module split into its gene part ("<module>|genes", tested in RNA / protein)
# and its metabolite part ("<module>|metab", tested in metabolomics), so each part is tested only where it
# is measured (run_cameraPR keeps a set only if >= 70% of its members are in that ome's background)
# where the GMT file goes
gmt <- file.path(OUT, "17_module_modules.gmt")
# the set name of each node: "<module>|genes" or "<module>|metab"
parts <- MOD[, .(set = paste0(module, fifelse(node_type == "protein", "|genes", "|metab")), node)]
# one GMT line per set: name, a description, then the members (tab-separated)
writeLines(parts[, paste(c(set[1], "network module", node), collapse = "\t"), by = set]$V1, gmt)

# ---- 3. module tests: MoTrPAC run_cameraPR with the modules as gene sets -----------------------------------
# run CAMERA-PR on every tissue for RNA, protein (MS and OLINK) and metabolomics; a set is tested only with
# >= 5 members and >= 70% of its members measured
cam <- as.data.table(MotrpacHumanPreSuspensionAnalysis::run_cameraPR(
  selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ol", "metab"), selected_tissues = "all",
  path_to_gmt = gmt, min_size = MIN_SIZE, overlap_cutoff = 0.7))
# keep the three contrast families used by the page and translate labels
cam <- cam[contrast_type %in% c("exercise_with_controls", "Endur_vs_Resist")]
# contrast labels as plain text
cam[, contrast_short := as.character(contrast_short)]
# arm code from the contrast label: EE, RE or ER (endurance vs resistance)
cam[, arm := fcase(grepl("^Endur.*Control", contrast_short), "EE", grepl("^Resist.*Control", contrast_short), "RE",
                   grepl("^Endur.*Resist", contrast_short), "ER", default = NA_character_)]
# time label from the contrast label (the first of the three time points it mentions)
cam[, time := names(TPS)[sapply(contrast_short, function(s) which(sapply(TPS, grepl, x = s))[1])]]
# ome code from the assay name (both protein platforms count as "prot")
cam[, ome := c("transcript-rna-seq" = "rna", "prot-pr" = "prot", "prot-ol" = "prot", "metab" = "metab")[as.character(assay)]]
# drop contrasts that are not one of the three arms at one of the three times
cam <- cam[!is.na(arm) & !is.na(time)]
# FDR across the modules of one network within each cell (the package adjusts per collection; redone here per network)
# split each set name back into its module and its part (genes / metab)
cam[, `:=`(module = sub("\\|.*$", "", as.character(set)), part = sub("^.*\\|", "", as.character(set)))]
# a gene part is only meaningful in RNA / protein, a metabolite part only in metabolomics
cam <- cam[(part == "genes" & ome %in% c("rna", "prot")) | (part == "metab" & ome == "metab")]
# network name from the module name (e.g. "joint_M01" -> "joint")
cam[, network := sub("_M[0-9]+$", "", module)]
# BH FDR per network within each tissue x ome x arm x time cell
cam[, fdr := p.adjust(p_value, "BH"), by = .(network, tissue, ome, arm, time)]
# the output columns
CAM <- cam[, .(network, module, tissue = as.character(tissue), ome, arm, time, n_members_tested = set_size,
               direction = as.character(direction), z = z.std, p = p_value, fdr)]
# save and report
fwrite(CAM, file.path(OUT, "17_module_camera.csv"))
message(sprintf("module tests: %d (modules x cells); FDR < 0.05: %d", nrow(CAM), sum(CAM$fdr < 0.05)))

# ---- 4. annotation layers: MoTrPAC phospho + mnet PTM (UniProt / OmniPath) + GlyGen extras ------------------
# PTM knowledge comes from the team's mnet resource (resources/mo_annotation): phosphosites (UniProt + OmniPath)
# with kinases, UniProt glycosites, and its isoform-aware bridge from MoTrPAC
# feature IDs to canonical sites. GlyGen (network/inventory) adds only what mnet lacks: glycan structures,
# protein-level O-GlcNAc evidence, O-glycosylation sites from other databases, and non-PTM fields.
# Where the mnet annotation files are (override with MNET_DIR; the default is outside the repo).
MNET_DIR <- Sys.getenv("MNET_DIR", unset = path.expand("~/Desktop/output/hackathon/resources/mo_annotation"))
# The GlyGen inventory folder; stop with a message if it has not been built.
INV <- file.path(OUT, "inventory")
if (!file.exists(file.path(INV, "protein_inventory.csv"))) stop("run network/inventory/ first (protein_inventory.csv missing)")
# Our genes with symbols and UniProt accessions (step 2).
genes <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol, uniprot)]
# our genes by UniProt accession (a gene can carry several accessions, ";"-separated)
g_acc <- genes[, .(acc = unlist(strsplit(uniprot, ";"))), by = .(entrez_gene, gene_symbol)]
# mnet tables: phosphosites (with kinases), glycosites, the MoTrPAC feature -> site bridge, per-protein PTM counts
mn_ph <- fread(file.path(MNET_DIR, "phosphosites.csv"))
mn_gl <- fread(file.path(MNET_DIR, "glycosites.csv"))
mn_map <- fread(file.path(MNET_DIR, "motrpac_feature_site_map.csv"))
mn_ptm <- fread(file.path(MNET_DIR, "proteins_ptm.csv"))
# 4a. MoTrPAC phosphosites of our proteins: logFC and adj. p per tissue x arm x time (muscle 0.5/4/24 h, adipose 4 h)
# MoTrPAC's feature -> gene table, phosphosite features only, limited to our genes
data("HUMAN_FEATURE_TO_GENE", package = PKG)
f2g <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay == "prot-ph", .(feature_id = as.character(feature_id), entrez_gene = as.character(entrez_gene))])
f2g <- genes[, .(entrez_gene, gene_symbol)][f2g, on = "entrez_gene", nomatch = 0]
# phospho results (EE, RE, EE-RE at 0.5 / 4 / 24 h where measured) in muscle and adipose
PH <- rbindlist(lapply(c("muscle", "adipose"), function(t) da_rows(paste0(toupper(t), "_PROT_PH_DA"))[, tissue := t]))
# attach the gene (features not on our genes are dropped); the platform column is not needed
PH <- f2g[PH, on = "feature_id", nomatch = 0][, platform := NULL]
# readable site from the feature ID (e.g. "O60664_S31s" -> "S31")
PH[, site := gsub("([sty])", "", sub("^[^_]*_", "", feature_id))]
# canonical site IDs per feature from the mnet bridge (all mapping statuses except residue mismatches)
br <- unique(mn_map[mapping_status != "residue_mismatch" & !is.na(site_id) & site_id != "", .(feature_id, site_id)])
# O-glycosylation sites (Ser / Thr / Tyr): mnet UniProt glycosites + GlyGen sites (other databases)
og_mnet <- mn_gl[residue %in% c("S", "T", "Y") & glyco_type %in% c("O-linked", "O-GlcNAc"), .(site_id, gly_type = glyco_type, gly_source = "UniProt (mnet)")]
# GlyGen O-linked sites on Ser / Thr / Tyr, rewritten into the same site-ID format
gg <- fread(file.path(INV, "glygen_glycosites.csv"))[type == "O-linked" & residue %in% c("Ser", "Thr", "Tyr")]
og_gg <- unique(gg[, .(site_id = paste0(sub("-.*$", "", glygen_ac), "_", substr(residue, 1, 1), position), gly_type = fifelse(subtype == "O-GlcNAcylation", "O-GlcNAc", "O-linked"),
                       gly_source = paste0("GlyGen: ", source))])
# one row per site: all glycosylation types and sources that report it
OG <- rbind(og_mnet, og_gg)[, .(gly_type = paste(sort(unique(gly_type)), collapse = ";"), gly_source = paste(sort(unique(gly_source)), collapse = "; ")), by = site_id]
# per feature: known in mnet phosphosites, kinases (mnet), crosstalk (any of its sites is an O-glycosylation site)
# the kinases known for each site (separators unified to ";")
kin_site <- mn_ph[kinases != "" & !is.na(kinases), .(site_id, kinases = gsub("[|,]", ";", kinases))]
fx <- br[, .(known = any(site_id %in% mn_ph$site_id), kinases = paste(sort(unique(unlist(strsplit(kin_site$kinases[match(site_id, kin_site$site_id)], ";")))), collapse = ";"),
             crosstalk = any(site_id %in% OG$site_id), sites = paste(site_id, collapse = ";")), by = feature_id]
# no kinase -> empty text
fx[kinases == "NA", kinases := ""]
# attach to every phospho row
PH <- fx[PH, on = "feature_id"]
# features without a mapped site: not known, no crosstalk, no kinases
PH[is.na(known), `:=`(known = FALSE, crosstalk = FALSE, kinases = "")]
# save (note: the column known_in_glygen holds "site listed in mnet phosphosites", i.e. UniProt / OmniPath)
fwrite(PH[, .(protein = gene_symbol, feature_id, site, tissue, arm, time, logFC, adj_p, known_in_glygen = known, kinases, crosstalk)],
       file.path(OUT, "17_phospho_site_stats.csv"))
# report
message(sprintf("phospho sites: %d features on %d proteins (%d mapped to mnet sites; %d crosstalk)", uniqueN(PH$feature_id), uniqueN(PH$gene_symbol),
                uniqueN(PH[feature_id %in% br$feature_id, feature_id]), uniqueN(PH[crosstalk == TRUE, feature_id])))
# crosstalk residues (replaces the GlyGen-only matching of step 16 for the pages)
# (same route as step 16: MoTrPAC feature -> our gene -> mnet canonical site -> O-glycosylation site)
# per feature: responds after EE / RE at any tissue and time (adj. p < 0.05); keep only features whose
# canonical site is an O-glycosylation site
xtab <- merge(PH[, .(responds_EE = any(adj_p < 0.05 & arm == "EE"), responds_RE = any(adj_p < 0.05 & arm == "RE")), by = .(feature_id, protein = gene_symbol)],
              br, by = "feature_id")[site_id %in% OG$site_id]
# one row per protein x residue
xtab <- xtab[, .(features = paste(sort(unique(feature_id)), collapse = ";"), responds_EE = any(responds_EE), responds_RE = any(responds_RE)), by = .(protein, site_id)]
# add the glycosylation type and source
xtab <- OG[xtab, on = "site_id"]
# residue label and one-word arm pattern
xtab[, `:=`(residue = sub("^[^_]*_", "", site_id), responds = fcase(responds_EE & responds_RE, "both", responds_EE, "endurance", responds_RE, "resistance", default = "no"))]
# save and report
fwrite(xtab[, .(protein, residue, site_id, features, responds, responds_EE, responds_RE, gly_type, gly_source)], file.path(OUT, "17_crosstalk_sites.csv"))
message(sprintf("crosstalk residues (MoTrPAC phosphosite = O-glycosylation site; mnet + GlyGen): %d on %d proteins; %d respond",
                nrow(xtab), uniqueN(xtab$protein), sum(xtab$responds != "no")))
# (no kinase -> substrate EDGES: network edges are only physical STRING / Rhea links weighted by the dot product)
# 4c. protein annotations: mnet PTM counts + GlyGen extras (one row per protein)
# mnet PTM counts per gene (if a gene has several accessions, the one with the most phosphosites is kept)
mp <- merge(g_acc, mn_ptm[, .(acc = sub("-[0-9]+$", "", node_id), n_phosphosites, n_phosphosites_with_kinase, is_kinase, n_substrate_sites,
                                n_glycosites, n_N_linked, n_O_linked, n_O_GlcNAc)], by = "acc")[order(-n_phosphosites)][!duplicated(gene_symbol)]
# the GlyGen inventory, with its columns renamed
P <- fread(file.path(INV, "protein_inventory.csv"))
ANN <- P[, .(protein = gene_symbol, glycans = glycans_at_sites, glyco_protein_level = gly_protein_level_no_site,
             mutations = n_snv, disease = n_disease, biomarkers = n_biomarkers, ptm_annotation = n_ptm_annotation,
             site_annotation = n_site_annotation, enzyme = n_enzyme_annotation, pathways = n_pathway, reactions = n_reactions,
             expression_tissues = n_expression_tissue, publications = n_publication)]
# combine the mnet counts and the GlyGen columns, one row per protein (proteins from either side kept)
ANN <- merge(mp[, .(protein = gene_symbol, phosphosites = n_phosphosites, kinase_sites = n_phosphosites_with_kinase, is_kinase = as.integer(is_kinase %in% c(TRUE, "True")),
                    substrate_sites = n_substrate_sites, glyco_sites = n_glycosites, glyco_N_sites = n_N_linked, glyco_O_sites = n_O_linked,
                    glyco_OGlcNAc_sites = n_O_GlcNAc)],
             ANN, by = "protein", all = TRUE)
# missing counts -> 0
for (k in setdiff(names(ANN), "protein")) set(ANN, which(is.na(ANN[[k]])), k, 0L)
# glycosylated = any known site or any protein-level evidence
ANN[, glycosylated := as.integer(glyco_sites + glyco_protein_level > 0)]
# add the number of crosstalk residues (0 if none)
ANN <- merge(ANN, xtab[, .(crosstalk_residues = .N), by = protein], by = "protein", all.x = TRUE)[is.na(crosstalk_residues), crosstalk_residues := 0L]
# save and report
fwrite(ANN, file.path(OUT, "17_glygen_protein_annotation.csv"))
message(sprintf("protein annotations: %d proteins x %d fields (PTM from mnet, extras from GlyGen)", nrow(ANN), ncol(ANN) - 1))

# ---- 5. module names: over-representation of pathways among each module's members ---------------------------
# Modules are gene / metabolite LISTS drawn from our 471-gene / 450-metabolite universe, so the right test is
# over-representation (hypergeometric ORA; GSEA needs a ranking), against THAT universe (not the genome), with
# MoTrPAC's run_ORA() and gene-set collections. MoTrPAC's 70% set-coverage rule is meant for genome-wide
# backgrounds and would remove almost every set here, so it is off (overlap_cutoff = 0); sets need >= 5 members
# in the universe (min_size) and >= 2 module members to name a module. BH within each collection and module.
# significance cut-off for naming a module
ALPHA_ORA <- 0.05
PATHWAY_DB <- c("REACTOME", "KEGG_MEDICUS", "WP", "PID", "BIOCARTA", "GOBP", "MITOCARTA")   # used for names
ORA_DB <- c(PATHWAY_DB, "GOCC", "GOMF", "CELLMARKER")                                     # also reported
# the two universes: our 471 genes and our measured metabolites
universe_genes <- genes$gene_symbol
universe_mets <- unique(M$node)
# Helper: ORA for one module part (needs >= 3 members in the universe); keeps sets hit by >= 2 members
ora_one <- function(members, background, db) {
  members <- intersect(members, background); if (length(members) < 3) return(NULL)
  r <- tryCatch(as.data.table(MotrpacHumanPreSuspensionAnalysis::run_ORA(input = members, background = background, database = db,
                                                                          min_size = MIN_SIZE, overlap_cutoff = 0)), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) return(NULL)
  r[set_size_in_input >= 2]
}
# for every module: ORA of its genes (gene collections) and of its metabolites (RefMet classes), labelled with
# its network and module
ORA <- rbindlist(lapply(split(MOD, MOD$module), function(m) {
  g <- ora_one(m[node_type == "protein", node], universe_genes, ORA_DB)
  k <- ora_one(m[node_type == "metabolite", node], universe_mets, "REFMET")
  r <- rbindlist(list(g, k), fill = TRUE); if (!nrow(r)) return(NULL)
  r[, `:=`(network = m$network[1], module = m$module[1])]
}), fill = TRUE)
# the module members found in each set (for tooltips / checking)
# every gene set / class by name, to list which module members fall in each
idx <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES
idx <- unlist(unname(idx[c(ORA_DB, "REFMET")]), recursive = FALSE)
ORA[, overlap := mapply(function(s, mo) paste(sort(intersect(idx[[s]], MOD[module == mo, node])), collapse = ";"), as.character(set), module)]
# readable set names: drop the collection prefix, sentence case, keep acronyms / gene symbols, add the source
# Readable collection names.
DB_LABEL <- c(REACTOME = "Reactome", KEGG_MEDICUS = "KEGG", WP = "WikiPathways", PID = "PID", BIOCARTA = "BioCarta", GOBP = "GO BP",
              GOCC = "GO CC", GOMF = "GO MF", MITOCARTA = "MitoCarta", CELLMARKER = "CellMarker", REFMET = "RefMet")
# Short words that are written in lower case when set names are made readable.
STOP <- c("of", "by", "to", "the", "in", "and", "via", "for", "on", "a", "an", "or", "with", "from", "into", "at", "as")
# Helper: e.g. "REACTOME_SIGNALING_BY_EGFR" -> "Signaling by EGFR"
pretty_set <- function(set, db) {
  # remove the collection prefix
  x <- sub("^(REACTOME|KEGG_MEDICUS|WP|PID|BIOCARTA|GOBP|GOCC|GOMF|MITOCARTA|CELLMARKER|REFMET)_", "", set)
  # RefMet class names are used as they are (only "Cer" is spelled out)
  if (db == "REFMET") return(if (x == "Cer") "Ceramides" else x)
  # split into words; keep words with digits, short vowel-less acronyms and listed acronyms in capitals
  w <- strsplit(x, "_")[[1]]
  keep <- grepl("[0-9]", w) | (nchar(w) <= 4 & !tolower(w) %in% STOP & !grepl("[AEIOU]", substr(w, 2, nchar(w)))) | w %in% c("ADME", "MHC", "NAD", "TNF", "ALPHA", "RNA", "DNA", "II", "III", "IV", "ER", "ATP", "GTP")
  # all other words in lower case
  w <- ifelse(keep, w, tolower(w)); w[w == "ALPHA"] <- "alpha"
  # join and capitalise the first letter
  out <- paste(w, collapse = " "); paste0(toupper(substr(out, 1, 1)), substring(out, 2))
}
# add the collection name in brackets
ORA[, set_label := mapply(pretty_set, as.character(set), as.character(database))]
ORA[, set_label := paste0(set_label, " [", DB_LABEL[as.character(database)], "]")]
# output columns, sorted by module and significance
ORA <- ORA[, .(network, module, database = as.character(database), set = as.character(set), set_label,
               set_size_in_universe = set_size, module_members_tested = input_size, overlap_n = set_size_in_input, overlap,
               p_value, adj_p_value)][order(module, adj_p_value, p_value)]
# save
fwrite(ORA, file.path(OUT, "17_module_ora.csv"))
# a readable name per module: its most significant pathway (Reactome, KEGG, WikiPathways, PID, BioCarta, GO BP,
# MitoCarta) if adj. p < 0.05, plus the metabolite class when a metabolite part is significant; otherwise
# "no significant pathway" with the module's best-connected members as a handle
# node degree (number of edges) in each network
deg <- rbindlist(lapply(names(nets), function(n) { g <- graph_from_data_frame(nets[[n]], directed = FALSE); data.table(network = n, node = V(g)$name, degree = degree(g)) }))
# the 3 best-connected members of each module (used when no pathway is significant)
hubs <- deg[MOD, on = c("network", "node")][order(-degree)][, .(hubs = paste(head(node, 3), collapse = ", ")), by = module]
# module sizes (all, proteins, metabolites)
nm <- MOD[, .(n = .N, n_prot = sum(node_type == "protein"), n_met = sum(node_type == "metabolite")), by = .(network, module)]
# most significant pathway per module (adj. p < 0.05)
best_path <- ORA[database %in% PATHWAY_DB & adj_p_value < ALPHA_ORA][order(adj_p_value, p_value)][, .SD[1], by = module]
# most significant RefMet class per module (adj. p < 0.05)
best_met <- ORA[database == "REFMET" & adj_p_value < ALPHA_ORA][order(adj_p_value, p_value)][, .SD[1], by = module]
# combine into one row per module
nm <- best_path[, .(module, path_name = set_label, path_db = database, path_fdr = adj_p_value, path_overlap = overlap_n)][nm, on = "module"]
nm <- best_met[, .(module, met_name = set_label, met_fdr = adj_p_value)][nm, on = "module"]
nm <- hubs[nm, on = "module"]
# number of significant pathway / class sets per module
nm[, n_sig_sets := sapply(module, function(mo) ORA[module == mo & database %in% c(PATHWAY_DB, "REFMET") & adj_p_value < ALPHA_ORA, .N])]
# the label: pathway + class, pathway only, class only, or "no significant pathway (hubs)"
nm[, label := fcase(!is.na(path_name) & !is.na(met_name), paste0(path_name, " + ", met_name),
                    !is.na(path_name), path_name, !is.na(met_name), met_name, default = paste0("no significant pathway (", hubs, ")"))]
# the name shown on the pages: short module ID + label
nm[, name := paste0(sub("^.*_", "", module), " · ", label)]
# column order, then save and report
setcolorder(nm, c("network", "module", "name", "label", "n", "n_prot", "n_met", "path_name", "path_db", "path_fdr", "path_overlap", "met_name", "met_fdr", "n_sig_sets", "hubs"))
fwrite(nm[order(network, module)], file.path(OUT, "17_module_names.csv"))
message(sprintf("module ORA: %d module x set rows; modules named by a significant pathway / class: %d of %d",
                nrow(ORA), sum(!grepl("^no significant", nm$label)), nrow(nm)))
# show the module names on screen
print(nm[order(network, module), .(module, name)])
