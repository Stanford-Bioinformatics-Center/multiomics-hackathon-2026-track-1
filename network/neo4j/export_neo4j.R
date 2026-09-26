#!/usr/bin/env Rscript
# =====================================================================================================
# network/neo4j/export_neo4j.R — WRITE THE EXERCISE NETWORKS AS NEO4J-READY CSV FILES
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   A teammate is building a Neo4j graph visualiser for our networks. This script gathers everything the
#   pipeline has produced about genes/proteins, metabolites, their exercise responses and their edges into
#   one tidy set of node and relationship CSV files that import.cypher loads into Neo4j. Re-run it whenever
#   the pipeline changes; new details are added here (one column = one new property).
#
# WHAT THIS SCRIPT DOES (plain language)
#   Reads the pipeline outputs (nothing is recomputed) and writes, to $NEO4J_IMPORT:
#     nodes_gene.csv, nodes_metabolite.csv, nodes_class.csv, nodes_contrast.csv, nodes_normalisation.csv,
#     nodes_dataset.csv, rel_string.csv, rel_class_link.csv, rel_rhea.csv, rel_in_class.csv,
#     rel_gene_response.csv, rel_metabolite_response.csv
#   The graph model (labels, relationship types, properties) is described in network/neo4j/README.md.
#   Lists (embeddings) are written as ";"-separated text, which import.cypher splits into number lists.
#
# HOW TO RUN
#   After the pipeline (steps 1-15; see network/README.md):   Rscript network/neo4j/export_neo4j.R
#   Then load into Neo4j with network/neo4j/import.cypher (see network/neo4j/README.md).
#
# DATA AND PROVENANCE
#   MoTrPAC human pre-suspension results (MotrpacHumanPreSuspensionAnalysis v0.2.4), normalised per ome
#   (steps 1 / 1b); STRING >= 700 (step 2); Rhea release 142 (step 5); RefMet / ChEBI (step 1c). The git
#   commit of the code that produced the export is recorded in nodes_dataset.csv.
#
# TECH STACK
#   R 4.4; data.table.
#
# INPUTS ($HACK_OUT; see network/neo4j/README.md, "Where each piece comes from")
#   01_nodes_{EE,RE}.csv, 01_nodes_{EE,RE}_raw_logFC.csv, 01_nodes_{EE,RE}_se.csv, 01_nodes_feature_provenance.csv,
#   01_scale_factors.csv, 01b_metab_nodes_{EE,RE}.csv (+ _raw_logFC, _se), 01b_metab_nodes_provenance.csv,
#   01b_metab_scale_factors.csv, 01c_metabolite_ids.csv, 02_edges.csv, 02_nodes_string.csv, 03_weighted_edges.csv,
#   05_metabolite_protein_links.csv, 06_metabolite_edges.csv, 10_layout_genes.csv, 10_layout_metabolites.csv,
#   14_joint_nodes.csv, 14_joint_edges.csv, 15_class_layout.csv
#
# OUTPUTS
#   $NEO4J_IMPORT (default $HACK_OUT/neo4j_import): the 12 CSV files above (never in the repo)
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   471 genes, 450 metabolites, 50 contrasts (per arm: 16 gene dimensions - RNA and protein in 3 tissues x
#   3 times, minus adipose protein at 0.5 h and 24 h - plus 9 metabolite dimensions), 431 STRING, 147 class-link and 186 Rhea relationships, 15,072 gene and 8,100 metabolite
#   response relationships. The script stops if any count differs from its source table, if an embedding has
#   a missing value outside the two always-empty adipose protein dimensions, or if a Rhea weight recomputed
#   from the exported embeddings differs from step 14.
#
# KNOWN LIMITS
#   Rhea reactions are stored as the example IDs kept by step 5, not as reaction nodes. Differences between
#   arms are descriptive (not tested). Layout coordinates are those of the static figures (0..1).
# =====================================================================================================

# Load packages quietly.
suppressMessages(library(data.table))

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where the Neo4j files go (override with NEO4J_IMPORT); outside the repo on purpose.
NEO <- Sys.getenv("NEO4J_IMPORT", unset = file.path(OUT, "neo4j_import"))
dir.create(NEO, recursive = TRUE, showWarnings = FALSE)
# Helper: read a pipeline table (Entrez IDs kept as text).
rd <- function(f) fread(file.path(OUT, f), colClasses = list(character = intersect(c("entrez_gene", "entrez_a", "entrez_b", "pubchem_cid"),
                                                                                   names(fread(file.path(OUT, f), nrows = 0)))))
# Helper: write one Neo4j file (empty cells for missing values; Neo4j reads them as absent properties).
wr <- function(x, f) { fwrite(x, file.path(NEO, f), na = ""); message(sprintf("-> %-30s %6d rows", f, nrow(x))) }
# Helper: join the numbers of each row into one ";"-separated text value (for list properties).
pack <- function(m) apply(m, 1, function(r) paste(format(r, digits = 15, trim = TRUE), collapse = ";"))
# Arm labels.
ARM <- c(EE = "endurance vs control", RE = "resistance vs control")

# ---- gene / protein nodes ---------------------------------------------------------------------------------
# Normalised embeddings per arm (18 columns: tissue x ome x time; adipose protein 0.5 h / 24 h always empty).
gE <- rd("01_nodes_EE.csv"); gR <- rd("01_nodes_RE.csv")
dims18 <- setdiff(names(gE), c("entrez_gene", "gene_symbol"))
# The 16 observed dimensions (columns that are not empty for every gene).
dims16 <- dims18[colSums(is.na(gE[, dims18, with = FALSE])) < nrow(gE)]
# Safety check: no missing value inside the 16 observed dimensions, in either arm.
stopifnot(length(dims16) == 16, !anyNA(gE[, dims16, with = FALSE]), !anyNA(gR[, dims16, with = FALSE]))
# STRING membership and the feature chosen for each tissue x ome (provenance).
sn <- rd("02_nodes_string.csv"); fp <- rd("01_nodes_feature_provenance.csv")
# Layouts: figure 10a (gene network) and figures 15a / 15b (joint network, class-grouped) and 14a / 14b.
lg <- rd("10_layout_genes.csv"); lj <- rd("15_class_layout.csv"); jn <- rd("14_joint_nodes.csv")
# One row per gene.
G <- gE[, .(entrez_gene, symbol = gene_symbol)]
G <- sn[, .(entrez_gene, uniprot, in_string, string_degree = degree, string_component = component, string_component_size = component_size)][G, on = "entrez_gene"]
G <- fp[, !"gene_symbol"][G, on = "entrez_gene"]
# Mean normalised response per arm (as the figure node colours).
G[, `:=`(mean_response_EE = rowMeans(as.matrix(gE[, dims16, with = FALSE])), mean_response_RE = rowMeans(as.matrix(gR[, dims16, with = FALSE])))]
# Joint-network strength and degree (step 14), and layout coordinates (empty if the gene has no edge there).
G <- jn[node_type == "protein", .(symbol = node, joint_strength_EE = strength_EE, joint_strength_RE = strength_RE, joint_degree = degree,
                                  x_joint14 = x, y_joint14 = y)][G, on = "symbol"]
G <- lj[node_type == "protein", .(symbol = node, x_joint = x, y_joint = y)][G, on = "symbol"]
G <- lg[, .(symbol = node, x_gene = x, y_gene = y)][G, on = "symbol"]
# Embeddings as number lists (16 observed dimensions, names in embedding_dims).
G[, `:=`(embedding_EE = pack(as.matrix(gE[match(G$entrez_gene, gE$entrez_gene), dims16, with = FALSE])),
         embedding_RE = pack(as.matrix(gR[match(G$entrez_gene, gR$entrez_gene), dims16, with = FALSE])),
         embedding_dims = paste(dims16, collapse = ";"))]
setcolorder(G, c("entrez_gene", "symbol", "uniprot"))
stopifnot(nrow(G) == 471, !anyDuplicated(G$entrez_gene), !anyDuplicated(G$symbol))
wr(G, "nodes_gene.csv")

# ---- metabolite nodes ---------------------------------------------------------------------------------------
mE <- rd("01b_metab_nodes_EE.csv"); mR <- rd("01b_metab_nodes_RE.csv")
dims9 <- setdiff(names(mE), "metabolite")
stopifnot(length(dims9) == 9, !anyNA(mE[, dims9, with = FALSE]), !anyNA(mR[, dims9, with = FALSE]))
ids <- rd("01c_metabolite_ids.csv"); pl <- rd("01b_metab_nodes_provenance.csv"); lm <- rd("10_layout_metabolites.csv")
M <- ids[, .(name = metabolite, refmet_name, refmet_id, super_class, main_class, pubchem_cid, inchi_key, chebi_id, chebi_all, lookup_status)]
M <- pl[, .(name = metabolite, platform_adipose, platform_blood, platform_muscle)][M, on = "name"]
M[, `:=`(mean_response_EE = rowMeans(as.matrix(mE[match(name, mE$metabolite), dims9, with = FALSE])),
         mean_response_RE = rowMeans(as.matrix(mR[match(name, mR$metabolite), dims9, with = FALSE])))]
M <- jn[node_type == "metabolite", .(name = node, joint_strength_EE = strength_EE, joint_strength_RE = strength_RE, joint_degree = degree,
                                     x_joint14 = x, y_joint14 = y)][M, on = "name"]
M <- lj[node_type == "metabolite", .(name = node, x_joint = x, y_joint = y)][M, on = "name"]
M <- lm[, .(name = node, x_metab = x, y_metab = y)][M, on = "name"]
# The metabolite embedding (9 dimensions) and the DOUBLED embedding aligned to the gene's 16 dimensions
# (each tissue x time value placed in both the RNA and the protein slot), so a metabolite - protein weight is
# one plain dot product: reduce(s = 0, i IN range(0, 15) | s + m.embedding_doubled_EE[i] * g.embedding_EE[i]).
dbl <- match(sub("_(rna|prot)_", "_metab_", dims16), dims9)
mEm <- as.matrix(mE[match(M$name, mE$metabolite), dims9, with = FALSE]); mRm <- as.matrix(mR[match(M$name, mR$metabolite), dims9, with = FALSE])
M[, `:=`(embedding_EE = pack(mEm), embedding_RE = pack(mRm), embedding_dims = paste(dims9, collapse = ";"),
         embedding_doubled_EE = pack(mEm[, dbl]), embedding_doubled_RE = pack(mRm[, dbl]), embedding_doubled_dims = paste(dims16, collapse = ";"))]
setcolorder(M, c("name", "refmet_name", "refmet_id", "super_class", "main_class"))
stopifnot(nrow(M) == 450, !anyDuplicated(M$name), !any(M$name %in% G$symbol))
wr(M, "nodes_metabolite.csv")

# ---- metabolite classes, contrasts, normalisation, dataset -------------------------------------------------
# One node per RefMet super class, with how many of our metabolites it holds.
C <- M[, .(n_metabolites = .N, n_in_joint_network = sum(!is.na(x_joint))), by = .(name = super_class)][order(-n_metabolites)]
wr(C[!is.na(name) & name != ""], "nodes_class.csv")
# One node per contrast (arm x tissue x ome x time); id e.g. "EE|blood|rna|4h".
ct <- function(d, arm) { p <- tstrsplit(d, "_"); data.table(contrast_id = paste(arm, p[[1]], p[[2]], p[[3]], sep = "|"), arm = arm, arm_label = ARM[[arm]],
                                                            tissue = p[[1]], ome = p[[2]], time = p[[3]], time_h = as.numeric(sub("h$", "", p[[3]]))) }
K <- unique(rbindlist(lapply(c("EE", "RE"), function(a) rbind(ct(dims16, a), ct(dims9, a)))))
stopifnot(nrow(K) == 50)
wr(K, "nodes_contrast.csv")
# The normalisation divisors (the critical QC step: each ome divided by its maximum |logFC|).
sf <- rbind(rd("01_scale_factors.csv"), rd("01b_metab_scale_factors.csv"))
wr(sf[, .(ome, divisor_max_abs_logFC = max_abs_logFC, set_by, n_values)], "nodes_normalisation.csv")
# Export provenance: date and the git commit of the code (run from inside the repo; empty otherwise).
commit <- suppressWarnings(tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) character(0)))
wr(data.table(name = "MoTrPAC human pre-suspension exercise networks", exported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
              code_commit = if (length(commit)) commit[1] else NA_character_, source_package = "MotrpacHumanPreSuspensionAnalysis v0.2.4",
              string_cutoff = 700, rhea_release = "142 (2026-09-02)", normalisation = "per-ome max |logFC| over tissues, times and both arms"),
   "nodes_dataset.csv")

# ---- relationships ------------------------------------------------------------------------------------------
# Gene - gene: STRING evidence (step 2) with per-arm weights (step 3).
e2 <- rd("02_edges.csv"); w3 <- rd("03_weighted_edges.csv")
S <- w3[e2[, .(entrez_a, entrez_b, uniprot_a, uniprot_b, cos_EE, cos_RE)], on = c("entrez_a", "entrez_b")]
stopifnot(nrow(S) == nrow(w3), !anyNA(S$w_EE))
wr(S[, .(entrez_a, entrez_b, symbol_a, symbol_b, combined_score, w_EE, w_RE, w_diff, sig_EE, sig_RE, cos_EE, cos_RE)], "rel_string.csv")
# Metabolite - metabolite: the class rule (step 6).
w6 <- rd("06_metabolite_edges.csv")
wr(w6[, .(metabolite_a, metabolite_b, class, class_level, link_type, n_shared_proteins, shared_proteins, string_protein_pairs,
          w_EE, w_RE, w_diff, sig_EE, sig_RE)], "rel_class_link.csv")
# Metabolite - protein: Rhea (step 5) with the doubled-embedding weights (step 14).
lk <- unique(rd("05_metabolite_protein_links.csv"), by = c("metabolite", "entrez_gene"))
je <- rd("14_joint_edges.csv")[edge_type == "metabolite - protein"]
R <- je[, .(metabolite = node_a, gene_symbol = node_b, w_EE, w_RE, w_diff)][lk, on = c("metabolite", "gene_symbol")]
stopifnot(nrow(R) == 186, !anyNA(R$w_EE))
# Safety check: the weights equal the dot product of the exported doubled metabolite and gene embeddings.
dot <- rowSums(mEm[match(R$metabolite, M$name), dbl] * as.matrix(gE[match(R$entrez_gene, gE$entrez_gene), dims16, with = FALSE]))
stopifnot(isTRUE(all.equal(unname(dot), R$w_EE)))
wr(R[, .(metabolite, entrez_gene, gene_symbol, uniprot, n_reactions, example_reactions, matched_via, w_EE, w_RE, w_diff)], "rel_rhea.csv")
# Metabolite - class membership.
wr(M[!is.na(super_class) & super_class != "", .(metabolite = name, super_class)], "rel_in_class.csv")
# Responses: one relationship per node x contrast with the raw logFC, the normalised value and its SE.
long <- function(id, norm, raw, se, arm, dims) {
  n <- melt(norm[, c(id, dims), with = FALSE], id.vars = id, variable.name = "dim", value.name = "normalised")
  r <- melt(raw[, c(id, dims), with = FALSE], id.vars = id, variable.name = "dim", value.name = "logFC")
  s <- melt(se[, c(id, dims), with = FALSE], id.vars = id, variable.name = "dim", value.name = "se")
  x <- n[r, on = c(id, "dim")][s, on = c(id, "dim")]
  p <- tstrsplit(as.character(x$dim), "_"); x[, contrast_id := paste(arm, p[[1]], p[[2]], p[[3]], sep = "|")][, dim := NULL]; x
}
GRsp <- rbind(long("entrez_gene", gE, rd("01_nodes_EE_raw_logFC.csv"), rd("01_nodes_EE_se.csv"), "EE", dims16),
              long("entrez_gene", gR, rd("01_nodes_RE_raw_logFC.csv"), rd("01_nodes_RE_se.csv"), "RE", dims16))
MRsp <- rbind(long("metabolite", mE, rd("01b_metab_nodes_EE_raw_logFC.csv"), rd("01b_metab_nodes_EE_se.csv"), "EE", dims9),
              long("metabolite", mR, rd("01b_metab_nodes_RE_raw_logFC.csv"), rd("01b_metab_nodes_RE_se.csv"), "RE", dims9))
stopifnot(nrow(GRsp) == 471 * 16 * 2, nrow(MRsp) == 450 * 9 * 2, !anyNA(GRsp$logFC), !anyNA(MRsp$logFC), all(c(GRsp$contrast_id, MRsp$contrast_id) %in% K$contrast_id))
wr(GRsp[, .(entrez_gene, contrast_id, logFC, normalised, se)], "rel_gene_response.csv")
wr(MRsp[, .(metabolite, contrast_id, logFC, normalised, se)], "rel_metabolite_response.csv")
message("Neo4j import files in ", NEO)
