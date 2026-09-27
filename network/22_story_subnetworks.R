#!/usr/bin/env Rscript
# =====================================================================================================
# 22_story_subnetworks.R — STEP 22: ONE PHYSICAL SUBNETWORK PER T2D / AGEING STORY (FOR SLIDES)
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   For each of the three top stories, which physical partners, metabolites and PTMs surround the proteins the
#   study found altered? One subnetwork per story, drawn the same way for all three, for slides and Cytoscape:
#     S1 muscle_t2d   — Öhman 2021 T2D muscle proteome (Amar et al. 2024 set), altered = p < 0.05 (as step 19)
#     S2 blood_t2d    — Gadd 2024 UK Biobank incident T2D plasma proteins, altered = p < 3.1e-6 (as step 20)
#     S3 blood_ageing — Sun 2023 UK Biobank plasma age proteins, altered = p < 1.7e-5 (as step 20)
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Seeds: the altered study proteins (gene symbols, one row per gene with the smallest p, as steps 19 / 20)
#      mapped to mnet protein nodes: exact symbol; else case-insensitive; else each part of a multi-gene Olink
#      assay ("A_B"); else a STRING v12 symbol alias (HGNC current / previous / alias symbol, UniProt gene name)
#      -> ENSP -> node, kept only if it resolves to one node. Direction (up / down in disease) and z are kept.
#   2. Physical layer (mnet): STRING ppi with score >= PPI_MIN, Rhea catalysis / transport; currency metabolites
#      and lipid-class edges are left out. Each ppi edge is labelled ppi_physical when the pair is in STRING's
#      physical-links file with the same score (within 0.05; ENSP -> node with interaction_db's own mapping,
#      string_network.load.load_mapping_with_fallback), else ppi_scaled (mnet's 0.9 x combined score).
#   3. Partners: non-seed nodes with >= MIN_SEED_LINKS seed neighbours. Kept edges: seed-seed and seed-partner.
#   4. PTMs per protein from mnet (phospho / glyco counts, upstream kinases, MoTrPAC sites). S1 only: the
#      MoTrPAC muscle phosphosites of the story's proteins as site nodes (protein -> site), known-kinase flag.
#   5. Tables and GraphML hold every altered protein. The PNG seeds are the altered proteins with a physical ppi
#      (>= 700) to another altered protein: up to 25 of them in our exercise network (MEASURED: the 353 nodes of
#      step 14 for now; STORY_MEASURED=genes471 switches to step 1's 471 genes), by |z|, then the top-|z| others
#      (15; S1: 5). Drawn: ppi_physical and Rhea edges, partners with >= 2 of those seeds (S1: >= 1 of them in our
#      network; raised to 3, then 4, while the graph has more than 80 nodes), metabolites only if measured or linked to >= 2 seeds, plus our measured
#      metabolites with a Rhea link to a seed (labelled), the largest connected component, and Louvain clusters laid
#      out one per box (no overlapping hulls), those with >= 5 nodes as hulls named by their 3 most connected nodes
#      (seeds first). Marks beside proteins: G = glycosylated (mnet); P (S1) = a MoTrPAC muscle phosphosite
#      responding to exercise (adj. p < 0.05, step 17s) when that file exists, else the interim "phosphosites
#      measured in MoTrPAC muscle"; S1 marks only labelled proteins. Two PNGs per story: a clean slide version
#      (network, hull names, labels, legend) and an annotated one (+ title, rules, full-network counts, caption);
#      the texts also go to <story>_slide_text.md / 22_slide_text.md, with a title read from step 20's tests.
#
# HOW TO RUN
#   After step 20 (and 14):   Rscript network/22_story_subnetworks.R   (about 1 minute; python3 must import
#   interaction_db's string_network package, i.e. pandas + pyyaml)
#
# DATA AND PROVENANCE
#   Study analytes: $HACK_OUT/20_disease_scores.csv, 20_disease_sets.csv (step 20). Physical layer: the team's
#   mnet resource at $MNET_DIR (STRING v12 ppi >= 500, Rhea, UniProt + OmniPath PTMs) and its raw STRING files at
#   $STRING_RAW (physical links, aliases, idmapping). Exercise layer (optional): step 14 nodes / edges.
#
# TECH STACK:  R 4.4; data.table, igraph, ggplot2, ggrepel, scales, patchwork; python3 (pandas, pyyaml) for the ENSP map.
#
# INPUTS:  $HACK_OUT/20_disease_scores.csv, 20_disease_sets.csv, [14_joint_nodes.csv, 14_joint_edges.csv,
#          01_nodes_EE.csv, 17_phospho_site_stats.csv, 20_tests.csv]; $MNET_DIR/edges.csv, nodes.csv, glycosites.csv (count check), proteins_ptm.csv,
#          motrpac_feature_site_map.csv; $STRING_RAW/9606.protein.physical.links.v12.0.txt.gz,
#          9606.protein.aliases.v12.0.txt.gz (+ the files interaction_db's mapping reads); interaction_db/config.yaml
# OUTPUTS: $HACK_OUT/22_story_networks/<story>_nodes.csv, <story>_edges.csv, <story>.graphml,
#          <story>.png (clean slide figure), <story>_annotated.png (with title, rules, counts), <story>_slide_text.md
#          $HACK_OUT/22_story_networks/22_summary.csv, 22_slide_text.md, 22_mapping_check.csv, 22_seed_mapping.csv, 22_ensp_to_node.csv
#
# KNOWN LIMITS
#   The subnetworks are descriptive (no test). Partners are chosen by seed links only, so large seed sets pull in
#   hubs. MEASURED = step 14's 353 network nodes marks "in our network", not "measured" in the wider sense, until
#   step 1's 471 genes are used. Blood stories: no plasma phospho data, so PTMs are node attributes only.
# =====================================================================================================

# Load packages quietly; C collation so pair keys and orderings do not depend on the machine's locale.
suppressMessages({ library(data.table); library(igraph); library(ggplot2); library(ggrepel); library(patchwork) })
invisible(Sys.setlocale("LC_COLLATE", "C"))
.args <- commandArgs(trailingOnly = FALSE); .here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", .args, value = TRUE))))

# Folders and inputs (override with environment variables).
OUT    <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
MNET   <- Sys.getenv("MNET_DIR", unset = path.expand("~/Desktop/string_database/data/output"))
RAW    <- Sys.getenv("STRING_RAW", unset = file.path(dirname(MNET), "raw"))
PY     <- Sys.getenv("PYTHON", unset = "python3")
IDB    <- file.path(.here, "interaction_db")          # read / imported only, never modified
MEASURED_SOURCE <- Sys.getenv("STORY_MEASURED", unset = "joint14")   # joint14 | genes471
DEST <- file.path(OUT, "22_story_networks"); dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

# Parameters.
PPI_MIN        <- 700      # STRING score for a ppi edge (mnet supplies >= 500)
MIN_SEED_LINKS <- 2        # a partner needs at least this many seed neighbours
PNG_MAX_NODES  <- 150      # PNG: raise the partner rule (to 3, then 4) while the drawn graph is larger
PNG_MAX_SEEDS  <- 40       # PNG: label all seeds up to this many (else the top 25 by |z|); seed quotas are per story
PNG_TARGET_CAPPED <- 80    # PNG: aim for <= 80 nodes (partner rule raised to 3, then 4)
PNG_MIN_CLUSTER <- 5       # PNG: Louvain clusters with at least this many nodes get a hull and a name
PNG_MAX_LABELS <- 35
PHOS_ALPHA     <- 0.05     # a MoTrPAC muscle phosphosite "responds" at adj. p < 0.05 after EE or RE (as step 19)
PHYS_TOL       <- 0.05     # |mnet score - physical-links score| for ppi_physical
SEED           <- 20260927
SITE_STORIES   <- "muscle_t2d"

# The three stories, in slide order: study set in step 20 and its "altered" threshold.
STORIES <- data.table(story = c("muscle_t2d", "blood_t2d", "blood_ageing"),
                      set   = c("ohman_2021", "gadd_2024_ukb_incident_T2D", "sun_2023_ukb_age"),
                      title = c("S1 · Muscle · T2D (Öhman 2021)", "S2 · Blood · future T2D (UK Biobank, Gadd 2024)",
                                "S3 · Blood · ageing (UK Biobank, Sun 2023)"),
                      alpha = c(0.05, 3.1e-6, 1.7e-5),
                      dir   = c("in T2D", "with future T2D", "with age"),   # wording of "up / down" on the slide
                      name  = c("Muscle · T2D", "Blood · future T2D", "Blood · ageing"),   # slide title stem
                      what  = c("T2D", "future-T2D", "ageing"),             # "... reverses the <what> proteins"
                      # PNG seeds: up to seeds_in in our exercise network + up to seeds_out outside it; partners need
                      # >= min_in_links links to in-network seeds (S1: the complex I subunits outside our network
                      # otherwise take over the figure)
                      seeds_in = c(25L, 25L, 25L), seeds_out = c(5L, 15L, 15L), min_in_links = c(1L, 0L, 0L))

# One palette for all three stories.
COL_ORIGIN <- c(study_up = "#D7301F", study_down = "#2B6CB0", physical_partner = "grey72", metabolite = "#2CA25F", ptm_site = "#E69F00")
COL_EDGE   <- c(ppi_physical = "grey55", rhea = "#2CA25F")

# ---- 1. inputs: print the columns of each ---------------------------------------------------------------------
show <- function(name, d) message(sprintf("%-32s %7d rows | %s", name, nrow(d), paste(names(d), collapse = ", ")))
need <- c("20_disease_scores.csv", "20_disease_sets.csv")
miss <- need[!file.exists(file.path(OUT, need))]
if (length(miss)) stop("missing in HACK_OUT: ", paste(miss, collapse = ", "),
                       " (S1 needs the Öhman table from step 20, which is only on Vidal's machine)")
DS   <- fread(file.path(OUT, "20_disease_scores.csv")); show("20_disease_scores.csv", DS)
SETS <- fread(file.path(OUT, "20_disease_sets.csv"));   show("20_disease_sets.csv", SETS)
EXE  <- all(file.exists(file.path(OUT, c("14_joint_nodes.csv", "14_joint_edges.csv"))))
if (EXE) { JN <- fread(file.path(OUT, "14_joint_nodes.csv")); show("14_joint_nodes.csv", JN)
           JE <- fread(file.path(OUT, "14_joint_edges.csv")); show("14_joint_edges.csv", JE) } else message("exercise layer (step 14): not found, skipped")
mn <- c("edges", "nodes", "glycosites", "proteins_ptm", "motrpac_feature_site_map")
MN <- setNames(lapply(mn, function(f) fread(file.path(MNET, paste0(f, ".csv")))), mn)
for (f in mn) show(paste0("mnet/", f, ".csv"), MN[[f]])
ND <- MN$nodes

# The thresholds must agree with step 20's set table.
chk <- merge(STORIES, SETS[, .(set, alpha_20 = alpha)], by = "set", sort = FALSE)
stopifnot(nrow(chk) == 3, isTRUE(all.equal(chk$alpha, chk$alpha_20)))

# MEASURED: the protein symbols the PNG may use as seeds.
MEASURED <- if (MEASURED_SOURCE == "genes471") {
  f <- file.path(OUT, "01_nodes_EE.csv"); if (!file.exists(f)) stop("STORY_MEASURED=genes471 needs ", f); fread(f)$gene_symbol
} else if (EXE) JN[node_type == "protein", node] else stop("STORY_MEASURED=joint14 needs 14_joint_nodes.csv")
message(sprintf("MEASURED (%s): %d protein symbols", MEASURED_SOURCE, length(MEASURED)))

# Exercise-responsive MoTrPAC muscle phosphosites (step 17s), for the S1 "P" marks; optional.
PHOS_FILE <- file.path(OUT, "17_phospho_site_stats.csv"); HAVE_PHOS <- file.exists(PHOS_FILE)
RESP_PHOS <- if (HAVE_PHOS) fread(PHOS_FILE)[tissue == "muscle" & arm %in% c("EE", "RE") & adj_p < PHOS_ALPHA, unique(protein)] else character(0)
message(if (HAVE_PHOS) sprintf("17_phospho_site_stats.csv: %d proteins with a responding muscle site", length(RESP_PHOS)) else "17_phospho_site_stats.csv: not found, no P marks")

# ---- 2. ENSP -> mnet node, with interaction_db's own mapping (imported, not re-implemented) ---------------------
ENSP_MAP <- file.path(DEST, "22_ensp_to_node.csv")
py <- sprintf(paste(
  "import sys, yaml",
  "sys.path.insert(0, %s)",
  "from string_network.load import load_mapping_with_fallback",
  "cfg = yaml.safe_load(open(%s)); cfg['paths']['raw'] = %s",
  "m = load_mapping_with_fallback(cfg, cfg.get('mapping_source', 'uniprot_idmapping'), int(cfg.get('taxon', 9606)),",
  "    alias_sources=cfg.get('mapping_alias_sources', ['UniProt_AC', 'Ensembl_UniProt']), return_source=True)",
  "m.rename(columns={'uniprot': 'node_id'}).to_csv(%s, index=False)", sep = "\n"),
  shQuote(file.path(IDB, "src"), "sh"), shQuote(file.path(IDB, "config.yaml"), "sh"), shQuote(RAW, "sh"), shQuote(ENSP_MAP, "sh"))
if (system2(PY, c("-c", shQuote(py)), stdout = "", stderr = "") != 0) stop("ENSP mapping via interaction_db failed (", PY, ")")
E2N <- fread(ENSP_MAP)
message(sprintf("ENSP -> node (interaction_db mapping): %d ENSPs; by source: %s", nrow(E2N),
                paste(E2N[, .N, by = source][, sprintf("%s %d", source, N)], collapse = ", ")))
e2n <- function(x) E2N$node_id[match(sub("^9606\\.", "", x), E2N$ensp)]

# ---- 3. seeds: map each story's analytes to mnet protein node_ids ------------------------------------------------
P <- ND[node_type == "protein" & !is.na(gene_symbol) & gene_symbol != "", .(node_id, gene_symbol)]
map_one <- function(g) {
  hit <- P[gene_symbol == g]; if (nrow(hit)) return(hit[, .(node_id, how = "exact")])
  hit <- P[toupper(gene_symbol) == toupper(g)]; if (nrow(hit)) return(hit[, .(node_id, how = "case-insensitive")])
  parts <- unlist(strsplit(g, "[;_,]")); if (length(parts) > 1) { hit <- P[gene_symbol %in% parts]
    if (nrow(hit)) return(hit[, .(node_id, how = "multi-gene part")]) }
  data.table(node_id = NA_character_, how = "unmapped")
}
MAP <- rbindlist(lapply(seq_len(nrow(STORIES)), function(i) {
  s <- STORIES[i]; d <- DS[set == s$set & !is.na(p)][order(p)][!duplicated(gene)]   # one row per gene, smallest p
  m <- rbindlist(lapply(d$gene, function(g) cbind(gene = g, map_one(g))))
  merge(d[, .(gene, logFC, p, z, altered = p < s$alpha)], m, by = "gene")[, story := s$story]
}))
# STRING symbol aliases for what is still unmapped: symbol -> ENSP -> node, only if exactly one mnet node.
ALIAS_SRC <- c("Ensembl_HGNC_symbol", "Ensembl_HGNC_prev_symbol", "Ensembl_HGNC_alias_symbol", "Ensembl_external_synonym_HGNC",
               "BioMart_HUGO", "Ensembl_HGNC", "UniProt_GN_Name", "UniProt_GN_Synonyms")
todo <- unique(MAP[how == "unmapped", gene])
AL <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(RAW, "9606.protein.aliases.v12.0.txt.gz"))), sep = "\t", header = TRUE,
            col.names = c("ensp", "alias", "source"))[alias %in% todo]
AL <- AL[vapply(strsplit(source, " "), function(x) any(x %in% ALIAS_SRC), logical(1))]
AL[, node_id := e2n(ensp)]
AL <- unique(AL[!is.na(node_id) & node_id %in% ND[node_type == "protein", node_id], .(gene = alias, node_id)])
AL <- AL[, if (.N == 1) .SD, by = gene]
MAP[how == "unmapped" & gene %in% AL$gene, `:=`(node_id = AL$node_id[match(gene, AL$gene)], how = "string alias")]
fwrite(MAP[altered == TRUE][order(match(story, STORIES$story), p)], file.path(DEST, "22_seed_mapping.csv"))

SUM_MAP <- MAP[altered == TRUE, .(altered = uniqueN(gene), altered_up = uniqueN(gene[z > 0]), altered_down = uniqueN(gene[z < 0]),
                                  mapped = uniqueN(gene[how != "unmapped"]), by_alias = uniqueN(gene[how == "string alias"]),
                                  unmapped = uniqueN(gene[how == "unmapped"]), seed_nodes = uniqueN(na.omit(node_id))), by = story]
SUM_MAP[, pct_mapped := round(100 * mapped / altered, 1)]
message("\nmapping check (altered study proteins -> mnet):"); print(SUM_MAP)
for (s in STORIES$story) { u <- MAP[story == s & altered & how == "unmapped", gene]; r <- MAP[story == s & altered & how == "string alias", gene]
  message(sprintf("  %s: recovered by alias %d (%s)\n  %s: still unmapped %d (%s)", s, length(r), paste(r, collapse = ", "),
                  s, length(u), paste(u, collapse = ", "))) }
fwrite(SUM_MAP, file.path(DEST, "22_mapping_check.csv"))

# ---- 4. physical layer: ppi (physical / scaled) and Rhea --------------------------------------------------------
cur <- ND[is_currency == TRUE, node_id]
E <- MN$edges[edge_type %in% c("ppi", "catalysis", "transport") & !(node1 %in% cur) & !(node2 %in% cur)]
PPI <- E[edge_type == "ppi" & score >= PPI_MIN, .(node1, node2, score)][, id := .I]
PH <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(RAW, "9606.protein.physical.links.v12.0.txt.gz"))), sep = " ")
PH <- PH[, .(a = e2n(protein1), b = e2n(protein2), phys = combined_score)][!is.na(a) & !is.na(b) & a != b]
PH <- unique(rbind(PH, PH[, .(a = b, b = a, phys)]))
hit <- merge(PPI, PH, by.x = c("node1", "node2"), by.y = c("a", "b"), allow.cartesian = TRUE)
phys_id <- hit[abs(score - phys) <= PHYS_TOL, unique(id)]
PPI[, edge_type := fifelse(id %in% phys_id, "ppi_physical", "ppi_scaled")]
# Consistency with the score value: a physical score is an integer; 0.9 x an integer is a non-integer or a multiple of 9.
PPI[, value_says := fifelse(abs(score - round(score)) > 1e-6, "non-integer", fifelse(round(score) %% 9 == 0, "integer, x9", "integer, not x9"))]
message(sprintf("\nppi >= %d: %d edges; ppi_physical %d, ppi_scaled %d (pair in physical file but score differs: %d)", PPI_MIN, nrow(PPI),
                PPI[edge_type == "ppi_physical", .N], PPI[edge_type == "ppi_scaled", .N], hit[!id %in% phys_id, uniqueN(id)]))
print(dcast(PPI[, .N, by = .(edge_type, value_says)], edge_type ~ value_says, value.var = "N", fill = 0))
RH <- E[edge_type %in% c("catalysis", "transport"), .(score = max(score), rhea_kind = paste(sort(unique(edge_type)), collapse = "+")), by = .(node1, node2)]
PE <- rbind(PPI[, .(node1, node2, edge_type, score, rhea_kind = NA_character_)], RH[, .(node1, node2, edge_type = "rhea", score, rhea_kind)])
PE[, key := fifelse(node1 < node2, paste(node1, node2), paste(node2, node1))]
PE <- PE[order(key, edge_type)][!duplicated(key)][, key := NULL]

# ---- 5. node attributes shared by all stories --------------------------------------------------------------------
PTM <- MN$proteins_ptm[, .(node_id, n_phospho = n_phosphosites, n_phospho_with_kinase = n_phosphosites_with_kinase, n_glyco = n_glycosites,
                           n_N_linked, n_O_linked, n_O_GlcNAc, upstream_kinases, is_kinase, n_motrpac_sites)]
stopifnot(!anyDuplicated(PTM$node_id))
SITES <- MN$motrpac_feature_site_map[tissue %in% c("muscle", "both") & mapping_status %in% c("canonical_ok", "isoform_mapped") &
                                     !is.na(node_id) & node_id != "", .(known_kinase = any(has_known_kinase)), by = .(site_id, node_id)]
NA_ <- ND[, .(node_id, node_type, label = fifelse(node_type == "protein" & !is.na(gene_symbol) & gene_symbol != "", gene_symbol, label),
              is_measured)]
if (EXE) { RESP <- JN[, .(key = tolower(node), resp_EE, resp_RE)]
           JE[, key := fifelse(tolower(node_a) < tolower(node_b), paste(tolower(node_a), tolower(node_b)), paste(tolower(node_b), tolower(node_a)))] }

# ---- 6. build one story --------------------------------------------------------------------------------------------
build <- function(st) {
  sd <- MAP[story == st & altered & how != "unmapped"][order(-abs(z))][!duplicated(node_id)]   # one row per node (largest |z|)
  S <- sd$node_id
  inc <- PE[node1 %in% S | node2 %in% S][, `:=`(s1 = node1 %in% S, s2 = node2 %in% S)]
  pl <- rbind(inc[s1 & !s2, .(partner = node2, seed = node1)], inc[s2 & !s1, .(partner = node1, seed = node2)])
  pc <- pl[, .(n_seed_links = uniqueN(seed)), by = partner][n_seed_links >= MIN_SEED_LINKS]
  ed <- inc[(s1 & s2) | (s1 & node2 %in% pc$partner) | (s2 & node1 %in% pc$partner)]
  ed[, edge_origin := fifelse(s1 & s2, "study_study", "study_partner")][, c("s1", "s2") := NULL]
  ss <- rbind(ed[edge_origin == "study_study", .(n = node1, m = node2)], ed[edge_origin == "study_study", .(n = node2, m = node1)])[, .(n_seed_links = uniqueN(m)), by = n]
  nodes <- rbind(sd[, .(node_id, study_gene = gene, z, logFC, p, map_how = how, node_origin = fifelse(z > 0, "study_up", "study_down"))],
                 data.table(node_id = pc$partner), fill = TRUE)
  nodes <- merge(nodes, NA_, by = "node_id", all.x = TRUE, sort = FALSE)
  nodes[is.na(node_origin), node_origin := fifelse(node_type == "metabolite", "metabolite", "physical_partner")]
  nodes[, n_seed_links := c(ss$n_seed_links, pc$n_seed_links)[match(node_id, c(ss$n, pc$partner))]][is.na(n_seed_links), n_seed_links := 0L]
  nodes[, in_motrpac := fifelse(node_type == "protein", label %in% MEASURED | study_gene %in% MEASURED, is_measured %in% TRUE)]
  nodes <- merge(nodes, PTM, by = "node_id", all.x = TRUE, sort = FALSE)
  if (EXE) { nodes[, `:=`(resp_EE = RESP$resp_EE[match(tolower(label), RESP$key)], resp_RE = RESP$resp_RE[match(tolower(label), RESP$key)])] }
  # S1: MoTrPAC muscle phosphosites of the story's proteins as site nodes
  sites <- if (st %in% SITE_STORIES) SITES[node_id %in% nodes[node_type == "protein", node_id]] else SITES[0]
  if (nrow(sites)) {
    nodes <- rbind(nodes, sites[, .(node_id = site_id, label = site_id, node_type = "ptm_site", node_origin = "ptm_site", known_kinase,
                                    in_motrpac = TRUE, n_seed_links = 0L)], fill = TRUE)
    ed <- rbind(ed, sites[, .(node1 = node_id, node2 = site_id, edge_type = "protein_site", score = NA_real_, rhea_kind = NA_character_,
                              edge_origin = "protein_site")], fill = TRUE)
  }
  lab <- setNames(nodes$label, nodes$node_id)
  if (EXE) { k <- fifelse(tolower(lab[ed$node1]) < tolower(lab[ed$node2]), paste(tolower(lab[ed$node1]), tolower(lab[ed$node2])),
                          paste(tolower(lab[ed$node2]), tolower(lab[ed$node1])))
             ed[, `:=`(w_EE = JE$w_EE[match(k, JE$key)], w_RE = JE$w_RE[match(k, JE$key)])] }
  ed[, `:=`(label1 = lab[node1], label2 = lab[node2])]
  g <- graph_from_data_frame(ed[, .(node1, node2)], directed = FALSE, vertices = nodes[, .(node_id)])
  cp <- components(g); nodes[, in_lcc := cp$membership[node_id] == which.max(cp$csize)]
  ed[, in_lcc := nodes$in_lcc[match(node1, nodes$node_id)]]
  list(nodes = nodes, edges = ed, g = g, lcc = max(cp$csize))
}

# Slide title from step 20's tests (never hardcoded): which arm reverses the story's proteins, with the EE - RE p.
TESTS20 <- if (file.exists(file.path(OUT, "20_tests.csv"))) fread(file.path(OUT, "20_tests.csv")) else NULL
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p %s", formatC(signif(p, 2), format = "fg"))
slide_title <- function(s) {
  t <- if (!is.null(TESTS20)) TESTS20[set == s$set & test == "A reversal"] else NULL
  if (is.null(t) || nrow(t) != 1 || anyNA(t[, .(diff, p_diff)])) return(s$name)
  if (t$p_diff >= 0.05) return(sprintf("%s: no clear difference between endurance and resistance (%s)", s$name, fmt_p(t$p_diff)))
  sprintf("%s: %s reverses the %s proteins (difference %s)", s$name, if (t$diff > 0) "endurance" else "resistance", s$what, fmt_p(t$p_diff))
}
message(if (is.null(TESTS20)) "20_tests.csv: not found, slide titles fall back to the story names" else "20_tests.csv: slide titles from the reversal tests")

# ---- 7. the slide figure --------------------------------------------------------------------------------------------
# Hull around a cluster (padded so small clusters still show a shape).
hull <- function(d, pad = 0.35) {
  a <- seq(0, 2 * pi, length.out = 13)[-13]
  pts <- d[, .(x = rep(x, each = 12) + pad * cos(a), y = rep(y, each = 12) + pad * sin(a))]
  pts[chull(pts$x, pts$y)]
}
# Cluster layout: Louvain first, then each cluster laid out on its own and the clusters packed in rows (largest first)
# into a 16:9 frame, so hulls never overlap. Box side ~ sqrt(cluster size); coordinates in box units.
pack_layout <- function(g, cl) {
  sz <- sort(table(cl), decreasing = TRUE); side <- 1.6 * sqrt(as.numeric(sz)) + 1.2
  W <- sqrt(sum(side^2) * 16 / 9) * 1.05; x0 <- 0; y0 <- 0; row_h <- 0; pos <- list()
  for (i in seq_along(sz)) {
    if (x0 > 0 && x0 + side[i] > W) { x0 <- 0; y0 <- y0 - row_h; row_h <- 0 }
    v <- names(cl)[cl == names(sz)[i]]; sg <- induced_subgraph(g, v)
    set.seed(SEED); xy <- if (length(v) > 1) layout_with_fr(sg, niter = 1500) else matrix(0, 1, 2)
    xy <- apply(xy, 2, function(u) if (diff(range(u)) > 0) (u - min(u)) / diff(range(u)) - 0.5 else u * 0)
    inner <- side[i] - 1.2   # 0.6 margin on each side for the hull and its label
    pos[[i]] <- data.table(node_id = V(sg)$name, x = x0 + side[i] / 2 + xy[, 1] * inner, y = y0 - side[i] / 2 + xy[, 2] * inner * 0.85)
    x0 <- x0 + side[i]; row_h <- max(row_h, side[i])
  }
  rbindlist(pos)
}
draw <- function(s, B) {
  N <- B$nodes; st_nodes <- N[node_origin %in% c("study_up", "study_down")]
  # seeds: altered proteins with >= 1 physical ppi (>= PPI_MIN) to another altered protein; up to s$seeds_in in our
  # exercise network (by |z|), then the top-|z| others up to s$seeds_in + s$seeds_out in all
  cand_ids <- B$edges[edge_origin == "study_study" & edge_type == "ppi_physical", unique(c(node1, node2))]
  cand <- st_nodes[node_id %in% cand_ids][order(-abs(z), node_id)]
  s_in <- head(cand[in_motrpac == TRUE, node_id], s$seeds_in)
  s_out <- head(cand[in_motrpac == FALSE, node_id], s$seeds_in + s$seeds_out - length(s_in))
  sd <- c(s_in, s_out); target <- PNG_TARGET_CAPPED
  pe <- B$edges[edge_type %in% c("ppi_physical", "rhea")]
  inc <- pe[node1 %in% sd | node2 %in% sd][, `:=`(s1 = node1 %in% sd, s2 = node2 %in% sd)]
  pl <- rbind(inc[s1 & !s2, .(partner = node2, seed = node1)], inc[s2 & !s1, .(partner = node1, seed = node2)])
  pl <- pl[, .(k = uniqueN(seed), k_in = uniqueN(seed[seed %in% s_in])), by = partner]
  pl <- merge(pl, N[, .(partner = node_id, node_type, is_measured)], by = "partner", sort = FALSE)
  for (kmin in sort(unique(c(MIN_SEED_LINKS, 3, 4)))) {
    # partners: >= kmin seeds, of which >= s$min_in_links in our exercise network
    cands <- pl[k >= kmin & k_in >= s$min_in_links]
    # metabolites only if measured in MoTrPAC or linked to >= 2 seeds
    met_drop <- cands[node_type == "metabolite" & !(is_measured %in% TRUE | k >= 2), .N]
    pt <- cands[node_type != "metabolite" | is_measured %in% TRUE | k >= 2, partner]
    e <- inc[(s1 & s2) | (s1 & node2 %in% pt) | (s2 & node1 %in% pt)]
    g <- graph_from_data_frame(e[, .(node1, node2)], directed = FALSE, vertices = data.table(node_id = unique(c(sd, pt))))
    cp <- components(g); keep <- names(cp$membership)[cp$membership == which.max(cp$csize)]
    if (length(keep) <= target) break
  }
  # our measured metabolites with a Rhea link to >= 1 seed (currency metabolites are already out of the edge list)
  mm <- setdiff(pl[node_type == "metabolite" & is_measured %in% TRUE, partner], pt)
  mm <- mm[mm %in% inc[edge_type == "rhea", c(node1, node2)]]
  if (length(mm)) {
    e <- inc[(s1 & s2) | (s1 & node2 %in% pt) | (s2 & node1 %in% pt) | (edge_type == "rhea" & (node1 %in% mm | node2 %in% mm) & (s1 | s2))]
    g <- graph_from_data_frame(e[, .(node1, node2)], directed = FALSE, vertices = data.table(node_id = unique(c(sd, pt, mm))))
    cp <- components(g); keep <- names(cp$membership)[cp$membership == which.max(cp$csize)]
  }
  g <- induced_subgraph(g, keep); e <- e[node1 %in% keep & node2 %in% keep]
  set.seed(SEED); cl <- setNames(as.integer(membership(cluster_louvain(g))), V(g)$name)
  V <- pack_layout(g, cl)[, `:=`(cluster = cl[node_id], deg = degree(g)[node_id])]
  V[, deg_in := vapply(node_id, function(v) sum(cl[names(neighbors(g, v))] == cl[[v]]), integer(1))]
  V <- merge(V, N[, .(node_id, label, node_type, node_origin, z, in_motrpac, is_measured, n_glyco)], by = "node_id", sort = FALSE)
  V[, `:=`(png_seed = node_id %in% sd, k = pl$k[match(node_id, pl$partner)], measured_met = node_type == "metabolite" & is_measured %in% TRUE)]
  V[, size := fifelse(png_seed, scales::rescale(abs(z), to = c(3.5, 10), from = range(abs(z[png_seed]))), fifelse(node_type == "metabolite", 3.2, 2.8))]
  V[, shape := fifelse(node_type == "metabolite", 24L, 21L)]
  es <- e[, .(edge_type, x = V$x[match(node1, V$node_id)], y = V$y[match(node1, V$node_id)],
              xend = V$x[match(node2, V$node_id)], yend = V$y[match(node2, V$node_id)])]
  # clusters with >= 5 nodes: light hull, named by their 3 most connected nodes (seeds first, proteins before metabolites)
  big <- V[, .N, by = cluster][N >= PNG_MIN_CLUSTER][order(-N, cluster)]
  H <- rbindlist(lapply(big$cluster, function(c) hull(V[cluster == c])[, cluster := c]))
  CN <- V[cluster %in% big$cluster][order(cluster, -png_seed, node_type != "protein", -deg_in, label)][, .(name = paste(head(label, 3), collapse = " · ")), by = cluster]
  CN <- merge(merge(big, CN, by = "cluster"), H[, .(x = mean(range(x)), y = max(y)), by = cluster], by = "cluster")[order(-N, cluster)]
  # labels: all seeds if <= 40 drawn (else top 25 by |z|), then partners linked to >= 3 seeds (at most PNG_MAX_LABELS),
  # plus every drawn measured metabolite
  sl <- V[png_seed == TRUE][order(-abs(z))]; if (nrow(sl) > PNG_MAX_SEEDS) sl <- head(sl, 25)
  L <- rbind(sl, V[png_seed == FALSE & k >= 3 & !measured_met][order(-k, label)])[seq_len(min(.N, PNG_MAX_LABELS))]
  L <- rbind(L, V[measured_met == TRUE])
  # PTM marks: G = glycosylated in mnet; P (S1, labelled proteins only) = exercise-responsive muscle phosphosite (step 17s,
  # as step 19) when that file exists, else (interim) = phosphosites measured in MoTrPAC muscle (the story's site nodes)
  np <- B$edges[edge_type == "protein_site", .N, by = node1]
  V[, mark_G := node_type == "protein" & !is.na(n_glyco) & n_glyco > 0]
  V[, has_P := node_type == "protein" & s$story %in% SITE_STORIES & (if (HAVE_PHOS) label %in% RESP_PHOS else node_id %in% np$node1)]
  V[, mark_P := has_P & node_id %in% L$node_id]

  # ---- the texts (annotated PNG and the slide-text sidecar) ----
  n_up <- V[png_seed & node_origin == "study_up", .N]; n_dn <- V[png_seed & node_origin == "study_down", .N]; ns <- n_up + n_dn
  partner_rule <- sprintf(">= %d seeds%s", kmin, if (s$min_in_links > 0) sprintf(" (>= %d in our exercise network)", s$min_in_links) else "")
  sub <- sprintf("Seeds: altered proteins with a physical link to another altered protein (%d); the top %d in our exercise network + the top %d outside it, by |z|; %d drawn, %d of them in our exercise network (%d up, %d down %s). Partners linked to %s: %d (%d unaltered proteins, %d other altered proteins, %d metabolites). %d physical PPI + %d Rhea edges; largest component, %d nodes%s; %d Louvain clusters >= %d nodes outlined. %d measured metabolites (Rhea). G: %d of %d drawn seeds glycosylated.%s",
                 nrow(cand), length(s_in), length(s_out), ns, V[png_seed & in_motrpac, .N], n_up, n_dn, s$dir, partner_rule, V[png_seed == FALSE, .N], V[node_origin == "physical_partner", .N],
                 V[png_seed == FALSE & node_origin %in% c("study_up", "study_down"), .N], V[node_type == "metabolite", .N],
                 e[edge_type == "ppi_physical", .N], e[edge_type == "rhea", .N], nrow(V),
                 if (nrow(V) > target) sprintf(" (above the %d target even at >= %d)", target, kmin) else "", nrow(big), PNG_MIN_CLUSTER,
                 V[measured_met == TRUE, .N], V[png_seed & mark_G, .N], ns,
                 if (s$story %in% SITE_STORIES) sprintf(" P: %d of %d drawn seeds with %s (marked on labelled proteins).", V[png_seed & has_P, .N], ns,
                                                        if (HAVE_PHOS) "exercise-responsive muscle phosphosites" else "phosphosites measured in MoTrPAC muscle") else "")
  wrap <- function(x, w) paste(strwrap(x, w), collapse = "\n")
  prot <- N[node_type == "protein"]
  n_sites <- if (s$story %in% SITE_STORIES) sprintf("%d MoTrPAC muscle phosphosites", N[node_origin == "ptm_site", .N]) else sprintf("%d known phosphosites (mnet)", sum(prot$n_phospho, na.rm = TRUE))
  strip <- sprintf("Full story network: %d study proteins (%d in our exercise network) · %d physical partners · %d metabolites (%d measured by us) · %s · %d glycosylated proteins — explore all of it in the interactive network.",
                   nrow(st_nodes), st_nodes[in_motrpac == TRUE, .N], N[node_origin == "physical_partner", .N], N[node_origin == "metabolite", .N],
                   N[node_origin == "metabolite" & is_measured %in% TRUE, .N], n_sites, prot[!is.na(n_glyco) & n_glyco > 0, .N])
  cap <- sprintf("Story statistics use only the outlined proteins; the others show the wider disease module from the study. %d seeds with no physical partner here. Edges: STRING v12 physical subnetwork >= %d, Rhea. Metabolites only if measured in MoTrPAC or linked to >= 2 seeds. Outlined = in %s. ppi_scaled edges, all altered proteins and PTM sites are in the tables / GraphML.",
                 length(sd) - sum(V$png_seed), PPI_MIN, if (MEASURED_SOURCE == "genes471") "the 471 measured genes" else "the step 14 network (metabolites: measured in MoTrPAC)")
  vst <- V[node_origin %in% c("study_up", "study_down")]
  oneline <- sprintf("%d study proteins (%d in our exercise network) · %d physical partners · %d measured metabolites · %d glycosylated",
                     nrow(vst), vst[in_motrpac == TRUE, .N], V[node_origin == "physical_partner", .N], V[measured_met == TRUE, .N], V[mark_G == TRUE, .N])
  top3 <- head(CN, 3)[, sprintf("C%d (%d nodes): %s", seq_len(.N), N, name)]

  # ---- the figure: CLEAN (network, hull names, labels, legend) or ANNOTATED (+ title, subtitle, strip, caption) ----
  render <- function(clean) {
    f <- if (clean) 1.3 else 1                      # text scale
    ux <- diff(range(V$x)) * 1.1 / 385; uy <- diff(range(V$y)) * 1.1 / (if (clean) 205 else 150)   # data units per mm of panel
    V[, `:=`(mdx = (size * f^0.5 / 2 + 1) * ux, mdy = (size * f^0.5 / 2.6 + 0.6) * uy)]
    cw <- 1.2 * f * ux; lh <- 4 * f * uy
    CN[, `:=`(x1 = x - nchar(name) * cw / 2, x2 = x + nchar(name) * cw / 2, y1 = y, y2 = y + lh)]
    ok <- logical(nrow(CN)); for (i in seq_len(nrow(CN))) ok[i] <- !any(ok[seq_len(i - 1)] & CN$x1[i] < CN$x2[seq_len(i - 1)] &
                                                                      CN$x2[i] > CN$x1[seq_len(i - 1)] & CN$y1[i] < CN$y2[seq_len(i - 1)] & CN$y2[i] > CN$y1[seq_len(i - 1)])
    HL <- CN[ok]
    lk <- V[1]; up_l <- paste("up", s$dir); dn_l <- paste("down", s$dir); out_l <- "outlined = in our exercise network"
    part_l <- paste("unaltered protein physically linked to", partner_rule)
    met_lv <- c(if (V[measured_met == TRUE, .N]) "measured by us (outlined)", if (V[node_type == "metabolite" & !measured_met, .N]) "other metabolite")
    ptm_key <- if (s$story %in% SITE_STORIES) sprintf("P = %s · G = glycosylated", if (HAVE_PHOS) "exercise-responsive phosphosites in MoTrPAC muscle" else "phosphosites measured in MoTrPAC muscle") else "G = glycosylated"
    gl <- function(order, ...) guide_legend(order = order, title.position = "left", nrow = 1, ...)
    p <- ggplot() +
      { if (nrow(H)) list(geom_polygon(data = H, aes(x, y, group = cluster), fill = "grey55", alpha = 0.07, colour = "grey70", linewidth = 0.3),
                          geom_text(data = HL, aes(x, y, label = name), vjust = -0.4, size = 3.3 * f, colour = "grey30", fontface = "italic")) } +
      geom_segment(data = es, aes(x, y, xend = xend, yend = yend, linetype = edge_type), colour = COL_EDGE[es$edge_type], linewidth = 0.35, alpha = 0.7) +
      geom_point(data = V, aes(x, y), fill = COL_ORIGIN[V$node_origin], colour = ifelse(V$in_motrpac, "black", "transparent"), shape = V$shape,
                 size = V$size * f^0.5, stroke = 0.5) +
      geom_text(data = V[mark_P == TRUE], aes(x + mdx, y + mdy, label = "P"), size = 2.9 * f, colour = "#E66100", fontface = "bold") +
      geom_text(data = V[mark_G == TRUE], aes(x + mdx, y - mdy, label = "G"), size = 2.9 * f, colour = "#7B3294", fontface = "bold") +
      geom_text_repel(data = L, aes(x, y, label = label), size = 3.6 * f, seed = SEED, max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0,
                      segment.colour = "grey40", segment.size = 0.25) +
      # legend-only layers (size 0: nothing drawn in the panel)
      geom_point(data = lk[rep(1, 3)][, g1 := factor(c(up_l, dn_l, out_l), levels = c(up_l, dn_l, out_l))], aes(x, y, fill = g1), size = 0) +
      { if (V[node_origin == "physical_partner", .N]) geom_point(data = lk, aes(x, y, colour = part_l), size = 0) } +
      { if (length(met_lv)) geom_point(data = lk[rep(1, length(met_lv))][, g3 := met_lv], aes(x, y, shape = g3), size = 0) } +
      geom_point(data = lk, aes(x, y, alpha = ptm_key), size = 0) +
      scale_fill_manual(values = setNames(c(COL_ORIGIN[["study_up"]], COL_ORIGIN[["study_down"]], "white"), c(up_l, dn_l, out_l)), name = "Study proteins",
                        guide = gl(1, override.aes = list(shape = 21, size = 4 * f, colour = c("transparent", "transparent", "black"), stroke = 0.8))) +
      scale_colour_manual(values = setNames(COL_ORIGIN[["physical_partner"]], part_l), name = "Physical partners",
                          guide = gl(2, override.aes = list(shape = 16, size = 3.5 * f))) +
      { if (length(met_lv)) scale_shape_manual(values = c("measured by us (outlined)" = 24, "other metabolite" = 24), name = "Metabolites",
                         guide = gl(3, override.aes = list(size = 3.8 * f, fill = COL_ORIGIN[["metabolite"]], colour = ifelse(met_lv == "measured by us (outlined)", "black", "transparent")))) } +
      scale_alpha_manual(values = 1, name = "PTM marks", guide = gl(4, override.aes = list(shape = NA))) +
      scale_linetype_manual(values = c(ppi_physical = "solid", rhea = "dashed"), labels = c(ppi_physical = "physical PPI (STRING)", rhea = "Rhea reaction"), name = "Edges",
                            guide = gl(5, override.aes = list(colour = unname(COL_EDGE[intersect(names(COL_EDGE), es$edge_type)])))) +
      coord_cartesian(clip = "off") + theme_void(base_size = 15) +
      theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom", legend.text = element_text(size = 11 * f),
            legend.title = element_text(size = 11 * f, face = "bold"), legend.box = "vertical", legend.box.just = "left",
            legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 0, 0))
    if (clean) return(p + theme(plot.margin = margin(25, 35, 15, 35)))
    p <- p + labs(title = s$title, subtitle = wrap(sub, 175)) +
      theme(plot.title = element_text(face = "bold", size = 22), plot.subtitle = element_text(size = 12, colour = "grey25", margin = margin(b = 10)),
            plot.margin = margin(20, 30, 4, 30), plot.title.position = "plot")
    # count strip and caption as their own rows below the plot (strip above caption), so they never overlap the legend
    txt <- function(x, size, face, col) wrap_elements(full = grid::textGrob(x, x = unit(30, "pt"), hjust = 0, gp = grid::gpar(fontsize = size, fontface = face, col = col)))
    p / txt(strip, 9.3, "bold", "grey15") / txt(wrap(cap, 235), 9.5, "plain", "grey35") + plot_layout(heights = c(1, 0.035, 0.06)) &
      theme(plot.background = element_rect(fill = "white", colour = NA))
  }
  f_clean <- file.path(DEST, paste0(s$story, ".png")); f_ann <- file.path(DEST, paste0(s$story, "_annotated.png"))
  ggsave(f_clean, render(TRUE), width = 16, height = 9, dpi = 200, bg = "white")
  ggsave(f_ann, render(FALSE), width = 16, height = 9, dpi = 200, bg = "white")

  # ---- slide text sidecar ----
  ttl <- slide_title(s)
  md <- c(sprintf("## %s", ttl), "", sprintf("**In this figure:** %s", oneline), "", sprintf("**%s**", strip), "",
          "**Top clusters**", "", paste0("- ", if (length(top3)) top3 else "none with >= 5 nodes"), "", sprintf("*%s*", cap), "",
          sprintf("<sub>Figure rules: %s</sub>", sub), "")
  writeLines(md, file.path(DEST, paste0(s$story, "_slide_text.md")))
  list(file = f_clean, md = md, png_seed_candidates = nrow(cand), png_seeds_in_network = length(s_in), png_seeds_outside = length(s_out),
       png_seeds_drawn = sum(V$png_seed), png_seeds_drawn_outlined = V[png_seed & in_motrpac, .N], png_rule = kmin, png_min_in_links = s$min_in_links,
       png_nodes = nrow(V), png_partners = V[png_seed == FALSE, .N], png_metabolites = V[node_type == "metabolite", .N], png_measured_metabolites = V[measured_met == TRUE, .N],
       png_metab_dropped = met_drop, png_edges = nrow(e), png_unconnected_seeds = length(sd) - sum(V$png_seed), png_clusters = nrow(big),
       png_labels = nrow(L), png_seeds_glyco = V[png_seed & mark_G, .N],
       png_seeds_P = if (s$story %in% SITE_STORIES) V[png_seed & has_P, .N] else NA_integer_, png_P_marks = V[mark_P == TRUE, .N],
       P_mark = if (s$story %in% SITE_STORIES) (if (HAVE_PHOS) "responsive" else "measured (interim)") else "",
       slide_title = ttl, top_clusters = paste(top3, collapse = " | "))
}

# ---- 8. run all three, write outputs and the summary ---------------------------------------------------------------
na2blank <- function(d) { d <- copy(d); for (j in names(d)) if (is.character(d[[j]]) || is.logical(d[[j]])) set(d, j = j, value = { v <- as.character(d[[j]]); v[is.na(v)] <- ""; v }); d }
SLIDE_MD <- list()
SUMMARY <- rbindlist(lapply(seq_len(nrow(STORIES)), function(i) {
  s <- STORIES[i]; B <- build(s$story)
  fwrite(B$nodes, file.path(DEST, paste0(s$story, "_nodes.csv"))); fwrite(B$edges, file.path(DEST, paste0(s$story, "_edges.csv")))
  g <- graph_from_data_frame(na2blank(B$edges), directed = FALSE, vertices = na2blank(B$nodes))
  write_graph(g, file.path(DEST, paste0(s$story, ".graphml")), format = "graphml")
  D <- draw(s, B); n <- B$nodes; e <- B$edges; SLIDE_MD[[s$story]] <<- D$md
  data.table(story = s$story, altered = SUM_MAP[story == s$story, altered], mapped = SUM_MAP[story == s$story, mapped],
             seeds = n[node_origin %in% c("study_up", "study_down"), .N], partners = n[node_origin == "physical_partner", .N],
             metabolites = n[node_origin == "metabolite", .N], ptm_sites = n[node_origin == "ptm_site", .N],
             e_ppi_physical = e[edge_type == "ppi_physical", .N], e_ppi_scaled = e[edge_type == "ppi_scaled", .N], e_rhea = e[edge_type == "rhea", .N],
             e_protein_site = e[edge_type == "protein_site", .N], lcc = B$lcc, as.data.table(D[!names(D) %in% c("file", "md")]))
}))
fwrite(SUMMARY, file.path(DEST, "22_summary.csv"))
writeLines(c("# Story subnetworks: slide text", "", unlist(SLIDE_MD)), file.path(DEST, "22_slide_text.md"))
message("\nsummary:"); print(SUMMARY)
message("outputs: ", DEST)
