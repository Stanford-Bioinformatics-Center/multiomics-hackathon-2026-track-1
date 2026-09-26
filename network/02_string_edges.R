#!/usr/bin/env Rscript
# =====================================================================================================
# 02_string_edges.R — STEP 2: WHICH GENES ARE CONNECTED (THE EDGES), USING THE STRING DATABASE
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Step 1 made the nodes (471 genes). This step decides which pairs of genes are joined by an edge.
#   An edge exists only if the STRING database says the two genes' proteins interact with high
#   confidence. STRING is used as an on/off GATE: the exercise data plays no part in whether an edge
#   exists, so the endurance and resistance networks have exactly the same edges. (How STRONG each
#   edge is in each arm comes from the exercise data, in step 3.)
#
# WHAT STRING IS
#   STRING (string-db.org) is a public database of known and predicted protein-protein associations.
#   Each pair gets a "combined_score" from 0 to 1000 that merges several kinds of evidence (lab
#   experiments, curated pathway databases, co-expression, text mining, etc.). 700 or more is STRING's
#   "high confidence" level.
#
# THE STRING FILE WE USE (and the filtering already applied to it before it reached us)
#   Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet, a curated file
#   supplied by the team. It contains only human protein pairs with combined_score >= 700 (so the 700
#   cutoff was applied upstream, not here), identified by UniProt protein IDs. It has 124,099 pairs and
#   13,855 proteins, each pair listed once, with no protein paired with itself. A second curated file is
#   expected later; point STRING_PARQUET at it and rerun.
#
# THE RULES WE FOLLOW: El-Kebir et al. 2015, "xHeinz", Bioinformatics 31:3147, section 3.3
#   1. The background network is STRING protein-protein interactions.
#   2. Edges have no direction, no protein is linked to itself, and each pair appears once.
#   3. Remove "outlier hubs": proteins with an extreme number of partners (degree above the 75th
#      percentile + 40 x the interquartile range of all degrees, computed on the FULL network before
#      restricting to our genes). In the paper this removed ubiquitin and ELAVL1, which link to
#      almost everything. In our curated file the cutoff is 741 partners and the busiest protein has 407,
#      so nothing is removed.
#   4. Keep only edges where both proteins belong to our 471 genes (the "induced subnetwork").
#   The paper's network has no edge weights; combined_score is kept only as extra information.
#
# PROTEIN -> GENE
#   STRING speaks UniProt protein IDs; our nodes are Entrez gene IDs. The package's lookup table
#   (HUMAN_FEATURE_TO_GENE) gives the UniProt IDs of each gene's measured proteins. Isoform suffixes
#   ("P12345-2") are removed so they match STRING's canonical IDs. If several protein pairs land on the
#   same gene pair, the one with the highest combined_score is kept.
#
# EXTRA COLUMNS (NOT PART OF EL-KEBIR)
#   cos_EE, cos_RE: the cosine similarity of the two genes' vectors in each arm (how alike their
#   response patterns are, ignoring size). Descriptive only; the edge weights proper are made in step 3.
#
# TECH STACK
#   R 4.4; data.table (tables), nanoparquet (reads the .parquet file), igraph (network components),
#   MotrpacHumanPreSuspensionAnalysis (gene/protein lookup table).
#
# INPUTS:  $STRING_PARQUET (default: the file above in ~/Downloads), $HACK_OUT/01_nodes_{EE,RE}.csv
# OUTPUTS: $HACK_OUT/02_edges.csv          one row per edge
#          $HACK_OUT/02_nodes_string.csv   per gene: in STRING or not, number of partners, component
#          $HACK_OUT/02_network_summary.csv  headline counts
# =====================================================================================================

# Load the packages quietly.
suppressMessages({
  library(MotrpacHumanPreSuspensionAnalysis); library(data.table); library(nanoparquet); library(igraph)
})

# Results folder (override with HACK_OUT); must be the same folder step 1 wrote to.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# The STRING file (override with STRING_PARQUET, e.g. when the second curated file arrives).
STRING_PARQUET <- Sys.getenv("STRING_PARQUET", unset = path.expand(
  "~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet"))
# The "40" in the hub rule (75th percentile + 40 x interquartile range), from El-Kebir 2015.
HUB_IQR_MULT <- 40

# ---- rules 1-2: background network, no direction, no self-links, one row per pair ------------------
# Read the STRING file and keep the two protein IDs (as text) and the score (as a number).
s <- as.data.table(read_parquet(STRING_PARQUET))[, .(a = as.character(protein1), b = as.character(protein2),
                                                     combined_score = as.numeric(combined_score))]
# Remember how many rows the file had.
n_raw <- nrow(s)
# Drop self-links (a == b); write every pair in alphabetical order so A-B and B-A become the same
# pair; if a pair appears more than once, keep its highest score.
s <- s[a != b][, `:=`(p1 = pmin(a, b), p2 = pmax(a, b))][, .(combined_score = max(combined_score)), by = .(p1, p2)]
# Report what the file contained.
message(sprintf("STRING: %d rows -> %d undirected edges, %d proteins, score %g-%g",
                n_raw, nrow(s), uniqueN(c(s$p1, s$p2)), min(s$combined_score), max(s$combined_score)))

# ---- rule 3: remove outlier hubs, judged on the full network --------------------------------------
# Count how many partners ("degree") every protein has.
deg <- table(c(s$p1, s$p2))
# The 25th and 75th percentiles of those counts.
q <- quantile(as.numeric(deg), c(.25, .75), names = FALSE)
# The hub cutoff: 75th percentile + 40 x (75th - 25th percentile).
hub_cut <- q[2] + HUB_IQR_MULT * (q[2] - q[1])
# Proteins above the cutoff are hubs.
hubs <- names(deg)[deg > hub_cut]
# Remove every edge that touches a hub (with this file, none are removed).
s <- s[!(p1 %in% hubs | p2 %in% hubs)]
# Report the rule's numbers and which proteins (if any) were removed.
message(sprintf("hub rule: Q75 %g + %d x IQR %g = %g; max degree %d; %d hubs removed%s",
                q[2], HUB_IQR_MULT, q[2] - q[1], hub_cut, max(deg), length(hubs),
                if (length(hubs)) paste0(" (", paste(hubs, collapse = ","), ")") else ""))

# ---- rule 4: translate proteins to our genes and keep edges inside the 471 --------------------------
# Helper to read one arm's node table from step 1 (gene IDs kept as text, not numbers).
read_nodes <- function(arm) fread(file.path(OUT, sprintf("01_nodes_%s.csv", arm)),
                                  colClasses = list(character = "entrez_gene"))
# Read both arms' node tables.
EE <- read_nodes("EE"); RE <- read_nodes("RE")
# Safety check: 471 genes, in the same order in both files.
stopifnot(nrow(EE) == 471, identical(EE$entrez_gene, RE$entrez_gene))

# Gene -> UniProt lookup, from the package's MS and OLINK protein entries.
acc <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay %in% c("prot-pr", "prot-ol")][
  # keep gene ID and UniProt ID; strip isoform suffixes like "-2" so IDs match STRING's
  , .(entrez_gene = as.character(entrez_gene), uniprot = sub("-[0-9]+$", "", as.character(uniprot)))][
  # only our 471 genes, and only rows that actually have a UniProt ID
  entrez_gene %in% EE$entrez_gene & !is.na(uniprot)])
# Safety check: every one of the 471 genes has at least one UniProt ID.
stopifnot(uniqueN(acc$entrez_gene) == 471)

# Attach a gene to the first protein of each STRING pair (pairs whose protein is not ours drop out).
e <- merge(s, acc[, .(p1 = uniprot, entrez_a = entrez_gene)], by = "p1", allow.cartesian = TRUE)
# Attach a gene to the second protein of each pair (same).
e <- merge(e, acc[, .(p2 = uniprot, entrez_b = entrez_gene)], by = "p2", allow.cartesian = TRUE)
# Drop pairs where both proteins belong to the same gene (a gene is not linked to itself).
e <- e[entrez_a != entrez_b]
# Write each gene pair in a fixed order (smaller gene ID first), carrying each gene's own protein ID
# with it, so the same pair can never appear twice in opposite orders.
e[, flip := entrez_a > entrez_b]
# Swap the flagged rows: g1/g2 are the ordered genes, u1/u2 their own protein IDs.
e[, `:=`(g1 = fifelse(flip, entrez_b, entrez_a), g2 = fifelse(flip, entrez_a, entrez_b),
         u1 = fifelse(flip, p2, p1),             u2 = fifelse(flip, p1, p2))]
# If several protein pairs give the same gene pair, keep the one with the highest STRING score.
e <- e[order(-combined_score)][, .(combined_score = combined_score[1], uniprot_a = u1[1], uniprot_b = u2[1]),
                               by = .(g1, g2)]
# Rename the gene columns.
setnames(e, c("g1", "g2"), c("entrez_a", "entrez_b"))

# ---- extra, descriptive: how alike are the two genes' response patterns in each arm -----------------
# Cosine similarity: +1 = same pattern, 0 = unrelated, -1 = opposite. Dimensions missing in either gene
# are ignored (set to 0 in both).
cosine <- function(M, i, j) {
  # the vectors of the first and second gene of every edge
  x <- M[i, , drop = FALSE]; y <- M[j, , drop = FALSE]
  # ignore dimensions that are missing in either gene
  ok <- !is.na(x) & !is.na(y); x[!ok] <- 0; y[!ok] <- 0
  # cosine = dot product divided by the product of the two vector lengths
  rowSums(x * y) / sqrt(rowSums(x^2) * rowSums(y^2))
}
# The 18 dimension columns (everything except gene ID and symbol).
dim_cols <- names(EE)[-(1:2)]
# Compute the cosine for every edge, once per arm.
for (arm in c("EE", "RE")) {
  # that arm's vectors as a matrix, one row per gene, named by gene ID
  N <- get(arm); M <- as.matrix(N[, ..dim_cols]); rownames(M) <- N$entrez_gene
  # store as cos_EE or cos_RE
  set(e, j = paste0("cos_", arm), value = cosine(M, e$entrez_a, e$entrez_b))
}

# Add gene symbols for readability.
sym <- setNames(EE$gene_symbol, EE$entrez_gene)
# Look up the symbol of each edge's first and second gene.
e[, `:=`(symbol_a = sym[entrez_a], symbol_b = sym[entrez_b])]
# Order the columns and rows (highest STRING score first), then save the edge list.
setcolorder(e, c("entrez_a", "symbol_a", "entrez_b", "symbol_b", "uniprot_a", "uniprot_b",
                 "combined_score", "cos_EE", "cos_RE"))
# Sort rows: highest STRING score first, then alphabetically.
setorder(e, -combined_score, symbol_a, symbol_b)
# Save the edge list.
fwrite(e, file.path(OUT, "02_edges.csv"))

# ---- per-gene table and headline counts -------------------------------------------------------------
# Build the network object (all 471 genes, including those with no edges).
g <- igraph::graph_from_data_frame(e[, .(entrez_a, entrez_b)], directed = FALSE,
                                   vertices = data.table(name = EE$entrez_gene))
# Find connected components (groups of genes that can reach each other through edges).
comp <- igraph::components(g)
# All protein IDs present in STRING (after hub removal).
string_ids <- unique(c(s$p1, s$p2))
# Per gene: its UniProt ID(s), whether STRING knows it, how many partners it has, and its component.
nodes <- data.table(entrez_gene = EE$entrez_gene, gene_symbol = EE$gene_symbol,
                    uniprot = acc[, .(u = paste(sort(unique(uniprot)), collapse = ";")), by = entrez_gene][
                      match(EE$entrez_gene, entrez_gene), u],
                    in_string = EE$entrez_gene %in% acc[uniprot %in% string_ids, entrez_gene],
                    degree = as.integer(igraph::degree(g)[EE$entrez_gene]),
                    component = comp$membership[EE$entrez_gene])
# The size of each gene's component.
nodes[, component_size := comp$csize[component]]
# Save the per-gene table.
fwrite(nodes, file.path(OUT, "02_nodes_string.csv"))

# Headline counts: genes, edges, isolated genes, components, hub cutoff, typical number of partners.
summ <- data.table(
  metric = c("nodes", "nodes_in_string", "nodes_not_in_string", "edges", "isolated_nodes",
             "components_size_ge2", "largest_component", "hub_cutoff", "hubs_removed", "median_degree"),
  value  = c(471, sum(nodes$in_string), sum(!nodes$in_string), nrow(e), sum(nodes$degree == 0),
             sum(comp$csize >= 2), max(comp$csize), hub_cut, length(hubs), median(nodes$degree)))
# Save and show them.
fwrite(summ, file.path(OUT, "02_network_summary.csv"))
# Show the headline counts on screen.
print(summ)
