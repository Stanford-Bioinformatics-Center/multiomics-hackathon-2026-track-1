#!/usr/bin/env Rscript
# 02_string_edges.R — STRING edges between the 471 nodes, built with the network rules of
# El-Kebir et al. 2015 (xHeinz, Bioinformatics 31:3147), section 3.3:
#   1. background = STRING protein-protein network (here: the curated combined_score >= 700 file)
#   2. undirected, no self-loops, one edge per pair
#   3. outlier hubs removed: degree > Q75 + 40 x IQR of the node-degree distribution,
#      computed on the FULL background network, before restricting to measured genes
#   4. network induced on the genes that pass the data filters (our 471 nodes)
# The paper's network is unweighted; combined_score is carried along as an attribute only.
#
# Protein -> gene: UniProt accessions from HUMAN_FEATURE_TO_GENE (prot-pr + prot-ol rows, isoform
# suffix stripped). If several protein pairs map to one gene pair, the highest combined_score is kept.
#
# Not part of El-Kebir: per-arm edge attributes cos_EE / cos_RE = cosine similarity of the two
# endpoints' 18-dim embeddings (over dimensions present in both; the adipose protein 0.5 h / 24 h
# dimensions are always NA). The topology is the same for both arms — STRING does not depend
# on the exercise data — so arm differences live in the node vectors and these edge attributes.
#
# Inputs:  $STRING_PARQUET (default ~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet)
#          $HACK_OUT/01_nodes_{EE,RE}.csv from 01_node_embeddings.R
# Outputs: $HACK_OUT/02_edges.csv, 02_nodes_string.csv, 02_network_summary.csv

suppressMessages({
  library(MotrpacHumanPreSuspensionAnalysis); library(data.table); library(nanoparquet); library(igraph)
})

OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
STRING_PARQUET <- Sys.getenv("STRING_PARQUET", unset = path.expand(
  "~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet"))
HUB_IQR_MULT <- 40

# ---- 1-2. background network: undirected, no self-loops, one edge per pair ----------------
s <- as.data.table(read_parquet(STRING_PARQUET))[, .(a = as.character(protein1), b = as.character(protein2),
                                                     combined_score = as.numeric(combined_score))]
n_raw <- nrow(s)
s <- s[a != b][, `:=`(p1 = pmin(a, b), p2 = pmax(a, b))][, .(combined_score = max(combined_score)), by = .(p1, p2)]
message(sprintf("STRING: %d rows -> %d undirected edges, %d proteins, score %g-%g",
                n_raw, nrow(s), uniqueN(c(s$p1, s$p2)), min(s$combined_score), max(s$combined_score)))

# ---- 3. hub removal on the full background ----------------------------------------------
deg <- table(c(s$p1, s$p2))
q <- quantile(as.numeric(deg), c(.25, .75), names = FALSE)
hub_cut <- q[2] + HUB_IQR_MULT * (q[2] - q[1])
hubs <- names(deg)[deg > hub_cut]
s <- s[!(p1 %in% hubs | p2 %in% hubs)]
message(sprintf("hub rule: Q75 %g + %d x IQR %g = %g; max degree %d; %d hubs removed%s",
                q[2], HUB_IQR_MULT, q[2] - q[1], hub_cut, max(deg), length(hubs),
                if (length(hubs)) paste0(" (", paste(hubs, collapse = ","), ")") else ""))

# ---- 4. map to the 471 genes and induce ----------------------------------------------------
read_nodes <- function(arm) fread(file.path(OUT, sprintf("01_nodes_%s.csv", arm)),
                                  colClasses = list(character = "entrez_gene"))
EE <- read_nodes("EE"); RE <- read_nodes("RE")
stopifnot(nrow(EE) == 471, identical(EE$entrez_gene, RE$entrez_gene))

acc <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay %in% c("prot-pr", "prot-ol")][
  , .(entrez_gene = as.character(entrez_gene), uniprot = sub("-[0-9]+$", "", as.character(uniprot)))][
  entrez_gene %in% EE$entrez_gene & !is.na(uniprot)])
stopifnot(uniqueN(acc$entrez_gene) == 471)

e <- merge(s, acc[, .(p1 = uniprot, entrez_a = entrez_gene)], by = "p1", allow.cartesian = TRUE)
e <- merge(e, acc[, .(p2 = uniprot, entrez_b = entrez_gene)], by = "p2", allow.cartesian = TRUE)
e <- e[entrez_a != entrez_b]
# orient each row so gene 1 < gene 2, carrying each gene's own accession with it
e[, flip := entrez_a > entrez_b]
e[, `:=`(g1 = fifelse(flip, entrez_b, entrez_a), g2 = fifelse(flip, entrez_a, entrez_b),
         u1 = fifelse(flip, p2, p1),             u2 = fifelse(flip, p1, p2))]
e <- e[order(-combined_score)][, .(combined_score = combined_score[1], uniprot_a = u1[1], uniprot_b = u2[1]),
                               by = .(g1, g2)]
setnames(e, c("g1", "g2"), c("entrez_a", "entrez_b"))

# ---- per-arm embedding similarity on each edge (not part of El-Kebir) ----------------------
cosine <- function(M, i, j) {
  x <- M[i, , drop = FALSE]; y <- M[j, , drop = FALSE]
  ok <- !is.na(x) & !is.na(y); x[!ok] <- 0; y[!ok] <- 0
  rowSums(x * y) / sqrt(rowSums(x^2) * rowSums(y^2))
}
dim_cols <- names(EE)[-(1:2)]
for (arm in c("EE", "RE")) {
  N <- get(arm); M <- as.matrix(N[, ..dim_cols]); rownames(M) <- N$entrez_gene
  set(e, j = paste0("cos_", arm), value = cosine(M, e$entrez_a, e$entrez_b))
}

sym <- setNames(EE$gene_symbol, EE$entrez_gene)
e[, `:=`(symbol_a = sym[entrez_a], symbol_b = sym[entrez_b])]
setcolorder(e, c("entrez_a", "symbol_a", "entrez_b", "symbol_b", "uniprot_a", "uniprot_b",
                 "combined_score", "cos_EE", "cos_RE"))
setorder(e, -combined_score, symbol_a, symbol_b)
fwrite(e, file.path(OUT, "02_edges.csv"))

# ---- node table + summary --------------------------------------------------------------
g <- igraph::graph_from_data_frame(e[, .(entrez_a, entrez_b)], directed = FALSE,
                                   vertices = data.table(name = EE$entrez_gene))
comp <- igraph::components(g)
string_ids <- unique(c(s$p1, s$p2))
nodes <- data.table(entrez_gene = EE$entrez_gene, gene_symbol = EE$gene_symbol,
                    uniprot = acc[, .(u = paste(sort(unique(uniprot)), collapse = ";")), by = entrez_gene][
                      match(EE$entrez_gene, entrez_gene), u],
                    in_string = EE$entrez_gene %in% acc[uniprot %in% string_ids, entrez_gene],
                    degree = as.integer(igraph::degree(g)[EE$entrez_gene]),
                    component = comp$membership[EE$entrez_gene])
nodes[, component_size := comp$csize[component]]
fwrite(nodes, file.path(OUT, "02_nodes_string.csv"))

summ <- data.table(
  metric = c("nodes", "nodes_in_string", "nodes_not_in_string", "edges", "isolated_nodes",
             "components_size_ge2", "largest_component", "hub_cutoff", "hubs_removed", "median_degree"),
  value  = c(471, sum(nodes$in_string), sum(!nodes$in_string), nrow(e), sum(nodes$degree == 0),
             sum(comp$csize >= 2), max(comp$csize), hub_cut, length(hubs), median(nodes$degree)))
fwrite(summ, file.path(OUT, "02_network_summary.csv"))
print(summ)
