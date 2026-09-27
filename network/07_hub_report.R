#!/usr/bin/env Rscript
# =====================================================================================================
# 07_hub_report.R — STEP 7: HOW MANY HUBS ARE THERE, AND HOW MANY ANALYTES HANG ON EACH? (NOTHING REMOVED)
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   Are there nodes that dominate our networks (hubs), which are they, and how much hangs on each? This is
#   a report to inform the team's decision on hub removal; nothing is removed.
#
# HUB DEFINITION (team decision, 2026-09-26): the HUB is the node with the highest summed edge weight
#   ("strength" = sum of |w| over its edges, per arm; the weights are signed dot products, so their sizes are
#   summed). "Strength hubs" are all nodes whose strength (the larger of the two arms) is above the Tukey fence
#   of the strength distribution (75th percentile + 1.5 x interquartile range). The degree-based rules below are
#   kept for reference (the El-Kebir rule is the one used to prune edges in step 2).
#
# WHAT THIS SCRIPT DOES (plain language)
#   A degree-"hub" is a node with far more connections than the rest. Hubs can dominate a network (one hub can
#   tie together everything around it), so the team wants to see them before deciding whether to remove
#   any. This script only REPORTS hubs; it removes nothing.
#   It covers the gene, metabolite and joint networks:
#     - gene networks, endurance (EE) and resistance (RE): same edges (STRING, step 2), different weights;
#     - metabolite networks, EE and RE: same edges (Rhea + class rule, step 6), different weights;
#     - the joint protein-metabolite network (step 14), so run_all.sh runs this step after step 14.
#   Because the two arms share their edges, the number of connections (degree) of every node is the same
#   in both arms; what differs by arm is how STRONG each node's connections are ("strength" = sum of the
#   sizes of its edge weights), which is reported per arm.
#   For the metabolite networks there is a second kind of hub: a PROTEIN that handles many metabolites
#   (Rhea, step 5). Such a protein links all its same-class metabolites to each other, so it is reported
#   too, with the number of metabolites it handles and the number of metabolite edges it supports.
#
# TWO HUB THRESHOLDS (both computed on the degree distribution of the network in question)
#   - El-Kebir rule (the one used in step 2): degree > 75th percentile + 40 x interquartile range.
#     Designed to catch only extreme outliers such as ubiquitin in the full human interactome.
#   - Tukey outlier fence: degree > 75th percentile + 1.5 x interquartile range, the standard boxplot
#     definition of an outlier. Far more lenient; it shows the MOST hubs anyone could reasonably flag.
#   Nodes with no connections are left out of the distribution (they would pull the percentiles to zero).
#
# TECH STACK
#   R 4.4; data.table.
#
# INPUTS:  $HACK_OUT/03_weighted_edges.csv, 06_metabolite_edges.csv, 05_metabolite_protein_links.csv, 14_joint_edges.csv
# OUTPUTS: $HACK_OUT/07_hub_summary.csv   one row per network x hub type: THE HUB (highest strength, per arm and
#                                         overall), strength-hub count, and the degree cutoffs / counts for reference
#          $HACK_OUT/07_hub_list.csv      every node above the strength fence or the degree Tukey fence, with its
#                                         strength per arm, degree, and what is attached to it
#
# HOW TO RUN
#   After steps 3, 5, 6 and 14:   Rscript network/07_hub_report.R   (run_all.sh runs it right after step 14)
#
# KNOWN LIMITS
#   - Gene nodes are named by gene symbol, so two genes sharing a symbol would be counted as one node.
#   - For mediating proteins, "strength" sums only the metabolite edges where the protein is a SHARED
#     protein; edges it supports only through a STRING-interacting partner are not counted.
#   - Both fences are conventions (El-Kebir 40 x IQR, Tukey 1.5 x IQR), not statistical tests.
# =====================================================================================================

# Load data.table quietly.
suppressMessages(library(data.table))

# Results folder (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))

# Helper: given a table of nodes with their degree, attached names and per-arm strength, compute both
# hub thresholds and return (a) a one-row summary and (b) the list of nodes above the Tukey fence.
hubs <- function(nodes, network, hub_type) {
  # the degrees of connected nodes only
  d <- nodes[degree > 0, degree]
  # 25th and 75th percentiles and the interquartile range
  q <- quantile(d, c(.25, .75), names = FALSE); iqr <- q[2] - q[1]
  # the two cutoffs
  cut_ek <- q[2] + 40 * iqr; cut_tk <- q[2] + 1.5 * iqr
  # the node with the most connections (ties: first alphabetically)
  top <- nodes[order(-degree, node)][1]
  # THE HUB (team definition): highest summed edge weight; per arm, and overall = the larger of the two arms
  nodes[, strength_max := pmax(strength_EE, strength_RE)]
  # the strength fence: over connected nodes, 75th percentile + 1.5 x interquartile range of strength_max
  sm <- nodes[degree > 0, strength_max]; qs <- quantile(sm, c(.25, .75), names = FALSE); cut_st <- qs[2] + 1.5 * (qs[2] - qs[1])
  # THE HUB itself: the node with the largest strength_max (ties: first alphabetically)
  hub <- nodes[order(-strength_max, node)][1]
  # the one-row summary
  s <- data.table(network = network, hub_type = hub_type, connected_nodes = length(d),
                  hub = hub$node, hub_strength_EE = round(hub$strength_EE, 3), hub_strength_RE = round(hub$strength_RE, 3), hub_degree = hub$degree,
                  hub_EE = nodes[order(-strength_EE, node)][1, node], hub_RE = nodes[order(-strength_RE, node)][1, node],
                  strength_cutoff = round(cut_st, 3), n_strength_hubs = sum(sm > cut_st),
                  median_degree = median(d), q75 = q[2], iqr = iqr,
                  elkebir_cutoff = cut_ek, n_hubs_elkebir = sum(d > cut_ek),
                  tukey_cutoff = cut_tk, n_hubs_tukey = sum(d > cut_tk),
                  top_hub = top$node, top_hub_degree = top$degree,
                  top_hub_attached = top$attached,
                  top_strength_EE = nodes[order(-strength_EE)][1, paste0(node, " (", round(strength_EE, 1), ")")],
                  top_strength_RE = nodes[order(-strength_RE)][1, paste0(node, " (", round(strength_RE, 1), ")")])
  # every node above the strength fence or the degree Tukey fence, strongest first
  l <- nodes[strength_max > cut_st | degree > cut_tk][order(-strength_max)][, .(network = network, hub_type = hub_type, node, degree,
                                                    strength_hub = strength_max > cut_st, is_the_hub = node == hub$node,
                                                    above_degree_tukey = degree > cut_tk, above_elkebir = degree > cut_ek,
                                                    strength_EE = round(strength_EE, 2), strength_RE = round(strength_RE, 2),
                                                    attached)]
  # return both
  list(summary = s, list = l)
}

# Helper: per-node degree, attached neighbours and per-arm strength from an edge list (columns a, b, w_EE, w_RE).
node_table <- function(e) {
  # every edge listed once from each end
  both <- rbind(e[, .(node = a, other = b, w_EE, w_RE)], e[, .(node = b, other = a, w_EE, w_RE)])
  # per node: number of neighbours, their names, and the summed size of its edge weights in each arm
  both[, .(degree = uniqueN(other), attached = paste(sort(unique(other)), collapse = ";"),
           strength_EE = sum(abs(w_EE)), strength_RE = sum(abs(w_RE))), by = node]
}

# ---- gene networks (EE and RE share edges) --------------------------------------------------------
# The weighted gene edges from step 3, using gene symbols as node names.
ge <- fread(file.path(OUT, "03_weighted_edges.csv"))[, .(a = symbol_a, b = symbol_b, w_EE, w_RE)]
# Hub report for genes.
r1 <- hubs(node_table(ge), "gene networks (EE & RE)", "gene (node degree)")

# ---- metabolite networks (EE and RE share edges) --------------------------------------------------
# The weighted metabolite edges from step 6.
me <- fread(file.path(OUT, "06_metabolite_edges.csv"))
# Hub report for metabolite nodes.
r2 <- hubs(node_table(me[, .(a = metabolite_a, b = metabolite_b, w_EE, w_RE)]),
           "metabolite networks (EE & RE)", "metabolite (node degree)")

# ---- mediating proteins behind the metabolite networks -------------------------------------------
# Metabolite -> protein links (step 5).
lk <- fread(file.path(OUT, "05_metabolite_protein_links.csv"))
# Per protein: how many of our metabolites it handles, and which ones.
pn <- lk[, .(degree = uniqueN(metabolite), attached = paste(sort(unique(metabolite)), collapse = ";")),
         by = .(node = gene_symbol)]
# Per protein: how many metabolite edges list it as a shared protein, and the summed edge sizes per arm.
sp <- me[, .(node = unlist(strsplit(shared_proteins, ";"))), by = .(metabolite_a, metabolite_b, w_EE, w_RE)]
# Per protein: number of metabolite edges it supports and their summed sizes in each arm.
ps <- sp[, .(edges_supported = .N, strength_EE = sum(abs(w_EE)), strength_RE = sum(abs(w_RE))), by = node]
# Combine (proteins that support no edge get zero strength).
pn <- ps[pn, on = "node"]
# Proteins that support no edge get zeros instead of missing values.
pn[is.na(edges_supported), `:=`(edges_supported = 0L, strength_EE = 0, strength_RE = 0)]
# Hub report for proteins (degree = number of metabolites handled).
r3 <- hubs(pn, "metabolite networks (EE & RE)", "mediating protein (metabolites handled)")
# Add how many metabolite edges each listed protein supports.
r3$list <- pn[, .(node, edges_supported)][r3$list, on = "node"]
# And in the summary, how many edges the top protein supports.
r3$summary[, top_hub_edges_supported := pn[node == r3$summary$top_hub, edges_supported]]

# ---- joint protein-metabolite network (step 14) -----------------------------------------------------
# The joint network's edges from step 14 (protein-protein, metabolite-metabolite and protein-metabolite).
je <- fread(file.path(OUT, "14_joint_edges.csv"))[, .(a = node_a, b = node_b, w_EE, w_RE)]
# Hub report for the joint network (proteins and metabolites together).
r4 <- hubs(node_table(je), "joint protein-metabolite network (EE & RE)", "protein or metabolite (joint)")

# ---- write --------------------------------------------------------------------------------------
# One summary table for all hub types.
summ <- rbind(r1$summary, r2$summary, r3$summary, r4$summary, fill = TRUE)
# Save it.
fwrite(summ, file.path(OUT, "07_hub_summary.csv"))
# One list of all flagged nodes.
fwrite(rbind(r1$list, r2$list, r3$list, r4$list, fill = TRUE), file.path(OUT, "07_hub_list.csv"))
# Show the summary without the long "attached" column.
print(summ[, !"top_hub_attached"])
