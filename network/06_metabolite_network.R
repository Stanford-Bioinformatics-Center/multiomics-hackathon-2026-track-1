#!/usr/bin/env Rscript
# =====================================================================================================
# 06_metabolite_network.R — STEP 6: THE METABOLITE-ONLY NETWORKS (ENDURANCE AND RESISTANCE)
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   Builds two metabolite-to-metabolite networks, one per exercise arm, with the SAME edges and different
#   edge weights, exactly like the gene networks (steps 2-3).
#   EDGE RULE (the team's rule, applied as an on/off gate that never uses the exercise data):
#     two metabolites are connected if BOTH
#       (a) they are handled by at least one SAME protein, where the protein must be one of our 471 genes
#           and "handled" means Rhea records the metabolite in a reaction that protein catalyses (step 5);
#       (b) they belong to the same RefMet SUPER class (14 broad families, e.g. both "Nucleic acids",
#           both "Fatty Acyls"; from step 1c). The team chose the super class over the 50 main classes
#           after step 8 showed the main-class rule left only 39 metabolites connected (super class: 44).
#           Set METAB_CLASS_LEVEL=main_class to use the narrower classes instead.
#     Only shared proteins count (no second step through interacting proteins).
#   EDGE WEIGHT (per arm, as in step 3): the dot product of the two metabolites' 9-number vectors from
#   step 1b (adipose, blood, muscle x 0.5 / 4 / 24 h), plus the 0-1 version sigmoid(w / s) with
#   s = median |w| over both arms.
#   These are INFERRED FUNCTIONAL links (shared enzyme + same chemical class), not physical interactions:
#   metabolites do not bind each other.
#
# NO HUB REMOVAL (the team's decision for now). A protein that handles many metabolites of one class
#   connects all of them to each other; step 7 reports these hubs.
#
# TECH STACK
#   R 4.4; data.table; igraph (connected components).
#
# INPUTS:  $HACK_OUT/05_metabolite_protein_links.csv, $HACK_OUT/01c_metabolite_ids.csv (classes),
#          $HACK_OUT/01b_metab_nodes_EE.csv and _RE.csv (the 9-number vectors)
# OUTPUTS: $HACK_OUT/06_metabolite_edges.csv    one row per edge: shared proteins, w_EE, w_RE, w_diff, sig_*
#          $HACK_OUT/06_metabolite_nodes.csv    one row per metabolite: class, proteins, degree, component
#          $HACK_OUT/06_metabolite_summary.csv  headline counts
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(igraph) })

# Results folder (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))

# Metabolite -> protein links from step 5.
link <- fread(file.path(OUT, "05_metabolite_protein_links.csv"), colClasses = list(character = "entrez_gene"))
# Which class level must match: the super class (default, the team's choice) or the main class.
CLASS_LEVEL <- Sys.getenv("METAB_CLASS_LEVEL", unset = "super_class")
# Safety check: only these two levels exist.
stopifnot(CLASS_LEVEL %in% c("super_class", "main_class"))
# Every metabolite's classes (step 1c); a blank class means "Unclassified", which never forms edges.
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))[, .(metabolite, super_class, main_class)]
# The class used by the rule, at the chosen level.
ids[, rule_class := get(CLASS_LEVEL)]

# ---- edge rule ------------------------------------------------------------------------------------
# All pairs of metabolites that share a protein: join the link table to itself on the protein.
pairs <- merge(link[, .(m1 = metabolite, entrez_gene, gene_symbol)],
               link[, .(m2 = metabolite, entrez_gene)], by = "entrez_gene", allow.cartesian = TRUE)
# Keep each unordered pair once (alphabetical order) and drop a metabolite paired with itself.
pairs <- pairs[m1 < m2]
# Per pair: how many proteins they share, and which.
pairs <- pairs[, .(n_shared_proteins = uniqueN(entrez_gene),
                   shared_proteins = paste(sort(unique(gene_symbol)), collapse = ";")), by = .(m1, m2)]
# Attach each metabolite's rule class (and its main class, for information).
pairs <- merge(pairs, ids[, .(m1 = metabolite, class1 = rule_class, main_class_a = main_class)], by = "m1")
# ...and of the second metabolite.
pairs <- merge(pairs, ids[, .(m2 = metabolite, class2 = rule_class, main_class_b = main_class)], by = "m2")
# Keep pairs in the same rule class (a missing class never matches).
e <- pairs[!is.na(class1) & class1 != "" & class1 == class2]
# Report how many shared-protein pairs the class rule keeps.
message(sprintf("pairs sharing a protein: %d; also same %s (edges): %d", nrow(pairs), CLASS_LEVEL, nrow(e)))

# ---- edge weights per arm (dot products of the 9-number vectors) --------------------------------
# Helper: read one arm's metabolite vectors as a matrix (row names = metabolite names).
vec <- function(arm) {
  x <- fread(file.path(OUT, sprintf("01b_metab_nodes_%s.csv", arm)))
  m <- as.matrix(x[, -1]); rownames(m) <- x$metabolite
  m
}
# Both arms.
ME <- vec("EE"); MR <- vec("RE")
# Safety check: same metabolites and dimensions in both arms, nothing missing.
stopifnot(identical(dimnames(ME), dimnames(MR)), !anyNA(ME), !anyNA(MR))
# The dot product of the two metabolites' vectors, per arm.
e[, w_EE := rowSums(ME[m1, , drop = FALSE] * ME[m2, , drop = FALSE])]
# The same for the resistance arm.
e[, w_RE := rowSums(MR[m1, , drop = FALSE] * MR[m2, , drop = FALSE])]
# The difference between arms (positive = the endurance weight is higher).
e[, w_diff := w_EE - w_RE]
# The 0-1 version: sigmoid(w / s), one s (median |w| over both arms) so the arms share a unit.
SIG_SCALE <- if (nrow(e)) median(abs(c(e$w_EE, e$w_RE))) else NA_real_
# Apply the sigmoid to both arms' weights.
e[, `:=`(sig_EE = plogis(w_EE / SIG_SCALE), sig_RE = plogis(w_RE / SIG_SCALE))]
# Tidy column order and names, strongest shared evidence first, then save.
setnames(e, c("m1", "m2", "class1"), c("metabolite_a", "metabolite_b", "class"))
# Record which class level the rule used.
e[, class_level := CLASS_LEVEL]
# The second rule-class column is now redundant (equal to class): drop it.
e[, class2 := NULL]
# Put the columns in a readable order.
setcolorder(e, c("metabolite_a", "metabolite_b", "class", "class_level", "main_class_a", "main_class_b",
                 "n_shared_proteins", "shared_proteins", "w_EE", "w_RE", "w_diff", "sig_EE", "sig_RE"))
# Sort: by class, most shared proteins first, then alphabetically.
setorder(e, class, -n_shared_proteins, metabolite_a, metabolite_b)
# Save the edge list.
fwrite(e, file.path(OUT, "06_metabolite_edges.csv"))

# ---- per-metabolite table and headline counts ----------------------------------------------------
# The network over all 450 metabolites (most will have no edge).
g <- graph_from_data_frame(e[, .(metabolite_a, metabolite_b)], directed = FALSE,
                           vertices = data.table(name = rownames(ME)))
# Connected components.
comp <- components(g)
# Per metabolite: class, how many of our proteins handle it, degree, component and its size.
nodes <- data.table(metabolite = rownames(ME))
# Add each metabolite's classes.
nodes <- ids[nodes, on = "metabolite"]
# Count how many of our proteins handle each metabolite (from step 5).
nodes[, n_proteins := link[, uniqueN(entrez_gene), by = metabolite][match(nodes$metabolite, metabolite), V1]]
# Metabolites with no protein get 0 instead of missing.
nodes[is.na(n_proteins), n_proteins := 0L]
# Number of metabolite neighbours and which connected group each metabolite is in.
nodes[, `:=`(degree = as.integer(degree(g)[metabolite]), component = comp$membership[metabolite])]
# The size of that group.
nodes[, component_size := comp$csize[component]]
# Save the per-metabolite table.
fwrite(nodes, file.path(OUT, "06_metabolite_nodes.csv"))

# Headline counts.
summ <- data.table(
  metric = c("class_level", "metabolites", "metabolites_with_a_protein", "pairs_sharing_a_protein", "edges_same_class",
             "metabolites_with_edges", "components_size_ge2", "largest_component", "median_degree_nonisolated",
             "sigmoid_scale", "cor_w_EE_w_RE", "edges_sign_change"),
  value = c(CLASS_LEVEL, nrow(nodes), sum(nodes$n_proteins > 0), nrow(pairs), nrow(e),
            sum(nodes$degree > 0), sum(comp$csize >= 2), max(comp$csize), median(nodes$degree[nodes$degree > 0]),
            SIG_SCALE, if (nrow(e) > 2) cor(e$w_EE, e$w_RE) else NA, sum(sign(e$w_EE) != sign(e$w_RE))))
# Save and show.
fwrite(summ, file.path(OUT, "06_metabolite_summary.csv"))
# Show the summary on screen.
print(summ)
