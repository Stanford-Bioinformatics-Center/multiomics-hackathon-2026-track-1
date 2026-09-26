#!/usr/bin/env Rscript
# =====================================================================================================
# 99_validate_outputs.R — CHECK THAT EVERY OUTPUT OF THE PIPELINE IS INTERNALLY CONSISTENT
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   After running steps 1-4, run this to confirm the results are sound. It does two kinds of checks:
#   1. HARD CHECKS (must always hold; the script stops with an error if one fails): table sizes, no
#      unexpected missing values, no duplicated or self-linked edges, edge weights that really are the dot
#      products of the node vectors, p-values in the valid range, class counts that add up, and so on.
#   2. EXPECTED NUMBERS (the headline results as of 2026-09-26, with the MoTrPAC package v0.2.4 and the
#      first curated STRING file): each is printed as "same" or "CHANGED". A change is not necessarily an
#      error (e.g. a new STRING file will change the edge counts), but it must be explained.
#
# TECH STACK
#   R 4.4; data.table.
#
# INPUT:  every CSV written by steps 1, 1b, 1c, 1d, 2, 3, 4 in $HACK_OUT
# OUTPUT: a report on screen; nothing is written
# =====================================================================================================

# Load data.table quietly.
suppressMessages(library(data.table))

# Results folder (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Helper: read one output table (gene and metabolite IDs kept as text).
rd <- function(f) fread(file.path(OUT, f), colClasses = list(character = intersect(
  c("entrez_gene", "entrez_a", "entrez_b", "pubchem_cid"), names(fread(file.path(OUT, f), nrows = 0)))))
# Counter of passed hard checks.
n_ok <- 0
# Helper: run one hard check; stop with its description if it fails, otherwise count it.
check <- function(ok, what) {
  # stop immediately with a clear message if the check fails
  if (!isTRUE(ok)) stop("FAILED: ", what, call. = FALSE)
  # otherwise count it as passed
  n_ok <<- n_ok + 1
}
# Table of expected headline numbers, filled in as the checks run.
expect <- list()
# Helper: record an observed number next to its expected value.
note <- function(name, observed, expected) expect[[name]] <<- c(observed = observed, expected = expected)

# ---- step 1: gene nodes ------------------------------------------------------------------------------
# The two arms' node tables.
EE <- rd("01_nodes_EE.csv"); RE <- rd("01_nodes_RE.csv")
# 471 genes, 2 ID columns + 18 dimensions, same genes in the same order in both arms.
check(nrow(EE) == 471 && ncol(EE) == 20 && identical(EE$entrez_gene, RE$entrez_gene), "step 1: 471 x 18, same order")
# The dimension columns, and the two that must be empty (adipose protein 0.5 h / 24 h were never sampled).
dims <- names(EE)[-(1:2)]; empty <- c("adipose_prot_0.5h", "adipose_prot_24h")
# Those two are entirely empty in both arms...
check(all(is.na(EE[, ..empty])) && all(is.na(RE[, ..empty])), "step 1: adipose protein 0.5/24 h empty")
# ...and no other value is missing.
check(!anyNA(EE[, setdiff(dims, empty), with = FALSE]) && !anyNA(RE[, setdiff(dims, empty), with = FALSE]),
      "step 1: no other missing values")
# Standard errors are positive wherever a value exists.
SEE <- rd("01_nodes_EE_se.csv"); check(all(as.matrix(SEE[, setdiff(dims, empty), with = FALSE]) > 0), "step 1: SE > 0")
# Correlations between arms are within -1..1.
RHO <- rd("01_nodes_arm_corr.csv"); r <- as.matrix(RHO[, setdiff(dims, empty), with = FALSE])
# Check the correlations are in range.
check(all(r >= -1 & r <= 1), "step 1: EE/RE correlation within -1..1")
# Scaling is exactly raw logFC / scale factor (check one block: muscle RNA 4 h, endurance).
sf <- rd("01_scale_factors.csv"); raw <- rd("01_nodes_EE_raw_logFC.csv")
# Check the scaling.
check(isTRUE(all.equal(EE$muscle_rna_4h, raw$muscle_rna_4h / sf[tissue == "muscle" & ome == "rna", rms_logFC])),
      "step 1: scaled = raw / scale factor")
# Expected headline: median EE/RE correlation of estimates.
note("step1_median_arm_corr", round(median(r), 2), 0.63)

# ---- step 1b: metabolite nodes -----------------------------------------------------------------------
# The two arms' metabolite tables.
ME <- rd("01b_metab_nodes_EE.csv"); MR <- rd("01b_metab_nodes_RE.csv")
# 450 metabolites, 1 name column + 9 dimensions, nothing missing, same order in both arms.
check(nrow(ME) == 450 && ncol(ME) == 10 && !anyNA(ME) && !anyNA(MR) && identical(ME$metabolite, MR$metabolite),
      "step 1b: 450 x 9, complete, same order")
# Metabolite names are unique.
check(!anyDuplicated(ME$metabolite), "step 1b: unique metabolite names")

# ---- step 1c / 1d: metabolite IDs and classes --------------------------------------------------------
# The ID table covers the same 450 metabolites.
ids <- rd("01c_metabolite_ids.csv"); check(setequal(ids$metabolite, ME$metabolite), "step 1c: same 450 metabolites")
# Every ChEBI ID is blank or has the form CHEBI:<digits>.
check(all(ids$chebi_id == "" | is.na(ids$chebi_id) | grepl("^CHEBI:[0-9]+$", ids$chebi_id)), "step 1c: ChEBI format")
# A method is recorded exactly when an ID was found.
has <- !(is.na(ids$chebi_id) | ids$chebi_id == ""); check(all((ids$chebi_method != "none") == has), "step 1c: method <-> ID")
# Expected headline: number with a ChEBI ID.
note("step1c_with_chebi", sum(has), 213)
# Class counts add up to 450 at both levels.
cls <- rd("01d_metabolite_class_counts.csv")
# Check the class counts add up.
check(cls[level == "super_class", sum(n_metabolites)] == 450 && cls[level == "main_class", sum(n_metabolites)] == 450,
      "step 1d: class counts sum to 450")
# Expected headline: number of lipids.
note("step1d_lipids", cls[level == "super_class" & is_lipid == "yes", sum(n_metabolites)], 320)

# ---- step 2: edges -----------------------------------------------------------------------------------
# The edge list.
e <- rd("02_edges.csv")
# No gene is linked to itself.
check(all(e$entrez_a != e$entrez_b), "step 2: no self-links")
# Each gene pair appears once (in either order).
check(!anyDuplicated(e[, .(pmin(entrez_a, entrez_b), pmax(entrez_a, entrez_b))]), "step 2: no duplicate pairs")
# Every endpoint is one of the 471 genes.
check(all(c(e$entrez_a, e$entrez_b) %in% EE$entrez_gene), "step 2: endpoints are nodes")
# Every edge meets the STRING high-confidence threshold.
check(all(e$combined_score >= 700), "step 2: combined_score >= 700")
# Expected headline numbers from the summary.
s2 <- rd("02_network_summary.csv"); v <- setNames(s2$value, s2$metric)
# Record edges and isolated genes.
note("step2_edges", v[["edges"]], 431); note("step2_isolated", v[["isolated_nodes"]], 185)
# Record the largest component and hubs removed.
note("step2_largest_component", v[["largest_component"]], 230); note("step2_hubs_removed", v[["hubs_removed"]], 0)

# ---- step 3: weights ---------------------------------------------------------------------------------
# The weighted edge list, with the same edges as step 2.
w <- rd("03_weighted_edges.csv"); check(nrow(w) == nrow(e), "step 3: same edges as step 2")
# Recompute every weight as the dot product of the two genes' vectors (16 dimensions) and compare.
M <- function(N) { m <- as.matrix(N[, setdiff(dims, empty), with = FALSE]); rownames(m) <- N$entrez_gene; m }
# The two arms' vectors as matrices.
ZE <- M(EE); ZR <- M(RE)
# Check the endurance weights.
check(isTRUE(all.equal(w$w_EE, unname(rowSums(ZE[w$entrez_a, ] * ZE[w$entrez_b, ])))), "step 3: w_EE = dot product")
# Check the resistance weights.
check(isTRUE(all.equal(w$w_RE, unname(rowSums(ZR[w$entrez_a, ] * ZR[w$entrez_b, ])))), "step 3: w_RE = dot product")
# The 0-1 weights are sigmoid(w / s), s = median |w| over both arms.
s <- median(abs(c(w$w_EE, w$w_RE)))
# Check both arms' sigmoid values.
check(isTRUE(all.equal(w$sig_EE, plogis(w$w_EE / s))) && isTRUE(all.equal(w$sig_RE, plogis(w$w_RE / s))), "step 3: sigmoid")
# Expected headline numbers.
note("step3_sigmoid_scale", round(s, 3), 2.667); note("step3_cor_EE_RE", round(cor(w$w_EE, w$w_RE), 2), 0.64)
# Record the sign changes.
note("step3_sign_changes", sum(sign(w$w_EE) != sign(w$w_RE)), 130)

# ---- step 4: comparison ------------------------------------------------------------------------------
# Per-edge and per-gene results.
ed <- rd("04_edge_diff.csv"); nd <- rd("04_node_diff.csv"); s4 <- rd("04_summary.csv"); v4 <- setNames(s4$value, s4$metric)
# Smallest possible bootstrap p-value.
pmin_ <- 2 / (v4[["bootstrap_draws"]] + 1)
# p-values within their valid range, and FDR-adjusted values never below the raw p.
check(all(ed$p_boot >= pmin_ - 1e-12 & ed$p_boot <= 1) && all(ed$fdr >= ed$p_boot - 1e-12), "step 4: edge p/FDR valid")
# Check the gene p-values and FDR.
check(all(nd$p_boot >= pmin_ - 1e-12 & nd$p_boot <= 1) && all(nd$fdr >= nd$p_boot - 1e-12), "step 4: gene p/FDR valid")
# Each 95% interval is ordered (lower end below upper end).
check(all(ed$w_diff_lo <= ed$w_diff_hi), "step 4: intervals ordered")
# Gene strength equals the sum of its edge weights (check endurance); as.numeric turns tapply's
# one-dimensional array into a plain vector so the comparison is like-for-like.
st <- tapply(c(w$w_EE, w$w_EE), c(w$entrez_a, w$entrez_b), sum)
# Check the strengths.
check(isTRUE(all.equal(as.numeric(st[nd$entrez_gene]), nd$strength_EE)), "step 4: strength = sum of edge weights")
# Expected headline numbers.
note("step4_edges_fdr_lt_0.1", v4[["edges_fdr_lt_0.1"]], 0); note("step4_edges_p_lt_0.05", v4[["edges_p_lt_0.05"]], 12)
# Record genes at FDR < 0.1.
note("step4_genes_fdr_lt_0.1", v4[["nodes_fdr_lt_0.1_signed"]], 1)

# ---- steps 5-7: metabolite-protein links (Rhea), metabolite networks, hub report ---------------------
# Metabolite-protein links: every metabolite is one of the 450 and every gene one of the 471.
lk <- rd("05_metabolite_protein_links.csv")
# Check the links only use our metabolites and genes.
check(all(lk$metabolite %in% ME$metabolite) && all(lk$entrez_gene %in% EE$entrez_gene), "step 5: links use our nodes")
# Metabolite edges: both ends are metabolites, no self-links, no duplicate pairs.
me <- rd("06_metabolite_edges.csv")
# Check the metabolite edge list is clean.
check(all(c(me$metabolite_a, me$metabolite_b) %in% ME$metabolite) && all(me$metabolite_a != me$metabolite_b) &&
      !anyDuplicated(me[, .(pmin(metabolite_a, metabolite_b), pmax(metabolite_a, metabolite_b))]), "step 6: clean edge list")
# The class rule holds: both ends of every edge have the edge's class, at the level step 6 used.
ids2 <- rd("01c_metabolite_ids.csv"); cl <- setNames(ids2[[me$class_level[1]]], ids2$metabolite)
# Check both ends of every edge are in the edge's class.
check(all(cl[me$metabolite_a] == me$class & cl[me$metabolite_b] == me$class), "step 6: same class")
# The protein rule holds: every listed shared protein handles both metabolites in step 5.
sp <- me[, .(g = unlist(strsplit(shared_proteins, ";"))), by = .(metabolite_a, metabolite_b)]
# Every metabolite-protein pair from step 5, as text keys.
key <- paste(lk$metabolite, lk$gene_symbol)
# Check each listed shared protein handles both metabolites.
check(all(paste(sp$metabolite_a, sp$g) %in% key & paste(sp$metabolite_b, sp$g) %in% key), "step 6: shared protein handles both")
# Metabolite edge weights are the dot products of the 9-number vectors.
MEm <- as.matrix(ME[, -1]); rownames(MEm) <- ME$metabolite; MRm <- as.matrix(MR[, -1]); rownames(MRm) <- MR$metabolite
# Check both arms' weights.
check(isTRUE(all.equal(me$w_EE, unname(rowSums(MEm[me$metabolite_a, ] * MEm[me$metabolite_b, ])))) &&
      isTRUE(all.equal(me$w_RE, unname(rowSums(MRm[me$metabolite_a, ] * MRm[me$metabolite_b, ])))), "step 6: weights = dot products")
# Record the headline numbers for steps 5-7.
s5 <- rd("05_rhea_summary.csv"); v5 <- setNames(s5$value, s5$metric)
# Record metabolites linked to our genes.
note("step5_metabolites_linked", as.numeric(v5[["metabolites_linked_to_our_genes"]]), 60)
# Record genes linked to our metabolites.
note("step5_genes_linked", as.numeric(v5[["genes_linked_to_our_metabolites"]]), 80)
# Record the number of metabolite edges.
note("step6_metabolite_edges", nrow(me), 122)
# Record the number of metabolites in the metabolite network.
note("step6_metabolites_in_network", uniqueN(c(me$metabolite_a, me$metabolite_b)), 44)
# Record gene hubs above the Tukey fence.
hs <- rd("07_hub_summary.csv"); note("step7_gene_hubs_tukey", hs[1, n_hubs_tukey], 13)
# Record hubs flagged by the El-Kebir rule in any network.
note("step7_hubs_elkebir_all_networks", sum(hs$n_hubs_elkebir), 0)

# ---- report ------------------------------------------------------------------------------------------
# All hard checks passed if we got here.
cat(sprintf("\n%d hard checks passed.\n\nExpected headline numbers (as of 2026-09-26):\n", n_ok))
# One line per headline number: observed, expected, and whether it changed.
for (k in names(expect)) cat(sprintf("  %-26s observed %-8s expected %-8s %s\n", k, expect[[k]][["observed"]],
                                     expect[[k]][["expected"]],
                                     if (expect[[k]][["observed"]] == expect[[k]][["expected"]]) "same" else "CHANGED"))
