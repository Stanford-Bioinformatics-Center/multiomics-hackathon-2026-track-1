#!/usr/bin/env Rscript
# =====================================================================================================
# 08_metabolite_rule_experiments.R — STEP 8: HOW MUCH DOES EACH PART OF THE METABOLITE EDGE RULE MATTER?
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The metabolite network (step 6) connects two metabolites if they share a protein that handles both
#   (Rhea, step 5) AND belong to the same chemical class. This script tests, one change at a time, how
#   many metabolites would be in the network under different rules. It does NOT change the outputs of
#   steps 5-6; it only reports counts. (After these experiments the team chose experiment 2, the super
#   class with our 471 genes, as the step 6 rule.)
#   Rules compared (everything else identical to step 6: shared protein only, no hub removal):
#     baseline      proteins = our 471 genes;                        class = RefMet MAIN class (50)  (original rule)
#     experiment 1  proteins = ALL human reviewed enzymes in Rhea;   class = MAIN class
#     experiment 2  proteins = our 471 genes;                        class = RefMet SUPER class (14) (chosen rule)
#     (side line)   proteins = all Rhea enzymes from ANY organism;   class = MAIN class
#   "Human reviewed enzymes" = UniProtKB/Swiss-Prot entries for Homo sapiens (organism 9606), downloaded
#   from the UniProt REST API. Rhea lists enzymes from every organism, so without this filter a bacterial
#   enzyme could link two human metabolites; the side line shows what that would do.
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis (gene -> UniProt). Internet on first run only
#   (the human Swiss-Prot list is cached in $HACK_EXT next to the Rhea files from step 5).
#
# INPUTS:  Rhea files cached by step 5 ($HACK_EXT); $HACK_OUT/01c_metabolite_ids.csv;
#          $HACK_OUT/01_nodes_EE.csv (our 471 genes)
# OUTPUT:  $HACK_OUT/08_metabolite_rule_experiments.csv   one row per rule: metabolites with a protein,
#          metabolite pairs sharing a protein, edges, metabolites in the network, change vs baseline
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table) })

# Results folder and external-file cache (same settings as step 5).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# External-file cache (Rhea files from step 5; override with HACK_EXT).
EXT <- Sys.getenv("HACK_EXT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/external/rhea"))
# Safety check: step 5 must have downloaded the Rhea files already.
stopifnot(file.exists(file.path(EXT, "rhea-kegg.reaction.gz")), file.exists(file.path(EXT, "rhea2uniprot_sprot.tsv")))

# ---- human reviewed proteins (UniProt), downloaded once ------------------------------------------
# Where the list of human Swiss-Prot accessions is cached.
HUMAN_SP <- file.path(EXT, "uniprot_human_swissprot_accessions.txt")
# Download it if it is not cached yet (one accession per line).
if (!file.exists(HUMAN_SP)) download.file(
  "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=list",
  HUMAN_SP, quiet = TRUE)
# Read the accessions.
human <- readLines(HUMAN_SP)
# Report how many (about 20,000 reviewed human proteins are expected).
message(sprintf("human Swiss-Prot accessions: %d", length(human)))

# ---- Rhea: reactions -> participants and reactions -> enzymes (as in step 5) ----------------------
# Read the reaction file and keep the ENTRY and EQUATION lines.
rx <- readLines(gzfile(file.path(EXT, "rhea-kegg.reaction.gz")))
# Keep only the reaction-start and participant lines.
rx <- rx[grepl("^(ENTRY|EQUATION)", rx)]
# Each EQUATION belongs to the ENTRY above it.
entry_idx <- cumsum(grepl("^ENTRY", rx))
# The reaction number from each ENTRY line.
entries <- sub("^ENTRY\\s+RHEA:", "", rx[grepl("^ENTRY", rx)])
# Each equation line with its reaction number.
eq <- data.table(rhea_id = entries[entry_idx[grepl("^EQUATION", rx)]], line = rx[grepl("^EQUATION", rx)])
# One row per reaction and ChEBI participant.
parts <- unique(eq[, .(chebi = regmatches(line, gregexpr("CHEBI:[0-9]+", line))[[1]]), by = rhea_id])
# One row per reaction and reviewed enzyme (any organism).
enz <- unique(fread(file.path(EXT, "rhea2uniprot_sprot.tsv"), colClasses = "character")[, .(rhea_id = RHEA_ID, uniprot = ID)])

# ---- our metabolites -> ChEBI (original + pH 7.3 form, as in step 5) ------------------------------
# Step 1c's identifier table (also holds the main and super class of every metabolite).
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
# One row per metabolite and ChEBI ID.
mchebi <- ids[chebi_all != "" & !is.na(chebi_all), .(chebi = unlist(strsplit(chebi_all, ";"))), by = metabolite]
# Rhea's pH 7.3 mapping.
ph <- fread(file.path(EXT, "chebi_pH7_3_mapping.tsv"), colClasses = "character")[
  , .(chebi = paste0("CHEBI:", CHEBI), chebi_ph73 = paste0("CHEBI:", CHEBI_PH7_3))]
# Use both the original IDs and their pH 7.3 forms.
mmatch <- unique(rbind(mchebi, merge(mchebi, ph, by = "chebi")[, .(metabolite, chebi = chebi_ph73)]))
# Every metabolite -> enzyme link through a shared Rhea reaction (enzymes from any organism).
m_enz <- unique(merge(merge(mmatch, parts, by = "chebi", allow.cartesian = TRUE), enz, by = "rhea_id",
                      allow.cartesian = TRUE)[, .(metabolite, uniprot)])

# ---- our 471 genes -> UniProt (as in step 5) -------------------------------------------------------
# The 471 gene IDs.
g471 <- fread(file.path(OUT, "01_nodes_EE.csv"), colClasses = list(character = "entrez_gene"))$entrez_gene
# Their UniProt accessions (isoform suffix removed).
ours <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay %in% c("prot-pr", "prot-ol")][
  , .(entrez_gene = as.character(entrez_gene), uniprot = sub("-[0-9]+$", "", as.character(uniprot)))][
  entrez_gene %in% g471 & !is.na(uniprot), uniprot])

# ---- the rule, parameterised ------------------------------------------------------------------------
# Helper: apply the step 6 rule (shared protein AND same class) for a given protein set and class level.
run_rule <- function(protein_set, class_col, label) {
  # metabolite -> protein links restricted to the chosen proteins
  mp <- m_enz[uniprot %in% protein_set]
  # pairs of metabolites that share a protein (each unordered pair once, no self-pairs)
  pr <- unique(merge(mp[, .(m1 = metabolite, uniprot)], mp[, .(m2 = metabolite, uniprot)],
                     by = "uniprot", allow.cartesian = TRUE)[m1 < m2, .(m1, m2)])
  # each metabolite's class at the chosen level (blank = Unclassified, never matches)
  cl <- setNames(ids[[class_col]], ids$metabolite)
  # keep pairs in the same class
  e <- pr[!is.na(cl[m1]) & cl[m1] != "" & cl[m1] == cl[m2]]
  # one row of counts
  data.table(rule = label, proteins_allowed = length(unique(protein_set)),
             metabolites_with_a_protein = uniqueN(mp$metabolite), pairs_sharing_a_protein = nrow(pr),
             edges = nrow(e), metabolites_in_network = uniqueN(c(e$m1, e$m2)))
}

# ---- run the rules ---------------------------------------------------------------------------------
res <- rbind(
  # the current step 6 network
  run_rule(ours, "main_class", "baseline: our 471 genes, main class"),
  # experiment 1: every human reviewed Rhea enzyme may link metabolites
  run_rule(intersect(unique(enz$uniprot), human), "main_class", "exp 1: all human Rhea enzymes, main class"),
  # experiment 2: broader chemical class
  run_rule(ours, "super_class", "exp 2: our 471 genes, super class"),
  # side line: enzymes from any organism (for reference only)
  run_rule(unique(enz$uniprot), "main_class", "side line: Rhea enzymes from any organism, main class"))
# Change in the number of metabolites in the network relative to the baseline.
res[, change_vs_baseline := metabolites_in_network - metabolites_in_network[1]]
# Safety check: the rule matching step 6's class level reproduces step 6 exactly.
s6 <- fread(file.path(OUT, "06_metabolite_summary.csv")); v6 <- setNames(s6$value, s6$metric)
# (row 1 = main class, row 3 = super class, both with our 471 genes)
row6 <- if (v6[["class_level"]] == "main_class") 1 else 3
# (edges and metabolites in the network must match step 6's summary)
stopifnot(res$edges[row6] == as.numeric(v6[["edges_same_class"]]),
          res$metabolites_in_network[row6] == as.numeric(v6[["metabolites_with_edges"]]))
# Save and show.
fwrite(res, file.path(OUT, "08_metabolite_rule_experiments.csv"))
# Show it.
print(res)
