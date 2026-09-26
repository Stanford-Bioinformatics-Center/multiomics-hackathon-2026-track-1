#!/usr/bin/env Rscript
# =====================================================================================================
# 08_metabolite_rule_experiments.R — STEP 8: HOW MUCH DOES EACH PART OF THE METABOLITE EDGE RULE MATTER?
# =====================================================================================================
#
# WHAT THIS SCRIPT DOES (plain language)
#   The metabolite network (step 6) connects two metabolites if they share a protein that handles both
#   (Rhea, step 5) AND belong to the same chemical class. This script tests, one change at a time, how
#   many metabolites would be in the network under different rules. It does NOT change the outputs of
#   steps 5-6; it only reports counts. (After these experiments the team chose experiment 3 as the step 6
#   rule: super class, our 471 genes, and shared OR STRING-interacting proteins.)
#   Rules compared (everything else identical to step 6: shared protein only, no hub removal):
#     baseline      proteins = our 471 genes;                        class = RefMet MAIN class (50)  (original rule)
#     experiment 1  proteins = ALL human reviewed enzymes in Rhea;   class = MAIN class
#     experiment 2  proteins = our 471 genes;                        class = RefMet SUPER class (14)
#     (side line)   proteins = all Rhea enzymes from ANY organism;   class = MAIN class
#     experiment 3  proteins = our 471 genes;                        class = SUPER class; two metabolites
#                   are also linked if their proteins are DIFFERENT but interact in STRING (>= 700),
#                   i.e. "same protein OR STRING-interacting proteins" instead of "same protein" only (chosen rule).
#   For every rule the table also reports how the two arms compare on that network: the correlation
#   between the arms' edge weights (dot products of the step 1b vectors) and how many edges change sign.
#   Terminology: none of these is a physical interaction BETWEEN metabolites. The metabolite-protein link
#   is enzyme-substrate (Rhea); the protein-protein link is STRING association from all evidence types.
#   "Human reviewed enzymes" = UniProtKB/Swiss-Prot entries for Homo sapiens (organism 9606), downloaded
#   from the UniProt REST API. Rhea lists enzymes from every organism, so without this filter a bacterial
#   enzyme could link two human metabolites; the side line shows what that would do.
#
# TECH STACK
#   R 4.4; data.table; MotrpacHumanPreSuspensionAnalysis (gene -> UniProt). Internet on first run only
#   (the human Swiss-Prot list is cached in $HACK_EXT next to the Rhea files from step 5).
#
# INPUTS:  Rhea files cached by step 5 ($HACK_EXT); $HACK_OUT/01c_metabolite_ids.csv;
#          $HACK_OUT/01b_metab_nodes_{EE,RE}.csv (vectors); $STRING_PARQUET (the step 2 STRING file);
#          $HACK_OUT/01_nodes_EE.csv (our 471 genes)
# OUTPUT:  $HACK_OUT/08_string_neighbour_extra_edges.csv   the edges experiment 3 adds to experiment 2, with
#          the protein pair(s) that link them and both arms' weights
#          $HACK_OUT/08_metabolite_rule_experiments.csv   one row per rule: metabolites with a protein,
#          metabolite pairs sharing a protein, edges, metabolites in the network, change vs baseline
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(MotrpacHumanPreSuspensionAnalysis); library(data.table); library(nanoparquet) })

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

# ---- STRING protein pairs and the metabolite vectors ------------------------------------------------
# The same STRING file as step 2 (override with STRING_PARQUET).
STRING_PARQUET <- Sys.getenv("STRING_PARQUET", unset = path.expand(
  "~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet"))
# Every interacting protein pair as a text key, in both orders (so the order of a pair never matters).
sp <- as.data.table(read_parquet(STRING_PARQUET))[, .(a = as.character(protein1), b = as.character(protein2))]
# (the keys, e.g. "P12345 Q67890", for both orders of every pair)
skey <- unique(c(paste(sp$a, sp$b), paste(sp$b, sp$a)))
# Helper: read one arm's normalised metabolite vectors (step 1b) as a matrix.
vec <- function(arm) { x <- fread(file.path(OUT, sprintf("01b_metab_nodes_%s.csv", arm)))
  m <- as.matrix(x[, -1]); rownames(m) <- x$metabolite; m }
# Both arms.
ME <- vec("EE"); MR <- vec("RE")

# ---- the rule, parameterised ------------------------------------------------------------------------
# Helper: apply the step 6 rule for a given protein set and class level. link = "shared": the two
# metabolites share a protein; link = "string": they share a protein OR their proteins interact in STRING.
# Returns the counts row and the edge list (with the protein pairs that justify each edge).
run_rule <- function(protein_set, class_col, label, link = "shared") {
  # metabolite -> protein links restricted to the chosen proteins
  mp <- m_enz[uniprot %in% protein_set]
  # metabolite pairs that share a protein (a join on the protein; each metabolite pair once)
  pp <- merge(mp[, .(m1 = metabolite, u1 = uniprot)], mp[, .(m2 = metabolite, u2 = uniprot)],
              by.x = "u1", by.y = "u2", allow.cartesian = TRUE)[m1 < m2][, u2 := u1]
  # for link = "string", also pairs whose proteins differ but interact in STRING: compare every
  # (metabolite, protein) with every other (only done for small protein sets, e.g. our 471 genes)
  if (link == "string") {
    # all pairs of (metabolite, protein) x (metabolite, protein), each metabolite pair once
    allp <- merge(mp[, .(k = 1L, m1 = metabolite, u1 = uniprot)], mp[, .(k = 1L, m2 = metabolite, u2 = uniprot)],
                  by = "k", allow.cartesian = TRUE)[m1 < m2 & u1 != u2]
    # keep those whose two proteins interact in STRING, and add them to the shared-protein pairs
    pp <- rbind(pp[, .(m1, m2, u1, u2)], allp[paste(u1, u2) %in% skey, .(m1, m2, u1, u2)])
  }
  # one row per metabolite pair, with the protein pairs that link it
  pr <- pp[, .(via = paste(sort(unique(fifelse(u1 == u2, u1, paste0(u1, "~", u2)))), collapse = ";")), by = .(m1, m2)]
  # each metabolite's class at the chosen level (blank = Unclassified, never matches)
  cl <- setNames(ids[[class_col]], ids$metabolite)
  # keep pairs in the same class
  e <- pr[!is.na(cl[m1]) & cl[m1] != "" & cl[m1] == cl[m2]]
  # both arms' edge weights (dot products of the normalised vectors)
  e[, `:=`(w_EE = rowSums(ME[m1, , drop = FALSE] * ME[m2, , drop = FALSE]),
           w_RE = rowSums(MR[m1, , drop = FALSE] * MR[m2, , drop = FALSE]))]
  # one row of counts, plus how the arms compare on this network
  row <- data.table(rule = label, proteins_allowed = length(unique(protein_set)),
                    metabolites_with_a_protein = uniqueN(mp$metabolite), pairs_linked = nrow(pr),
                    edges = nrow(e), metabolites_in_network = uniqueN(c(e$m1, e$m2)),
                    cor_w_EE_w_RE = if (nrow(e) > 2) round(cor(e$w_EE, e$w_RE), 3) else NA_real_,
                    edges_sign_change = sum(sign(e$w_EE) != sign(e$w_RE)))
  # return both
  list(row = row, edges = e)
}

# ---- run the rules ---------------------------------------------------------------------------------
runs <- list(
  # the original main-class rule
  run_rule(ours, "main_class", "baseline: our 471 genes, main class"),
  # experiment 1: every human reviewed Rhea enzyme may link metabolites
  run_rule(intersect(unique(enz$uniprot), human), "main_class", "exp 1: all human Rhea enzymes, main class"),
  # experiment 2: broader chemical class (the adopted step 6 rule)
  run_rule(ours, "super_class", "exp 2: our 471 genes, super class"),
  # side line: enzymes from any organism (for reference only)
  run_rule(unique(enz$uniprot), "main_class", "side line: Rhea enzymes from any organism, main class"),
  # experiment 3: the step 6 rule, but STRING-interacting proteins also link metabolites
  run_rule(ours, "super_class", "exp 3: exp 2 + STRING-interacting proteins (>= 700) (step 6 rule)", link = "string"))
# The counts table, one row per rule.
res <- rbindlist(lapply(runs, `[[`, "row"))
# The edges experiment 3 adds to experiment 2 (same proteins and class; only the STRING step differs).
e2 <- runs[[3]]$edges; e3 <- runs[[5]]$edges
extra <- e3[!paste(m1, m2) %in% paste(e2$m1, e2$m2)][, class := ids$super_class[match(m1, ids$metabolite)]]
# (UniProt accessions shown as gene symbols for readability, using the package's lookup table)
u2s <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[!is.na(uniprot), .(u = sub("-[0-9]+$", "", as.character(uniprot)), g = as.character(gene_symbol))])
sym <- function(x) { for (i in seq_len(nrow(u2s))) x <- gsub(paste0("\\b", u2s$u[i], "\\b"), u2s$g[i], x); x }
extra[, via_genes := sym(via)]
# Save the extra edges.
fwrite(extra[, .(metabolite_a = m1, metabolite_b = m2, class, linked_by = via_genes, w_EE, w_RE, w_diff = w_EE - w_RE)],
       file.path(OUT, "08_string_neighbour_extra_edges.csv"))
# Change in the number of metabolites in the network relative to the baseline.
res[, change_vs_baseline := metabolites_in_network - metabolites_in_network[1]]
# Safety check: the rule matching step 6's class level reproduces step 6 exactly.
s6 <- fread(file.path(OUT, "06_metabolite_summary.csv")); v6 <- setNames(s6$value, s6$metric)
# (row 1 = main class, row 3 = super class, row 5 = super class + STRING-interacting proteins; our 471 genes)
row6 <- if (v6[["class_level"]] == "main_class") 1 else if (v6[["link_rule"]] == "shared") 3 else 5
# (edges and metabolites in the network must match step 6's summary)
stopifnot(res$edges[row6] == as.numeric(v6[["edges_same_class"]]),
          res$metabolites_in_network[row6] == as.numeric(v6[["metabolites_with_edges"]]))
# Save and show.
fwrite(res, file.path(OUT, "08_metabolite_rule_experiments.csv"))
# Show it.
print(res)
