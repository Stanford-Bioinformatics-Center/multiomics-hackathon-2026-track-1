# Regression: the engine must reproduce the pipeline's own edge weights (steps 3, 6, 14) exactly.
# Skipped when the pipeline outputs are not present (e.g. on a machine that has not run it).
out <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
have <- all(file.exists(file.path(out, c("14_joint_edges.csv", "01_nodes_EE.csv", "01_nodes_RE.csv", "01b_metab_nodes_EE.csv", "01b_metab_nodes_RE.csv"))))
test_that("engine weights equal the joint network of step 14", {
  skip_if_not(have, "pipeline outputs not found (run network/run_all.sh first)")
  rd <- function(f, id) { x <- data.table::fread(file.path(out, f)); m <- as.matrix(x[, setdiff(names(x), c("entrez_gene", "gene_symbol", "metabolite")), with = FALSE]); rownames(m) <- x[[id]]; m }
  G <- list(EE = rd("01_nodes_EE.csv", "gene_symbol"), RE = rd("01_nodes_RE.csv", "gene_symbol"))
  M <- list(EE = rd("01b_metab_nodes_EE.csv", "metabolite"), RE = rd("01b_metab_nodes_RE.csv", "metabolite"))
  J <- data.table::fread(file.path(out, "14_joint_edges.csv"))
  run <- function(type, EE, RE) { j <- J[edge_type == type]; build_network(physical_edges(j$node_a, j$node_b, type), EE, RE, type) }
  gdim <- colnames(G$EE)
  pp <- run("protein - protein", embedding(G$EE, "EE"), embedding(G$RE, "RE"))
  mm <- run("metabolite - metabolite", embedding(M$EE, "EE", "metabolite"), embedding(M$RE, "RE", "metabolite"))
  both <- function(arm) { g <- embedding(G[[arm]], arm); bind_embeddings(double_embedding(embedding(M[[arm]], arm, "metabolite"), g), g) }
  mp <- run("metabolite - protein", both("EE"), both("RE"))
  e <- edge_table(combine_networks(pp, mm, mp))
  chk <- merge(J, e, by.x = c("node_a", "node_b", "edge_type"), by.y = c("a", "b", "edge_type"), suffixes = c("_pipeline", "_engine"))
  expect_equal(nrow(chk), nrow(J))
  expect_lt(max(abs(chk$w_EE_pipeline - chk$w_EE_engine)), 1e-12)
  expect_lt(max(abs(chk$w_RE_pipeline - chk$w_RE_engine)), 1e-12)
})
