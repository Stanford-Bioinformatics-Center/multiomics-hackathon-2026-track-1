# Bad inputs must fail early with a classed error that says what is wrong.
test_that("invalid embeddings are rejected", {
  m <- matrix(1:4 / 10, 2, dimnames = list(c("A", "B"), c("d1", "d2")))
  expect_error(embedding(m, "XX"), class = "exnet_input_error")                           # unknown arm
  expect_error(embedding(matrix(letters[1:4], 2), "EE"), class = "exnet_input_error")     # not numeric
  bad <- m; bad[1, 1] <- NA
  expect_error(embedding(bad, "EE"), "partly missing")                                     # a gap inside a column
  inf <- m; inf[2, 2] <- Inf
  expect_error(embedding(inf, "EE"), "infinite")
  dup <- m; rownames(dup) <- c("A", "A")
  expect_error(embedding(dup, "EE"), "unique")
})
test_that("invalid physical edges are rejected", {
  expect_error(physical_edges("A", "A", "STRING"), "self-loops")
  expect_error(physical_edges(c("A", "B"), c("B", "A"), "STRING"), "duplicate")
  expect_error(physical_edges("A", c("B", "C"), "STRING"), class = "exnet_input_error")
})
test_that("mismatched arms, dimensions and missing nodes are caught", {
  t <- toy_example()
  expect_error(build_network(t$hard_pp, t$gene_RE, t$gene_RE), class = "exnet_arm_error")
  other <- embedding(t$gene_RE@values[, 4:1], "RE")
  expect_error(build_network(t$hard_pp, t$gene_EE, other), class = "exnet_dimension_error")
  hard <- physical_edges("A", "Z", "STRING")
  gE <- embedding(t$gene_EE@values, "EE"); gR <- embedding(t$gene_RE@values, "RE")
  expect_error(build_network(hard, gE, gR), class = "exnet_missing_node_error")
  liver <- embedding(matrix(0.1, 1, 1, dimnames = list("A", "liver_rna_4h")), "EE")
  expect_error(double_embedding(t$metab_EE, liver), class = "exnet_dimension_error")
  expect_error(double_embedding(t$metab_EE, t$gene_RE), class = "exnet_arm_error")
})
