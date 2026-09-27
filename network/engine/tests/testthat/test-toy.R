# The toy example: every weight must equal the hand calculation, and the hard gate must hold.
test_that("toy weights equal the hand-computed dot products", {
  t <- toy_example()
  pp <- build_network(t$hard_pp, t$gene_EE, t$gene_RE, "protein - protein")
  # Metabolite - protein: double the metabolite into the gene dimensions, then one dot product.
  dims <- colnames(t$gene_EE@values)
  bindE <- function(m, g) bind_embeddings(double_embedding(m, g), g)
  mp <- build_network(t$hard_mp, bindE(t$metab_EE, t$gene_EE), bindE(t$metab_RE, t$gene_RE), "metabolite - protein")
  net <- combine_networks(pp, mp); e <- edge_table(net)
  ex <- t$expected
  got <- merge(ex, e, by = c("a", "b"), suffixes = c("_hand", ""))
  expect_equal(nrow(got), 5L)
  expect_equal(got$w_EE, got$w_EE_hand, tolerance = 1e-12)
  expect_equal(got$w_RE, got$w_RE_hand, tolerance = 1e-12)
  expect_equal(e$w_diff, e$w_EE - e$w_RE)
})

test_that("the hard gate: similar molecules without a physical link get no edge", {
  t <- toy_example(); net <- build_network(t$hard_pp, t$gene_EE, t$gene_RE)
  e <- edge_table(net)
  # A and E respond almost identically (dot product 0.57 in EE) but are not physically linked.
  expect_equal(sum(t$gene_EE@values["A", ] * t$gene_EE@values["E", ]), 0.57, tolerance = 1e-12)
  expect_false(any((e$a == "A" & e$b == "E") | (e$a == "E" & e$b == "A")))
  expect_false("E" %in% c(e$a, e$b))
})

test_that("arm-specific edges with an explicit threshold", {
  t <- toy_example(); dims <- colnames(t$gene_EE@values)
  bindE <- function(m, g) bind_embeddings(double_embedding(m, g), g)
  net <- combine_networks(build_network(t$hard_pp, t$gene_EE, t$gene_RE), build_network(t$hard_mp, bindE(t$metab_EE, t$gene_EE), bindE(t$metab_RE, t$gene_RE)))
  s <- arm_specific_edges(net, tau = 0.3); lab <- setNames(s$specificity, paste(s$a, s$b))
  expect_equal(unname(lab[c("A B", "C D", "M D", "B C", "A C")]), c("both", "resistance-specific", "endurance-specific", "neither", "neither"))
})

test_that("node strength is the sum of |w| over a node's edges", {
  t <- toy_example(); s <- node_strength(build_network(t$hard_pp, t$gene_EE, t$gene_RE))
  expect_equal(s[node == "C", strength_EE], 0.13 + 0.13 + 0.19, tolerance = 1e-12)
})
