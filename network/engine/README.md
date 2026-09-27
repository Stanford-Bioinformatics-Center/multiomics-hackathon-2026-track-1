# exnet — the exercise network engine

The core rule of the project, implemented once, typed and tested:

```
edge(u, v) exists  <=>  a physical database links u and v            HARD  (STRING >= 700, Rhea)
w_arm(u, v)        =    sum_d z_u,d(arm) * z_v,d(arm)                 SOFT  (dot product of exercise responses)
w_diff             =    w_EE - w_RE
```

## Classes (S4, with validity rules)

| Class | Holds | Guaranteed by its validity rule |
|---|---|---|
| `Embedding` | one arm's response vectors: molecules x dimensions (tissue x ome x time) | numeric; unique molecule and dimension names; arm is `EE` or `RE`; no infinite values; missing values only as whole empty dimensions |
| `PhysicalEdges` | the hard layer: `a`, `b`, `source`, `evidence` | no self-loops; no duplicate undirected edges; every edge names its database |
| `WeightedNetwork` | every physical edge with `w_EE`, `w_RE`, `w_diff`, `n_dims` | finite weights; `w_diff == w_EE - w_RE` |

## Functions

| Function | Does |
|---|---|
| `embedding(values, arm, kind)` | validated response vectors |
| `physical_edges(a, b, source, evidence)` | validated hard layer |
| `double_embedding(metab, genes)` | metabolite vectors into the gene dimension space (RNA and protein slots), empty where no gene is measured |
| `bind_embeddings(...)` | stack embeddings of one arm (e.g. doubled metabolites + genes) |
| `build_network(hard, EE, RE, edge_type)` | the rule above, for one edge type |
| `combine_networks(...)` | protein-protein + metabolite-metabolite + metabolite-protein in one network |
| `edge_table(net)`, `node_strength(net)` | edges as a table; strength = sum of \|w\| per node (the hub definition) |
| `arm_specific_edges(net, q = 0.75, tau = NULL)` | strong (>= tau, default the q-quantile of \|w\| over both arms) after one arm only |
| `toy_example()` | 5 proteins + 1 metabolite with hand-computed weights (tests, method figure) |

Every function checks its inputs and fails with a classed error that says what is wrong:
`exnet_input_error`, `exnet_arm_error`, `exnet_dimension_error`, `exnet_missing_node_error` (all inherit `exnet_error`).

## Example

```r
pkgload::load_all("network/engine")          # or: install.packages("network/engine", repos = NULL, type = "source")
t   <- toy_example()
net <- build_network(t$hard_pp, t$gene_EE, t$gene_RE, "protein - protein")
edge_table(net)                               # A-B: w_EE 0.50, w_RE 0.38 (see the hand calculation in toy.R)
arm_specific_edges(net, tau = 0.3)
```

## Tests

`Rscript network/engine/run_tests.R` (also `run_all.sh` step 14t) — 25 tests: the toy example against hand-computed
weights, the hard gate, arm-specific edges, node strength, every bad-input path, and a regression test that rebuilds
the pipeline's joint network (step 14, 704 edges) and matches every weight to within 1e-12.
