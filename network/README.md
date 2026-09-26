# Exercise network: endurance vs resistance

Goal: two networks over the same genes, one for endurance vs control (EE) and one for resistance
vs control (RE), to compare against each other.

- **Step 1 — nodes** (`01_node_embeddings.R`): one node per gene. Each node has **two** 18-value
  exercise-response vectors, one per arm, built independently from that arm's contrast.
- **Step 2 — edges** (`02_string_edges.R`): STRING links between the nodes, built with the
  network rules of El-Kebir et al. 2015. STRING gates whether an edge exists.
- **Step 3 — edge weights** (`03_edge_weights.R`): each STRING edge is weighted by the dot product
  of its two genes' embeddings, once per arm, giving the EE and RE weighted networks.

## Run

```bash
Rscript network/01_node_embeddings.R
Rscript network/02_string_edges.R
Rscript network/03_edge_weights.R
```

Requires R with `data.table`, `igraph`, `nanoparquet` and `MotrpacHumanPreSuspensionAnalysis`
(the MoTrPAC human pre-suspension package, which ships the differential analysis results and the
feature-to-gene map). Step 2 reads the STRING file from `$STRING_PARQUET` (default
`~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet`).
Outputs go to `$HACK_OUT` (default `~/Desktop/output/hackathon-2026-track1/network`); code and
results are kept in separate trees, so nothing is written inside the repo.

## Node universe: 471 genes

A gene is included only if it is measured in **all six** tissue × ome combinations, in **both**
exercise arms, at every post-exercise timepoint that combination has:

| Tissue  | RNA                          | Protein                                  |
|---------|------------------------------|------------------------------------------|
| Adipose | RNA-seq, 0.5 / 4 / 24 h      | MS proteomics, **4 h only** (see below)  |
| Blood   | RNA-seq, 0.5 / 4 / 24 h      | OLINK (targeted panel), 0.5 / 4 / 24 h   |
| Muscle  | RNA-seq, 0.5 / 4 / 24 h      | MS proteomics, 0.5 / 4 / 24 h            |

Genes are matched across omes by Entrez ID (`HUMAN_FEATURE_TO_GENE`). The count is bounded by
OLINK, which is a fixed panel of ~1,400 proteins: with MS proteomics only (blood RNA only) the
same rule gives 4,878 genes. Phosphoproteomics is not included.

## Embedding: 18 dimensions

Order is tissue (adipose → blood → muscle), then ome (RNA → protein), then timepoint
(0.5 → 4 → 24 h):

```
adipose_rna_0.5h  adipose_rna_4h  adipose_rna_24h  adipose_prot_0.5h  adipose_prot_4h  adipose_prot_24h
blood_rna_0.5h    blood_rna_4h    blood_rna_24h    blood_prot_0.5h    blood_prot_4h    blood_prot_24h
muscle_rna_0.5h   muscle_rna_4h   muscle_rna_24h   muscle_prot_0.5h   muscle_prot_4h   muscle_prot_24h
```

**Value.** logFC of the delta-delta contrast: the exercise arm's change from pre-exercise minus
the control group's change from pre-exercise at the same timepoint. There is one file per arm,
because each arm has its own vector:

- `EE`: endurance vs control (`contrast_category == "EE-CON"`)
- `RE`: resistance vs control (`contrast_category == "RE-CON"`)

**Scaling across omes.** Raw logFC magnitudes differ by ome (the RMS of muscle protein is 0.10,
of blood OLINK 0.38), so unscaled vectors would be dominated by blood protein. Each tissue × ome
block is divided by its root-mean-square logFC, computed over the 471 genes, all timepoints that
block has, and both arms. Because there is **one factor per block**:

- all six omes end up on a comparable scale (block RMS = 1);
- differences between timepoints and between arms *within* a block are preserved;
- zero still means "no change" and signs are unchanged (no centering).

**Adipose protein at 0.5 h and 24 h is empty (NA) for every gene.** Those samples do not exist:
adipose proteomics is 22 participants sampled only pre-exercise and at 4 h. The columns are kept
so that every file has the same 18-dimension layout; downstream similarity measures must ignore
NAs (or drop these two columns).

**Several features per gene.** When more than one transcript or protein maps to a gene within a
tissue × ome, the feature with the highest `AveExpr` (mean abundance) is kept. Abundance does not
depend on the exercise response, so this does not select for responsive features. This affects
12 genes in muscle protein, 9 in adipose protein and 2 in each tissue's RNA.

## Outputs

| File | Contents |
|------|----------|
| `01_nodes_EE.csv`, `01_nodes_RE.csv` | **The node tables.** 471 rows: `entrez_gene`, `gene_symbol`, 18 scaled dimensions |
| `01_nodes_EE_raw_logFC.csv`, `01_nodes_RE_raw_logFC.csv` | Same layout, unscaled logFC |
| `01_scale_factors.csv` | The RMS divisor for each tissue × ome block |
| `01_nodes_feature_provenance.csv` | The feature id behind each gene in each block, and how many candidates it was chosen from |

The script checks that the universe is exactly 471 genes and that each table is 471 × 18.

## Step 2: edges

Rules follow El-Kebir et al. 2015 (xHeinz, *Bioinformatics* 31:3147), section 3.3:

1. **Background network:** STRING protein–protein interactions. We use the curated file with
   `combined_score >= 700` (124,099 edges, 13,855 proteins, UniProt IDs).
2. **Undirected, no self-loops, one edge per pair.** The file already satisfies this.
3. **Outlier hubs removed:** nodes with degree above the 75th percentile + 40 × IQR of the degree
   distribution, computed on the full background network before restricting to our genes. Here the
   cutoff is 21 + 40 × 18 = 741; the highest degree is 407, so **no hubs are removed** (the paper
   removed ELAVL1 and ubiquitin from the raw STRING, and this file has already been curated).
4. **Induced subnetwork:** only edges where both ends are among the 471 genes.

The paper's network is **unweighted**; `combined_score` is kept as an attribute. Proteins are
mapped to genes through the UniProt accessions in `HUMAN_FEATURE_TO_GENE` (isoform suffix
stripped); where several protein pairs map to one gene pair, the highest score is kept.

**Same topology for both arms.** STRING does not depend on the exercise data, so EE and RE share
one edge set. The arms differ in their node vectors, and in two per-edge attributes added here
(not part of El-Kebir): `cos_EE` and `cos_RE`, the cosine similarity of the two endpoints'
embeddings in that arm, over the dimensions present in both.

**Result:**

| | |
|---|---|
| Nodes | 471 (444 in STRING; 27 have no edge ≥ 700 anywhere in the file) |
| Edges | 431 |
| Isolated nodes | 185 |
| Largest connected component | 230 nodes |
| Other components | 20 small ones (sizes 2–7) |
| Median degree | 1 (top hubs: ITGB1 21, CD34 16, NT5E 16, ITGAM 14) |

The network is small and sparse because OLINK limits the node set to secreted and cell-surface
proteins (integrins, adhesion molecules, cytokines).

| File | Contents |
|------|----------|
| `02_edges.csv` | One row per edge: Entrez, symbol and UniProt for both ends, `combined_score`, `cos_EE`, `cos_RE` |
| `02_nodes_string.csv` | Per node: UniProt, whether it is in STRING, degree, component id and size |
| `02_network_summary.csv` | The counts above |

## Step 3: edge weights

STRING decides whether an edge exists; the exercise data decides how strong it is in each arm.
Following the encoder–decoder framing of node embeddings (similarity(u, v) ≈ z_uᵀz_v, Stanford
CS224W), each edge gets the **dot product of its endpoints' embeddings**, computed per arm:

```
w_EE(u,v) = z_u(EE) · z_v(EE)        w_RE(u,v) = z_u(RE) · z_v(RE)        w_diff = w_EE − w_RE
```

- **Whole vector, one number.** The dot product runs over the entire embedding (all tissues ×
  RNA and protein × every timepoint), not per tissue or ome. The two adipose protein columns that
  exist for no gene (0.5 h, 24 h) are left out, so every edge uses the same 16 dimensions.
- **Signed.** Positive = the two genes respond in the same direction; negative = opposite
  directions; near 0 = at least one barely responds. Large values need both genes to respond strongly.
- **Untransformed**, so both arms are on the same scale (the embeddings share scale factors) and
  `w_diff` is directly interpretable.

| | EE | RE |
|---|---|---|
| Median weight | 1.24 | 1.89 |
| Range | −14.8 to 66.2 | −54.0 to 66.8 |
| Negative edges | 145 of 431 | 137 of 431 |

The two arms' weights correlate at r = 0.64, and 130 of 431 edges change sign between arms (e.g.
HSPA1A–DNAJB1: +22.9 in EE, −3.1 in RE).

| File | Contents |
|------|----------|
| `03_weighted_edges.csv` | One row per edge: Entrez and symbol for both ends, `combined_score`, `w_EE`, `w_RE`, `w_diff` |
