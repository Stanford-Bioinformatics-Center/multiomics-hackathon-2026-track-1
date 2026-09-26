# Exercise network: endurance vs resistance

Goal: an honest representation of two networks over the same genes, one for endurance vs control
(EE) and one for resistance vs control (RE), and a test of whether they are more similar or more
different. Nothing in the pipeline is designed to make the two networks look alike: they are built
independently from each arm's data, measured on one shared scale, and compared against measurement
noise.

Every script is commented line by line in plain language, with a header explaining what it does,
the upstream QC, the methods and the tools.

- **Step 1 — nodes** (`01_node_embeddings.R`): one node per gene. Each node has **two** 18-value
  exercise-response vectors, one per arm, built independently from that arm's contrast.
- **Step 1b — metabolite nodes** (`01b_metabolite_embeddings.R`): one node per metabolite (450).
  Metabolomics is one layer, so each vector has 9 values (3 tissues × 3 times), half a gene's 18;
  again two vectors per metabolite, one per arm, built independently.
- **Step 2 — edges** (`02_string_edges.R`): STRING links between the nodes, built with the
  network rules of El-Kebir et al. 2015. STRING gates whether an edge exists.
- **Step 3 — edge weights** (`03_edge_weights.R`): each STRING edge is weighted by the dot product
  of its two genes' embeddings, once per arm, giving the EE and RE weighted networks.
- **Step 4 — compare arms** (`04_compare_arms.R`): which edges and genes are wired differently in
  EE vs RE, with bootstrap uncertainty from each estimate's standard error.

## Run

```bash
Rscript network/01_node_embeddings.R
Rscript network/01b_metabolite_embeddings.R
Rscript network/02_string_edges.R
Rscript network/03_edge_weights.R
Rscript network/04_compare_arms.R      # ~10 s; N_BOOT (default 10000) sets the draws
```

Requires R with `data.table`, `igraph`, `nanoparquet` and `MotrpacHumanPreSuspensionAnalysis`
(the MoTrPAC human pre-suspension package, which ships the differential analysis results and the
feature-to-gene map). Step 2 reads the STRING file from `$STRING_PARQUET` (default
`~/Downloads/Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet`).
Outputs go to `$HACK_OUT` (default `~/Desktop/output/hackathon-2026-track1/network`); code and
results are kept in separate trees, so nothing is written inside the repo.

## Tech stack

| Component | Version | Used for |
|---|---|---|
| R | 4.4.3 | everything |
| `MotrpacHumanPreSuspensionAnalysis` | 0.2.4 | the data: differential-analysis results and the feature-to-gene map |
| `data.table` | 1.18 | tables |
| `nanoparquet` | 0.4 | reading the STRING `.parquet` file |
| `igraph` | 2.2 | connected components |
| `Matrix` | 1.7 | sparse gene × edge matrix (step 4) |
| STRING (curated file) | combined_score ≥ 700 | which genes are connected |

## Upstream QC and filtering (done by the consortium, before our code)

We read published results; no statistics are re-run on raw data.

- **Samples:** first supervised exercise bout only (visit `ADU_BAS`); samples flagged in the
  consortium's outlier review (`OUTLIERS`: failed genotype/sex checks, etc.) removed.
- **RNA-seq (adipose, blood, muscle):** low-expression genes removed by a log-CPM cutoff; TMM
  normalisation (edgeR); voom precision weights; linear mixed model with dream (variancePartition):
  `~ 0 + group_timepoint + Batch + BMI + calculatedAge + codedsiteid + pct_umi_dup + RIN + Sex + (1 | pid)`.
- **MS proteomics (adipose, muscle):** QC-normalised; a protein is kept only if every group ×
  timepoint has ≥ 3 participants measured both before and after exercise (`filter_paired_n`);
  model `~ 0 + group_timepoint + BMI + calculatedAge + Sex + (1 | pid)` (dream).
- **OLINK proteomics (blood):** targeted ~1,400-protein panel, QC-normalised; same model as MS.
- **Metabolomics (adipose, blood/plasma, muscle; 10–12 platforms per tissue):** features missing in
  > 20% of samples removed; KNN imputation; log2; median/MAD normalisation (untargeted); names
  standardised to RefMet; where a metabolite was measured on several platforms, the lowest-CV
  (most reproducible) platform is kept (`METABOLOMICS_CVS`). Model: targeted
  `~ 0 + group_timepoint + BMI + calculatedAge + codedsiteid + Sex + (1 | pid)`; untargeted adds
  `raw_intensity_post` (each sample's overall signal level).
- **STRING file:** already restricted to combined_score ≥ 700, UniProt IDs, one row per pair.

## Filtering and choices made in this pipeline

| Step | Choice | Why |
|---|---|---|
| 1 | Delta-delta contrasts (EE-CON, RE-CON) at 0.5 / 4 / 24 h | removes changes that happen to controls too |
| 1 | Genes measured in all six tissue × layer combinations (471) | every node has a complete vector |
| 1 | One feature per gene: the most abundant (`AveExpr`) | independent of the exercise response |
| 1 | One scale factor per tissue × layer, shared by both arms | puts RNA and protein on one footing without erasing arm differences |
| 1 | Adipose protein 0.5 / 24 h left empty | those samples do not exist |
| 1b | Metabolites measured in all three tissues, matched by exact RefMet name (450) | complete vectors; exact names because lipid punctuation is meaningful |
| 1b | One scale factor per tissue, shared by both arms | same common-unit logic as genes |
| 2 | El-Kebir 2015 network rules (hub removal removes nothing here) | published, standard construction |
| 3 | Dot product over the whole 16-dimension vector | encoder–decoder node-similarity framework |
| 3 | 0-1 weights: σ(w / s), s = median \|w\| over both arms | temperature scaling with the median heuristic (see step 3) |
| 4 | Parametric bootstrap from standard errors, EE/RE correlation included | error bars for every difference |

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
| `01_nodes_EE_se.csv`, `01_nodes_RE_se.csv` | Standard error of every embedding value, same scale (from the model's 95% CI) |
| `01_nodes_arm_corr.csv` | Correlation between each gene's EE and RE estimate per dimension (shared controls; median ~0.5–0.65) |
| `01_nodes_feature_provenance.csv` | The feature id behind each gene in each block, and how many candidates it was chosen from |

The script checks that the universe is exactly 471 genes and that each table is 471 × 18.

## Step 1b: metabolite nodes (450 metabolites)

The metabolite counterpart of step 1. A metabolite is included if it is measured in all three tissues
(both arms, all three post-exercise times), matched across tissues by **exact** RefMet name: no
lower-casing or punctuation stripping, because in lipid names punctuation distinguishes different
molecules. 450 metabolites qualify (71% lipids by RefMet class; 442 measured on untargeted platforms,
8 on targeted panels in muscle). No duplicate handling is needed: upstream QC already keeps one
(lowest-CV) platform per metabolite per tissue.

**Embedding: 9 dimensions per arm**

```
adipose_metab_0.5h  adipose_metab_4h  adipose_metab_24h
blood_metab_0.5h    blood_metab_4h    blood_metab_24h
muscle_metab_0.5h   muscle_metab_4h   muscle_metab_24h
```

Values are the same delta-delta logFC as for genes (EE-CON, RE-CON), each tissue divided by one
root-mean-square factor pooled over metabolites, times and both arms (adipose 0.249, blood 0.231,
muscle 0.307). Nothing is missing by design. Standard errors and the EE/RE estimate correlation
(median 0.65) are saved for later comparison steps, exactly as for genes. Median |value| / SE is
0.76, similar to the genes.

| File | Contents |
|------|----------|
| `01b_metab_nodes_EE.csv`, `01b_metab_nodes_RE.csv` | **The metabolite node tables.** 450 rows: `metabolite`, 9 scaled dimensions |
| `01b_metab_nodes_{EE,RE}_raw_logFC.csv` | Same layout, unscaled logFC |
| `01b_metab_nodes_{EE,RE}_se.csv` | Standard error of every value, same scale |
| `01b_metab_nodes_arm_corr.csv` | Correlation between each metabolite's EE and RE estimate per dimension |
| `01b_metab_scale_factors.csv` | The divisor for each tissue |
| `01b_metab_nodes_provenance.csv` | The platform that measured each metabolite in each tissue |

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
- **Sigmoid version for positive-only methods:** `sig_EE`, `sig_RE` = σ(w / s) = 1 / (1 + e^(−w/s)),
  in 0–1 (negative w → below 0.5; 0.5 = no co-response), s = 2.667.

**Why divide by the median, and what is the precedent.**

- *Why a sigmoid at all:* some network methods (random walks, community detection) need positive
  weights. σ(dot product) is the decoder used by DeepWalk / node2vec-style embeddings.
- *Why divide first:* those methods *learn* vectors whose dot products sit in the few-units range
  where the sigmoid is informative. Ours are measured, and their dot products run from −54 to +67, so
  plain σ(w) pinned 116 EE / 146 RE of 431 edges at ~0 or ~1 (with s: 19 / 31). Dividing a sigmoid's
  input by a constant is **temperature scaling** (Hinton, Vinyals & Dean 2015; Guo et al. 2017). It
  changes the steepness only; the order and sign of every edge are unchanged.
- *Why the median:* setting that constant from the median of the data is the **median heuristic**,
  the standard default for the width of similarity kernels (Schölkopf & Smola 2002; Gretton et al.
  2012, *JMLR*). It is robust to the few very large weights. It is a heuristic, not a derived
  optimum: it puts the typical edge at σ(±1) = 0.27 / 0.73.
- *Why one s for both arms (and why this does not make the networks similar):* s is a unit of
  measurement. With one s, the same weight maps to the same 0–1 value in either arm, so a real
  difference between arms passes through unchanged. A separate s per arm would rescale each network
  to its own typical edge and could *hide* a genuine overall difference. The same reasoning applies to
  step 1's scale factors, which are also shared by the two arms.
- *Caveat:* the sigmoid is non-linear and compresses large differences, so similarity-vs-difference
  conclusions are drawn from the raw w (step 4), not from the 0–1 version.

| | EE | RE |
|---|---|---|
| Median weight | 1.24 | 1.89 |
| Range | −14.8 to 66.2 | −54.0 to 66.8 |
| Negative edges | 145 of 431 | 137 of 431 |

The two arms' weights correlate at r = 0.64, and 130 of 431 edges change sign between arms (e.g.
HSPA1A–DNAJB1: +22.9 in EE, −3.1 in RE).

| File | Contents |
|------|----------|
| `03_weighted_edges.csv` | One row per edge: Entrez and symbol for both ends, `combined_score`, `w_EE`, `w_RE`, `w_diff`, `sig_EE`, `sig_RE` |

## Step 4: compare EE vs RE

Same topology, two sets of weights, so the comparison is per edge (`w_diff = w_EE − w_RE`) and per
gene (**strength** = sum of a gene's edge weights; `delta_strength` = EE − RE, signed, plus the same
on the rescaled sigmoid weights σ(w / s)).

**Uncertainty: parametric bootstrap.** Every embedding value is an estimate with a standard error.
In each of 10,000 draws, each gene's EE and RE value in every dimension is resampled from a
bivariate normal centred on the estimate, using the two SEs and their correlation (the arms share
one control group; `EE-RE` in the package equals `EE-CON − RE-CON` exactly, so
ρ = (se_EE² + se_RE² − se_EE−RE²) / (2·se_EE·se_RE)). All weights are recomputed per draw.
Two-sided p = 2 × min(share of draws ≤ 0, share ≥ 0); BH within edges and within genes; 95%
percentile intervals reported. Errors are treated as independent across genes.

An arm-label swap permutation was tried first and **rejected**: swapping a gene's EE and RE vectors
only flips the sign of an edge's `w_diff` when both ends swap, so half the draws reproduce |w_diff|
exactly and no edge can reach p < ~0.5.

**Result.** The two arms' edge weights correlate at r = 0.64 (bootstrap 95% interval 0.34–0.69;
noise pulls it down), and 130 of 431 edges change sign between arms. Whether 0.64 counts as "similar"
or "different" needs a reference point that has not been computed yet: the correlation expected if
the arms were truly identical and differed only by measurement noise. Most responses are small
relative to their error (median |value| / SE = 0.78), so the per-edge test has little power, and
"not significant" here is **not** evidence that the arms are the same:

| | |
|---|---|
| Edges at FDR < 0.1 | 0 (12 at nominal p < 0.05) |
| Genes at FDR < 0.1 | 1: **HSPB1** (FDR 0.057; stronger in RE) |

The nominal edges form small, coherent modules (exploratory tier in `04_diff_subnetworks.csv`):

- **Heat-shock chaperones:** HSPA1A–DNAJB1 stronger in EE; HSPA1A–HSPB1 and HSPB1–BAG3 stronger in RE.
- **Lactate transport:** SLC16A1 (MCT1)–BSG (basigin) stronger in EE.
- **CEBPB–FOXO1** stronger in EE; **LPL–SORT1**, **ITGAV–ITGB1/MFGE8**, **CD55–CD59** stronger in RE.

| File | Contents |
|------|----------|
| `04_edge_diff.csv` | Per edge: weights, `w_diff` with 95% interval, `p_boot`, `fdr`, `stronger_in`, `sign_change` |
| `04_node_diff.csv` | Per gene with ≥1 edge: degree, signed and sigmoid strength per arm, deltas with interval, p, FDR |
| `04_diff_subnetworks.csv` | Connected modules of FDR < 0.1 edges (none yet) and of p < 0.05 edges (exploratory) |
| `04_summary.csv` | The counts above |
