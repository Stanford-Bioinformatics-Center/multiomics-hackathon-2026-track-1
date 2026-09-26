# Exercise network: node embeddings

Step 1 of the network build: one node per gene, carrying an 18-value exercise-response
vector. The next step links the nodes with STRING.

## Run

```bash
Rscript network/01_node_embeddings.R
```

Requires R with `data.table` and `MotrpacHumanPreSuspensionAnalysis` (the MoTrPAC human
pre-suspension package, which ships the differential analysis results and the feature-to-gene map).
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

## Next

Map the 471 genes to STRING protein IDs and pull the edges between them (step 2).
