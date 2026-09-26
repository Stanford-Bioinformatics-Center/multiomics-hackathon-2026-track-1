# MotrpacHumanPreSuspensionAnalysis 2.0.8

## Changes

- `plot_enrich_heatmap()` gains `return_drawing`: with `TRUE` it returns `draw`, `width`,
  `height` and `n_sets` instead of writing a PDF, and `filename` is optional. It also names
  any `set_ids` it drops.

## Dependencies

- MotrpacBicQC is required at `>= 2.0.0`, the `v2.0.0` release tag (2026-09-23), and
  `DESCRIPTION` declares `Remotes: MoTrPAC/MotrpacBicQC@v2.0.0` so `pak` resolves it from
  GitHub. 2.0.0 drops the `inspectdf` import, which is archived on CRAN and made
  MotrpacBicQC unresolvable from a clean library; `assay_codes` is unchanged from 1.9.0.

- `gridtext` is added to `Imports`: `plot_feature_heatmap()` renders its column title with
  `ComplexHeatmap::gt_render()`, which needs it, and ComplexHeatmap only suggests it.

- `cmapR` (Bioconductor) is added to `Suggests`. `PTMSEA_INPUT` is a list of `cmapR::GCT`
  S4 objects, and on a machine without `cmapR` `R CMD check` failed its data inspection
  with a WARNING ("unable to load required package 'cmapR'") because no field declared
  the package. The object still loads and its slots are reachable without `cmapR`;
  printing it or using the `cmapR` accessors needs it installed.

## Internals

- A GitHub Actions workflow, `R-CMD-check.yaml`, now runs `R CMD check` — tests and
  vignettes included — on every pull request and on pushes to `main`. Until now only the
  pkgdown site build ran, and only after a merge, so a failing test could reach `main`
  unseen.

- The roxygen2 hook in `R/zzz.R` calls `utils::getFromNamespace()` with its namespace,
  clearing the "no visible global function definition for 'getFromNamespace'" NOTE that
  `R CMD check` has reported since before 2.0.

- `R-CMD-check.yaml` runs the check under a virtual display (Xvfb): Mfuzz loads Tk, which
  otherwise warns at install time on a headless runner.

# MotrpacHumanPreSuspensionAnalysis 2.0.7

## New data

- `ORA_COLORS` — the white-to-`#543483` ramp for ORA heatmaps, as `c(low, high)`; pass it to
  `TMSig::enrichmap(colors = )`. `plot_cluster_enrichment()` now reads it.

- `PTMSEA_RESULTS` — PTM-SEA results for the prot-ph EE-CON and RE-CON contrasts (muscle
  506 signatures x 6 contrasts, adipose 437 x 2), in the long layout of `CAMERA_RESULTS` with
  `NES` in place of `t`, `df` and `z.std`. Provenance in `data-raw/PTMSEA/README.md`.

## Removed data

- The vendored `assay_codes` object is removed, along with its man page and
  `data-raw/assay_codes.R`. It was a 44-row snapshot of `MotrpacBicQC::assay_codes`, taken
  when `inspectdf` was archived on CRAN and MotrpacBicQC could not be installed from a
  current snapshot. The snapshot has been retired in favour of reading upstream live. Code
  that referenced `MotrpacHumanPreSuspensionAnalysis::assay_codes` should read
  `MotrpacBicQC::assay_codes`.

## Data objects

- `HUMAN_FEATURE_TO_GENE` drops `confident_site` and has 12 columns where 2.0.3 gave it 13.
  The flag is measured per tissue, and this table is keyed on `(assay, feature_id)` with no
  tissue column, so the only value it could carry was the collapse across tissues — not the
  measurement for either tissue on the 859 of 7,865 shared prot-ph sites where muscle and
  adipose disagree. Read the per-tissue value from `*_PROT_PH_QC$feature_metadata` in
  `MotrpacHumanPreSuspensionData`, which is what `preprocess_PTMSEA()` does; `PTMSEA_INPUT`
  never read this table and is unaffected. All 1,920,618 rows and the other 11 columns are
  unchanged from 2.0.3.

## Changes

- `plot_single_feature()` reads `MotrpacBicQC::assay_codes` directly. Upstream 1.9.0 adds
  rows for `metab-t-clinical`, `prot-clinical`, `metab-t-conv` and `metab-t-imm-crt`, none
  of which existed in the vendored snapshot, so the hard-coded label fallback for the
  clinical omes is deleted — those facet strips are now labelled from upstream and read
  `Clin. Chem` rather than `Clin. Chem.`. The `metab-t-conv` label is still overridden to
  `Conv. Metab (log2)`: upstream labels it `Conv(T)` in the `LAB` family, which does not
  distinguish it from its log2 twin `metab-t-clinical`.

- Upstream also revises the 44 pre-existing rows in two columns this package does not
  display: `assay_name` punctuation on six rows (comma to hyphen) and `cas_code` on four
  (`transcript-rna-seq` and `transcript-rna-seq-splicing` from `mssm` to `stanford`,
  `prot-pr` and `prot-ph` from `pnnl` to `broad_prot`). Every `assay_short_text` is
  unchanged, so no existing figure label moves.

- `load_differential_analysis(epigen = TRUE)` reads the c2.0 epigenomics DA
  from the public CloudFront release again, with no bucket access or local cache.
  `load_differential_analysis()` and `plot_single_feature()` drop the `gsutil` and `bucket`
  arguments; `repo_local_dir` is kept but ignored, with a message.

- `plot_enrich_heatmap()` accepts `PTMSEA_RESULTS` again, plotting NES; `set_ids` takes
  PTMsigDB signature IDs for PTM-SEA input. `n_top` breaks p-value ties by the absolute
  statistic, so PTM-SEA's permutation-floor p-values no longer pull in every tied set.

- `data-raw/google_cloud_bucket_checks/` is removed; the bucket validation pipeline lives in
  motrpac-human-presuspension-repro.

## Documentation

- The 2.0.0 entry below gains the collection-level provenance behind the v2.0 regeneration
  and the versioning policy that governs it. Those objects have shipped since 2.0.0; only
  the record of them is new.

# MotrpacHumanPreSuspensionAnalysis 2.0.6

## Changes

- `plot_feature_heatmap()`: `multi_tissue_clust_rows = TRUE` works (it referenced an undefined
  object); new `right_annotation`, `heatmap_args`, `draw_args` and `return_drawing` arguments.

## Data objects

- `COVARIATES_FILE` drops the `BMI` and `codedsiteid` rows for `epigen-methylcap-seq`, which
  the MALAX model (`~0 + group_timepoint + age + sex + (1 | pid)`) never included.

# MotrpacHumanPreSuspensionAnalysis 2.0.5

## New data

- `PTMSEA_INPUT` — the prot-ph differential-analysis z-statistics as a PTM-SEA input, one
  `GCT` per tissue, from the confidently localized sites. `GCT` is `cmapR`'s S4 class;
  attach `cmapR` to use its accessors.

## Data objects

- The thirteen `*_DA` objects gain `CI.L_calculated` and `CI.R_calculated`, the 95%
  confidence interval on `logFC`, computed against each contrast's own residual degrees of
  freedom. No existing column or row changed. `topTable`'s `CI.L`/`CI.R` remain unshipped:
  for a `dream` fit they bound every contrast by the first contrast's degrees of freedom.

- `UTORONTO_TFs` is now the prot-ph TF regulator pool — 1,381 rows, `feature_id` and
  `gene_symbol` — where it was the raw UToronto extract, 2,765 rows and 28 columns keyed by
  Ensembl gene ID. The annotation columns are gone and a row is a phosphosite,
  not a gene.

# MotrpacHumanPreSuspensionAnalysis 2.0.4

## Data objects

- `FCM_CLUSTERS` clusters blood into 12 clusters where it had 13; adipose (13) and muscle
  (12) are unchanged and reproduce bit for bit. `FCM_CAMERA` and `FCM_ORA` follow the new
  clustering and have 434,337 rows where they had 444,806. Cluster numbers are a trajectory
  order, not identities, so blood cluster N in 2.0.4 is not blood cluster N in 2.0.3 —
  against the outgoing objects blood scores an adjusted Rand index of 0.711.

# MotrpacHumanPreSuspensionAnalysis 2.0.3

## Data objects

- `HUMAN_FEATURE_TO_GENE` gains three columns and has 13 where it had 10. No row and no
  existing column changed: all 1,920,618 rows and all ten previously shipped columns are
  identical to 2.0.2, so nothing that reads the table today reads anything different. Code
  selecting columns by position has to be updated; code selecting by name does not.

  - `custom_annotation` (factor) and `relationship_to_gene` (numeric) say where an
    ATAC-seq or MethylCap-seq peak sits relative to the gene it was assigned to — the
    region it falls in (`"Promoter (<=1kb)"`, `"Intron"`, `"Distal Intergenic"`, and seven
    others) and the signed base-pair distance to that gene, `0` where the peak overlaps it.
    Both were produced by the pipeline all along and dropped before the table was built, so
    an epigenomics feature arrived carrying only a gene: a peak in a promoter and a peak
    40 kb into an intron were indistinguishable once mapped. They are populated for all
    1,852,716 epigenomics rows and `NA` everywhere else.

    Unlike `confident_site` they are not a per-tissue measurement — they are derived from
    the peak coordinates in the `feature_id` — so they take the same value in every tissue
    a peak appears in and are not collapsed. This was checked rather than assumed: none of
    the 306,788 ATAC or 1,545,930 MethylCap feature identifiers shared between tissues
    disagree.

  - `confident_site` (logical) is the phosphosite localization flag, `NA` outside
    `prot-ph`. It is **collapsed across tissues** — `TRUE` only where a site is confidently
    localized in every tissue that measured it — because this table is keyed on
    `(assay, feature_id)` and has no tissue column, and 859 sites disagree between muscle
    and adipose. Read `*_PROT_PH_QC$feature_metadata` in
    `MotrpacHumanPreSuspensionData` when tissue-specific localization matters;
    `preprocess_PTMSEA()` already does, and is unaffected by this addition.
    **Removed again in 2.0.7** for the reason given there: a collapse across tissues is not
    the measurement for either tissue.

## Documentation

- `?HUMAN_FEATURE_TO_GENE` documented `assay` as a factor. It is a character vector, and
  has been for as long as the table has been built this way.

- The documented `assay` values did not include `"prot-clinical"`, which the 2.0 split of
  clinical chemistry into a metabolomics and a proteomics assay introduced. All eight
  values the column actually takes are now listed.

# MotrpacHumanPreSuspensionAnalysis 2.0.2

## Changes

- `load_differential_analysis(epigen = TRUE)` reads the epigenomics tables from Google
  Cloud Storage via gsutil instead of the public CloudFront release, and defaults to the
  current motrpac-human-presuspension-repro staging bucket. It now requires
  `repo_local_dir`, and gains `gsutil` and `bucket` arguments; `plot_single_feature()`
  passes all three through.
- The unexported `load_DA_from_AWS()` and `.load_single_ome_tissue_AWS()` are removed.
  Its pinned `version = "1.2"` no longer matched the atac-seq tables, which are at v2.0.

# MotrpacHumanPreSuspensionAnalysis 2.0.1

## Data objects

- The `*_SUM_STATS` objects are named, ordered and keyed the way the `*_DA` objects are.
  Every metabolomics platform is now labelled `assay = "metab"` with the platform
  in its own `platform` column, where before the platform was written into `assay` and
  there was no `platform` column. 

  There are 17 objects where there were 46. The research metabolomics platforms are no
  longer one object each: they are stacked into a single `{TISSUE}_METAB_SUM_STATS` per
  tissue 

- `load_summary_stats()` returns the research metabolomics platforms as a single `"metab"`
  element per tissue rather than one element per platform — the nesting
  `load_differential_analysis()` returns, so the two tiers can be walked together. Naming
  one platform still loads them all, as before.

## Bug fixes

- `plot_single_feature()` draws the same legend for every tissue, whether or not that
  tissue has a timepoint below the p threshold. Combining plots with
  `patchwork::plot_layout(guides = "collect")` previously produced a repeated p
  threshold legend, because collection only merges guides that are identical and a
  tissue with nothing significant contributed a one-key legend. A single plot with no
  significant timepoints now shows both p threshold keys rather than only `adj p >=`.

- `plot_single_feature()` plots clinical chemistry only when a clinical ome is
  requested. It previously appended the clinical rows whenever the feature name
  matched an analyte, so `selected_omes = "transcript-rna-seq"` with `"Glucose"`
  returned a clinical chemistry plot. `selected_omes` now accepts the omes in
  `clinical_ome_list()` by name and `"all"` includes them, matching how
  `load_differential_analysis()` treats clinical chemistry.

  A request that names another ome no longer returns clinical chemistry
  alongside it, and `"metab"` no longer implies `"metab-t-clinical"`. Analytes
  measured both clinically and on a research platform, such as Cortisol and
  Lactate, return only what was asked for.

- `plot_single_feature()` loads every non-epigenetic tissue and ome once and filters
  afterwards, rather than assembling the request ome by ome. Clinical chemistry is no
  longer a special case appended after the load, and the differential analysis and the
  summary statistics are put in one vocabulary before either is filtered.

  `plot_single_feature("VEGFA")` works again. The default `selected_omes = "all"` was
  broken for every non-metab feature by `filter(platform != "metab-t-conv")`: `platform`
  is NA on non-metab rows and `filter` drops NA, so the filter deleted the whole
  non-metab payload and the feature was reported as absent from the data.

- `plot_single_feature()` no longer excludes the `metab-t-conv` platform, which is now
  plotted and labelled `Conv. Metab (log2)`. It has no `assay_codes` row, so without
  that fallback its facet strip read `NA`. Note that it is the `metab-t-clinical`
  measurement on a log2 scale, so Glucose, Glycerol, KET and NEFA in blood now return a
  panel from each.


# MotrpacHumanPreSuspensionAnalysis 2.0.0

Data objects regenerated by the motrpac-human-presuspension-repro pipeline for the v2.0 data
collection.

## Versioning

**Package versions and data collection versions are not the same thing.** This package is
versioned independently of the freeze it carries: the 2.0.x series all ships the v2.0
collection, and a package release may change nothing about the data at all. Cite the
collection version, not the package version, when describing which data an analysis used.

Within the collection, versioning is per file: a file is bumped to v2.0 only where its
content actually changed, so a v2.0 collection legitimately contains files carrying earlier
version suffixes. See the MoTrPAC Knowledge Center for the release and versioning policy:
<https://motrpac-data.org/knowledge-center>.

## Why the v2.0 objects differ from v1.3

- **Sample misalignment in QC-norm batch correction, affecting transcriptomics and Olink.**
  `limma::removeBatchEffect()` pairs covariate row *i* with matrix column *i* positionally.
  The covariate table was built with `merge()`, which returns rows sorted by `vialLabel`,
  and was passed against a matrix whose columns were in count-file order (transcriptomics)
  or pivot order (Olink). The two orders are not the same, so each sample was
  batch-corrected using another sample's batch, site and plate assignment. Median
  per-feature correlation against v1.3 is 0.925 for transcriptomics (0.908-0.944 across the
  three tissues) and 0.900 for `prot-ol`. Every differential-analysis and summary-statistic
  object on those two omes moves with it: counted at each feature's best post-exercise
  timepoint on the endurance exercise-vs-control contrasts, blood `prot-ol` goes from 103
  to 146 significant features of 1,417 Olink targets.

- **Replicate averaging in muscle proteomics (`prot-pr`, `prot-ph`).** Muscle samples
  measured twice are meant to be merged by averaging the pair and dropping the now-redundant
  column. For intra-site pairs the mean was written into the column that was about to be
  deleted, so the value that survived was the first measurement on its own rather than the
  mean of the two. Inter-site pairs were averaged correctly, which is what made this easy to
  miss. Batch correction and replicate handling also ran in a different order in `prot-pr`
  than in `prot-ph`, so the two omes were not processed identically.

- **Feature metadata rebuilt.** Each ome's `metadata_features` is now a self-contained
  feature-to-gene mapping that matches its QC-norm matrix exactly. This is why every `*_QC`
  object in `MotrpacHumanPreSuspensionData` differs from v1.3 while only six have a
  `qc_norm` matrix whose values differ.

- **ATAC differential analysis refit.** The released DA tables were built against an earlier
  feature set than the QC-norm matrices they accompany; the two are reconciled here. Muscle
  ATAC moves from 1,584 to 1,521 significant features on the delta-delta contrasts — 87
  gained, 149 lost, `logFC` correlating 0.992 with the release. Blood (`t05-pbmc`) has no
  significant peaks, as in v1.3.

- **Six ATAC samples removed as sample mix-ups.** `OUTLIERS` has 160 rows where v1.3 had
  154; the six added are `epigen-atac-seq` samples, four blood and two muscle, now excluded
  from the analysis. No row was dropped.

- **One muscle `transcript-rna-seq` sample restored.** It was missing from
  `metadata_samples` and was therefore dropped at the modelling stage. Differences in the
  refit differential analysis are minimal.

- **Metabolite names follow a frozen RefMet snapshot with updated overrides.** 78
  metabolomics `feature_id`s in the `*_METAB_DA` objects are renamed from v1.3, so joins on
  v1.3 names will miss them. The curated overrides gain `13 HODE` -> `13-HODE` and `20-HETE`
  entries that previously resolved to `NA` (adipose `13-HODE` now enters the DA). Other
  renames: `*` suffixes dropped (`FA 18:1;O*` -> `FA 18:1;O`), `NAGly` -> `NA-Gly`, `NAD+`
  -> `NAD`, `DHA` -> `Docosahexaenoic acid`, `Edetic acid` -> `EDTA`, `13-Oxo-ODE` ->
  `13-OxoODE`.

## Reference and enrichment objects

- `CAMERA_RESULTS` has 1,016,991 rows, 1,869 fewer than v1.3, the losses concentrated in
  blood `prot-ol` (-1,419), adipose `prot-pr` (-270) and muscle `prot-pr` (-189): the gene
  universe moved with the feature-metadata rebuild, and a set that no longer intersects it
  is not scored. The ranking is largely intact — Spearman 0.957 to
  1.000 on the signed statistic.

- `METABOLOMICS_CVS` has 4,000 rows where the v1.3 object had 3,878, and the 122 rows are
  mostly a v1.3 defect rather than a v2.0 change: the v1.3 `.rda` was staler than the v1.3
  `.txt` it was meant to mirror, which already published 4,000 keys. The change that
  actually propagates downstream is 69 `lowest_CV` flips and 93 `refmet_name` corrections,
  which together decide which copy of a duplicated RefMet name survives de-duplication.

- `OME_TISSUE_CODE` has 53 rows where it had 49: it gains the two clinical omes
  (`metab-t-clinical`, `prot-clinical`) and the five `lab-*` plasma-chemistry tiers, and
  loses `metab-meta-reg` in all three tissues. The `lab-*` codes are inputs to the clinical
  omes and are deliberately absent from `ome_available_list()`; do not pass them to a
  loader.

- `"metab-meta-reg"` is gone from `ome_available_list()`, `OME_TISSUE_CODE`
  and `HUMAN_OME_COLORS` (31 entries where there were 32). No object was ever built under
  that name, so every accessor offered it as a choice that returned nothing. Passing it to
  a loader now errors at `match.arg()`, and `HUMAN_OME_COLORS[["metab-meta-reg"]]` returns
  `NULL`.

## Split data

- `BLOOD_CLINICAL_CHEMISTRY_SUM_STATS` — replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_SUM_STATS and BLOOD_PROT_CLINICAL_SUM_STATS.
- `BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS` — no blood ATAC feature is significant at adj_p_value < 0.05 this cycle, and epigenomics summary statistics carry significant features only, so no object is built.
- `CLIN_CHEMISTRY_DA` — replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_DA and BLOOD_PROT_CLINICAL_DA.

## New data

- BLOOD_METAB_T_CLINICAL_DA, BLOOD_METAB_T_CLINICAL_SUM_STATS, BLOOD_PROT_CLINICAL_DA, BLOOD_PROT_CLINICAL_SUM_STATS

  Clinical chemistry, one assay in v1.3, is split into a metabolomics and a proteomics assay.

## Changes to data objects

- 11 objects drop the `CI.L`, `CI.R` columns; code that selects them will error: ADIPOSE_METAB_DA, ADIPOSE_PROT_PH_DA, ADIPOSE_PROT_PR_DA, ADIPOSE_TRNSCRPT_DA, BLOOD_METAB_DA, BLOOD_PROT_OL_DA, BLOOD_TRNSCRPT_DA, MUSCLE_METAB_DA, MUSCLE_PROT_PH_DA, MUSCLE_PROT_PR_DA, MUSCLE_TRNSCRPT_DA.

- HUMAN_FEATURE_TO_GENE gains the `flanking_sequence` column.


# MotrpacHumanPreSuspensionAnalysis 0.2.4

This release addresses installation failures on R 4.6 and removes the `Remotes:`
field that was breaking dependency resolution.

Verified with clean installs on R 4.4 (Bioconductor 3.20), R 4.5 (3.22), and
R 4.6 (3.23).

# MotrpacHumanPreSuspensionAnalysis 0.2.3

## User-facing functions

- clinical chemistry is labeled `clinical-chemistry` instead of `clinical_chemistry` 
- `plot_single_feature()` now clarifies clinical chemistry features correctly; the y-axis
  label now clarifies that clinical chemistry features are shown on an absolute scale
  while all other features are shown as `log2(normalized value)`.
- Updated `HUMAN_FEATURE_TO_GENE` following changes described in 0.2.2., which now maps to the features in qc-norm properly.   

## Internals and package checks

- Namespace hygiene updates - proper imports have been labeled throughout
- Trimmed `globalVariables()`, consolidated `@importFrom` tags, added the `grid.rect`
  import, filled in the package `Description`, and added `RColorBrewer`, `doParallel`,
  `parallel`, and `randomForest` to `Suggests`. 

## Backend: data generation and release tooling (data-raw)

These changes only affect analysts with sample-level data from `MotrpacHumanPreSuspensionData`.

- Added a Google Cloud bucket validation pipeline under
  `data-raw/google_cloud_bucket_checks/`. The numbered scripts snapshot the production
  GCS bucket, copy it to an isolated staging folder, diff staging against the local
  updated files, upload changed files (removing superseded versions), and validate the
  staging structure and values against the installed package. A `run_validation_pipeline.R`
  driver orchestrates the steps; see the directory `README.md` for setup.
- Moved ATAC peak annotation into `generate_atac_qc_norm.R`
  (`.annotate_atac_features()` / `pre_cawg_get_peak_annotations_hs()`) and refined the
  gene-mapping logic. Also updated the methylcap, clinical, and metabolomics generation
  scripts.


# MotrpacHumanPreSuspensionAnalysis 0.2.2


## Backend: Feature metadata gene annotation (data-raw)

These changes will only make functional differences for analysts with sample level data from `MotrpacHumanPreSuspensionData`

- Added gene-level annotation to feature_metadata outputs for all proteomics omes. Previously,
  feature_metadata for Prot-PR, Prot-PH, and Prot-OL contained only raw provenance columns (UniProt
  accessions, PTM identifiers, redundant IDs) with no standardized gene symbol or Ensembl mapping.
  Feature metadata files now include `gene_symbol`, `ensembl_gene`, and `entrez_gene` columns
  derived from a three-round BioMart lookup strategy.
- Added `.annotate_olink()` to `generate_prot_ol_qc_norm.R`: resolves Olink protein metadata
  (UniProt accession from `uniprot_entry`) to gene symbols and Ensembl IDs. Lookup proceeds via
  UniProt → gene symbol → Entrez ID, with three progressive BioMart passes to maximize coverage.
- Added `.annotate_prot_pr()` to `generate_prot_pr_qc_norm.R`: extends the per-tissue
  `feature_metadata_output` from `PR@rdesc` with gene annotations. Uses `feature_id` (= `protein_id`,
  a UniProt accession) taken directly from the file as the BioMart lookup key.
- Added `.annotate_prot_ph()` to `generate_prot_ph_qc_norm.R`: same annotation for
  phosphoproteomics. Uses `protein_id` (not `feature_id`, which is the PTM site identifier) taken
  directly from the file as the UniProt lookup key, preserving the full PTM `feature_id` in the output.
- Prot-PR and Prot-PH annotation now resolves UniProt accessions to genes via a single BioMart
  `getBM` lookup (mirroring the Olink `.annotate_olink()` implementation), using the UniProt
  accession read directly from the source files rather than a downloaded UniProt ID mapping table.


# MotrpacHumanPreSuspensionAnalysis 0.2.1

- Added clinical chemistry differential analysis using the same structure as the molecular differential analysis;
  see: `CLIN_CHEMISTRY_DA`

# MotrpacHumanPreSuspensionAnalysis 0.2.0

## Documentation and website

- Added pkgdown GitHub Actions workflow for automated website deployment.
- Added pkgdown configuration updates for GitHub Pages publishing.
- Clarified consortium-only optional functionality in the README.

## Dependency handling

- Made `MotrpacHumanPreSuspensionData` an optional consortium-only runtime dependency
  instead of a public package dependency for CI/pkgdown resolution.

# MotrpacHumanPreSuspensionAnalysis 0.1.0

## Initial release

- Public release of summary statistics and modeling outputs from the
  MoTrPAC human pre-COVID suspension phase.
- Added functions for loading differential analysis and summary statistics.
- Added clustering, enrichment, and visualization utilities.
- Added vignettes and package website scaffolding.
