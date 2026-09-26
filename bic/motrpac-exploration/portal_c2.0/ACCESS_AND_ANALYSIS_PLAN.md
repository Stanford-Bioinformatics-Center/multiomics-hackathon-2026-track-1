# Human acute-exercise collection: access and analysis audit

Verified 2026-09-25.

## Exact study and release

- Portal title: **Acute Exercise in Human Sedentary Adults**.
- Study identifier: `human-precovid-sed-adu`.
- Current portal analysis collection: **c2.0**, with **131 files** visible in the file browser.
- Cohort: 175 pre-suspension sedentary adults, randomized to endurance exercise, resistance exercise, or non-exercising control; acute sampling in blood/plasma, skeletal muscle, and adipose through 24 hours, with timing and coverage varying by tissue and assay.
- Portal: https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0

The existing Table 3 analysis uses `MotrpacHumanPreSuspensionAnalysis` **2.0.8**, pinned to commit `535b4044e7417413de471104c619120337602b77`. The official package NEWS explicitly says all 2.0.x package releases carry the v2.0 data collection. Thus the current screen is already based on this human study and collection. Package version, collection version, and individual file suffixes are distinct identifiers. Numerical equality to freshly downloaded portal tables has not been independently checked.

Package provenance: https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/NEWS.md

## Actual access boundary

Public summary results support effect sizes, uncertainty, significance, time-course comparisons, gene/pathway overlays, and comparisons with independent disease signatures. The portal labels human Quant-ID and Phenotype collections restricted and points to dbGaP. Individual-level normalized matrices are also restricted according to the release documentation. A generic page description saying Analysis includes normalized tables does not make human participant-level tables public.

The portal download button opens a name/email request form; this was inspected and closed without submission. Existing public R-package files remain available locally and usable without submitting that form. The release-notes PDF was downloaded directly and saved here; it still describes v1.3 and is outdated as a guide to the latest collection version. Use the live portal and current package NEWS for version selection.

## Confirmed public file classes

- RNA-seq differential results in blood, muscle, and adipose.
- Plasma Olink protein differential results and feature annotations.
- Global proteomics and phosphoproteomics in muscle and adipose.
- Targeted and untargeted metabolomics, including per-platform files and clinical chemistry.
- ATAC-seq results in PBMCs and muscle; MethylCap-seq results in blood, muscle, and adipose.
- Feature-to-gene mapping, RefMet mapping/provenance, metabolite CV annotations, feature metadata, and assay QC reports.

Example filenames actually visible in c2.0:

```
human-precovid-sed-adu_t02-plasma_prot-ol_da_dream-acute_v2.1.txt
human-precovid-sed-adu_t02-plasma_prot-ol_metadata_features_v2.0.txt
human-precovid-sed-adu_t02-plasma_metab-t-oxylipneg_da_dream-acute_v2.1.txt
human-precovid-sed-adu_t02-plasma_metab-t-amines_metadata_features_v2.0.txt
human-precovid-sed-adu_t06-muscle_transcript-rna-seq_da_dream-acute_v2.1.txt
human-precovid-sed-adu_t06-muscle_epigen-atac-seq_da_dream-acute_v2.1.txt
motrpac-mappings-human-feature-to-gene_v2.2.txt
motrpac_human-precovid_refmet-map_provenance_v2.0.json
```

## Recommended project using available data

**Question:** Which of the 28 cardiometabolic exerkines in Chow et al. Table 3 show an acute response in humans, and which responding molecules, receptors, or downstream pathways connect to type 2 diabetes through independent human evidence?

1. Retain all 28 candidates, including non-significant, absent, and ambiguous identities. Use exercise-versus-control difference-in-changes estimates; retain direct endurance-versus-resistance contrasts separately.
2. Use plasma results to establish a circulating response; use tissue RNA, protein, and phosphosite results as biological context. Do not infer secretion or tissue origin from RNA alone.
3. Extend identity checking with per-platform metabolite annotations/results. The existing package metabolomics screen uses the consortium's lowest-CV selected measurement per metabolite; it is not an exhaustive interrogation of every platform measurement. Keep BAIBA unresolved until its conflicting chemical identifiers are reconciled, and distinguish 12,13-diHOME from 9,10-diHOME.
4. Overlay a versioned human T2D gene/variant evidence resource and a separately chosen disease expression dataset. Specify the tested feature universe, disease resource version, overlap/enrichment test, and multiple-testing family before reporting enrichment. Receptor and downstream-pathway evidence should be labeled separately from evidence for the exerkine gene itself.
5. Prioritize transparent hypotheses showing the exercise effect, disease evidence, directional compatibility or conflict, and a falsifiable follow-up. An acute rise does not itself imply benefit; opposite disease/exercise expression does not establish causal reversal.

Already completed: the 28-entry lookup, assay coverage audit, 1,171 candidate/assay/contrast result rows, and 13 differential-result objects searched. See `../table3/REPORT.md`. Sixteen candidates have resolved circulating matches; CX3CL1 and lactate have at least one exercise-versus-control result passing the original source FDR threshold. These are current screening results, not established mediators of disease benefit.

Not yet performed: independent disease-signature comparison, disease genetic enrichment/colocalization, exhaustive portal per-platform metabolite re-screen, or epigenomic integration. Subject-level correlation, response heterogeneity, mediation, or new covariate-adjusted models would require approved individual-level data.

## Sources and saved records

- Portal access overview: https://motrpac-data.org/data-download
- Public package overview: https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis
- Saved `package_NEWS.md`: pinned version history with explicit v2.0 collection provenance.
- Saved `release_notes.pdf` and `release_notes.txt`: portal-linked document, useful for study design/access levels but still labeled v1.3.
