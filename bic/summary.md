# MoTrPAC exercise–disease exploration

Work completed on 25–26 September 2026.

The project now lives in `/Users/acheron/BIC`. Start with the executed [fractalkine differential-analysis notebook](fractalkine_differential_analysis.ipynb); the original exploration remains in `motrpac-exploration/`. The supplied review PDF is also available locally in `references/`.

## Purpose and current conclusion

Explore relationships between exercise-responsive human biology and disease using data that are actually available, and generate transparent, testable hypotheses for the [Stanford challenge](https://stanford.bioinformatics-center.org/).

The work progressed from an initial disease-hypothesis shortlist to a systematic screen of the 28 cardiometabolic exerkines in Table 3 of Chow et al., followed by a detailed significance audit of the 16 candidates with identifiable plasma measurements.

**Main finding:** Of those 16 candidates, **fractalkine/CX3CL1 and lactate** have exercise-versus-control plasma responses passing the original MoTrPAC adjusted p-value threshold of 0.05. **FGF21 and IL-15** have secondary, exploratory evidence under a separate correction restricted to the candidate tests. Several candidates with nonsignificant plasma results have significant tissue RNA responses.

CX3CL1 is the leading circulating protein candidate for a subsequent disease-evidence analysis. Type 2 diabetes is the proposed disease focus, but a systematic disease-genetics or disease-expression integration has not yet been performed.

## Data used and access verified

- **Study:** Acute Exercise in Human Sedentary Adults.
- **Study identifier:** `human-precovid-sed-adu`.
- **Collection:** c2.0 / v2.0, as verified on the portal on 25 September 2026.
- **Cohort:** 175 pre-suspension sedentary adults randomized to endurance exercise (EE), resistance exercise (RE), or non-exercise control (CON). Sample availability varies by tissue, assay, participant, and timepoint.
- **Sampling:** Blood/plasma, skeletal muscle, and adipose, around a single acute exercise bout and through 24 hours afterward. There are no during-resistance samples.
- **Public source used for calculations:** `MotrpacHumanPreSuspensionAnalysis` version **2.0.8**, pinned to commit `535b4044e7417413de471104c619120337602b77`.
- **Portal inventory:** 131 public analysis and annotation files in c2.0. The package version history explicitly states that its 2.0.x series carries the v2.0 data collection. Package versions, collection versions, and individual file suffixes are distinct.

The calculations use public consortium summary statistics and differential-analysis results. They do not use individual participants' omics matrices or phenotype records. The portal directs access requests for restricted human measurements and phenotypes to dbGaP. No access application or name/email download request was submitted.

The portal-linked release-notes PDF still describes v1.3; current collection identification was checked against the live portal and current package version history. Numerical equality between the package objects and newly downloaded portal tables has not been independently tested.

Sources: [portal collection](https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0), [public package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis), [pinned version history](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/NEWS.md).

## Candidate sources and screening method

Candidate nomination used [Exerkine Atlas](https://exerkineatlas.org/) and Chow et al., [*Exerkines in health, resilience and disease*](https://www.nature.com/articles/s41574-022-00641-2). The supplied local PDF was read, and both pages of Table 3 were visually checked. The table contains 28 examples affecting the cardiometabolic system. Its arrows describe plasma changes and combine human, animal, and cell evidence; not every listed action is established in humans.

All 28 entries were retained, including negative and unresolved results. Gene aliases and metabolite identities were mapped explicitly. Thirteen public differential-result objects were searched: RNA in blood, muscle, and adipose; plasma Olink proteins; muscle/adipose proteins and phosphoproteins; plasma/muscle/adipose metabolomics; and available blood clinical chemistry/protein results. Epigenomic integration was not performed.

The primary comparisons are **differences in changes**: change from baseline in an exercise arm minus the corresponding change in resting controls. Direct EE–RE comparisons were retained separately. No new participant-level statistical model was fitted.

Time labels were translated by tissue:

- `post_15_30_45_min`: muscle 15 minutes, plasma/blood 30 minutes, adipose 45 minutes.
- `post_3.5_4_hr`: muscle/plasma/blood 3.5 hours, adipose 4 hours.

## Coverage of the 28 Table 3 candidates

- **16** had resolved matches in the analyzed plasma protein/metabolite tables.
- **11** had no match in those analyzed plasma tables.
- **1**, BAIBA, had unresolved chemical identity.

“Circulating measurement” means an assay of the molecule in **blood plasma**. It is distinct from measuring that molecule's gene expression in blood cells, muscle, or adipose. The public results summarize assays performed on participant specimens; the underlying individual measurements were not accessed.

The 16 candidates comprise **15 protein candidates represented by 21 Olink assay IDs**, plus **lactate measured by two assays**. IL-6 and IL-8 each have four Olink IDs; these were retained separately rather than merged or selected for significance.

Important identity limitations:

- **BAIBA:** raw beta-aminoisobutyric-acid labels, generic RefMet naming, and a KEGG annotation identifying alpha/2-aminoisobutyric acid conflict. BAIBA was excluded from confirmed plasma matches and significance counts.
- **12,13-diHOME:** no exact analyzed match was found; 9,10-diHOME was not substituted.
- **Irisin:** FNDC5 RNA or precursor protein does not establish mature circulating irisin.
- **HSP72:** the plasma match is HSPA1A. Tissue context also includes explicitly labeled HSPA1B measurements.
- **VEGF:** operationally mapped to VEGFA, not the whole VEGF family.
- Package metabolomics results select the lowest-CV measurement for a mapped metabolite. An exhaustive re-screen of every portal platform file remains outstanding.

Absence from an analyzed table can reflect assay coverage, QC, or analysis eligibility. It does not establish biological absence.

## Significance of the 16 plasma candidates

The audit covered **230 exercise-versus-control plasma tests**: 23 assay features, each with six EE–CON and four RE–CON comparisons. It also checked **92 direct EE–RE plasma comparisons**.

The table below shows each candidate's **smallest original MoTrPAC adjusted p-value** across tested conditions. These minima are descriptive, not candidate-level omnibus p-values. Directions indicate model estimates, including nonsignificant estimates. “After” means after the end of the exercise bout.

| Candidate | Comparison with resting controls | Estimated direction | Original adjusted p |
|---|---|---|---:|
| **Lactate** | RE, 10 min after | Increase | **3.62 × 10⁻⁸¹** |
| **Fractalkine / CX3CL1** | EE, 40 min during | Increase | **0.000227** |
| FGF21 | RE, 3.5 h after | Decrease | 0.0714 |
| IL-15 | RE, 30 min after | Decrease | 0.0717 |
| VEGF / VEGFA | EE, 40 min during | Increase | 0.101 |
| IL-8 / CXCL8 | EE, 40 min during | Increase | 0.127 |
| GDF15 | EE, 40 min during | Increase | 0.157 |
| ANGPTL4 | RE, 30 min after | Decrease | 0.164 |
| IL-6 | RE, 3.5 h after | Decrease | 0.177 |
| SPARC | RE, 30 min after | Decrease | 0.184 |
| Follistatin / FST | RE, 24 h after | Decrease | 0.220 |
| Angiopoietin-1 / ANGPT1 | RE, 30 min after | Decrease | 0.263 |
| HSP72 / HSPA1A | RE, 10 min after | Increase | 0.552 |
| IL-7 | EE, 20 min during | Increase | 0.555 |
| SDC4 | EE, 40 min during | Increase | 0.588 |
| Myostatin / MSTN | EE, 40 min during | Increase | 0.632 |

### Primary findings

**CX3CL1:** Plasma increases during endurance at both 20 minutes (adjusted p = 0.000364) and 40 minutes (0.000227). At 40 minutes, the source model effect is +0.561, with pointwise 95% CI [0.332, 0.789]. No direct post-exercise EE–RE plasma contrast is significant; the smallest adjusted p is 0.472. The lack of during-RE sampling prevents an endurance-specific claim about the circulating pulse.

**Lactate:** Both clinical chemistry and research metabolomics show increases during endurance at 20/40 minutes and at 10/30 minutes afterward, and at 10/30 minutes after resistance. The minimum adjusted p above comes from clinical chemistry; the research assay minimum is 2.83 × 10⁻⁶⁴. Direct EE–RE contrasts support a larger increase after resistance at 10 and 30 minutes. The assays corroborate the response within this cohort, not in independent cohorts.

There are **14 source-FDR-significant plasma test rows**, representing these two candidates: two CX3CL1 tests and twelve lactate tests across two assays.

### Secondary findings and correction families

Original MoTrPAC adjusted p-values use Benjamini–Hochberg correction within tissue, assay, platform, and contrast. They do not provide a separate omnibus correction over all times and assays in this exploratory screen.

A separate BH calculation over all **230 resolved plasma candidate tests** reproduces the earlier Table 3 analysis. Under this smaller test family, **FGF21 and IL-15 each have adjusted p ≈ 0.0435** for the decreases listed above. Their raw p-values are approximately 0.003. They remain exploratory findings because they fail the original source correction and the data have already been examined. Follistatin's secondary adjusted p is approximately 0.054 and does not pass 0.05.

Pointwise confidence intervals excluding zero can coexist with nonsignificant adjusted p-values. A nonsignificant result is not an equivalence test or proof of no response.

## Tissue responses

Strong tissue responses occur even when plasma results are nonsignificant. Selected muscle RNA results are:

| Gene | Response versus controls | Original adjusted p |
|---|---|---:|
| CX3CL1 | Increase, 15 min after EE | 5.50 × 10⁻⁵³ |
| SDC4 | Increase, 3.5 h after RE | 3.03 × 10⁻³² |
| VEGFA | Increase, 3.5 h after EE | 6.64 × 10⁻²⁷ |
| HSPA1A | Increase, 3.5 h after RE | 1.04 × 10⁻¹² |
| MSTN | Decrease, 3.5 h after RE | 1.43 × 10⁻¹⁰ |
| IL15 | Decrease, 3.5 h after RE | 7.66 × 10⁻⁵ |
| ANGPTL4 | Increase, 3.5 h after EE | 0.00135 |

CX3CL1 also increases in muscle RNA after resistance and in adipose RNA after endurance. ANGPTL4 has an earlier muscle RNA decrease before its later increase. Blood-cell RNA responses and protein/phosphosite context are retained in the detailed report. These source-FDR results did not receive an additional cross-tissue omnibus correction.

RNA responses do not establish translation, secretion, source tissue, or target-organ activity. The early CX3CL1 plasma increase precedes the first post-exercise muscle biopsy; new muscle transcription cannot be assumed to cause that initial pulse.

## Disease hypotheses considered

1. **CX3CL1 and type 2 diabetes:** The exercise-associated pulse may connect to islet survival or glucagon regulation. External human-islet experiments reported reduced beta-cell apoptosis and glucagon secretion, without increased human-islet glucose-stimulated insulin secretion. Exercise-related islet benefit and the source of circulating CX3CL1 remain untested. [Primary human-islet study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4209359/).
2. **BCAAs and diabetes-related biomarkers:** The initial broader screen found plasma leucine/isoleucine/valine decreases after resistance exercise. These could reflect uptake, redistribution, synthesis, or oxidation. Concentrations and later pathway RNA changes do not establish flux or clinical benefit.
3. **Kynurenine handling and depression:** The broader screen found an early endurance-associated rise in kynurenic acid and fall in kynurenine. This provides a possible connection to prior muscle–kynurenine experiments, but no brain or depression outcomes were measured here. [Mechanistic study](https://pubmed.ncbi.nlm.nih.gov/25259918/).
4. **Neprilysin and heart failure:** A plasma MME increase was used as an interpretive counterexample: an exercise-associated abundance increase does not necessarily imply that increasing the protein's activity benefits disease.

BCAAs, kynurenines, and MME were part of the initial broader exploration, not the subsequent 16-candidate Table 3 significance audit. The initial report contains their quantitative results and supporting references.

## Verification and deliverables

The full Table 3 screen produced **1,171 candidate/assay/contrast rows**. For the 16-candidate follow-up, **897 extracted rows** were independently matched against the pinned R objects, verifying effects, confidence intervals, raw p-values, and original adjusted p-values. Python checks confirmed unique keys, p-value ranges, interval bounds, contrast coverage, and the 230-test BH calculation. The significance figure was visually inspected.

All analysis files are under [motrpac-exploration](motrpac-exploration/).

| Deliverable | Location |
|---|---|
| Initial disease-hypothesis shortlist | [README.md](motrpac-exploration/README.md) |
| Dataset and access audit | [ACCESS_AND_ANALYSIS_PLAN.md](motrpac-exploration/portal_c2.0/ACCESS_AND_ANALYSIS_PLAN.md) |
| All 28 Table 3 candidates and findings | [Table 3 report](motrpac-exploration/table3/REPORT.md) |
| Candidate mapping and metadata | [table3_candidates.csv](motrpac-exploration/table3/table3_candidates.csv) |
| Coverage of all 28 candidates | [coverage_summary.csv](motrpac-exploration/table3/coverage_summary.csv) |
| Full Table 3 result rows | [all_candidate_results.csv](motrpac-exploration/table3/all_candidate_results.csv) |
| Plasma results and correction values | [plasma_results.csv](motrpac-exploration/table3/plasma_results.csv) |
| Detailed 16-candidate significance report | [REPORT.md](motrpac-exploration/sixteen_candidates/REPORT.md) |
| All 16-candidate statistics and provenance | [results.json](motrpac-exploration/sixteen_candidates/results.json) |
| All-timepoint plasma significance map | [plasma_significance.png](motrpac-exploration/sixteen_candidates/plasma_significance.png) |

Reproduction scripts: `reproduce_screen.R`, `plot_results.py`, `table3_manifest.py`, `table3_screen.R`, `table3_report.py`, `verify_sixteen.R`, and `explore_sixteen.py`, all in `motrpac-exploration/`. Input checksums and pinned source versions are saved with the results.

## Outstanding work

The focused CX3CL1 notebook now extracts all source contrasts from the cached human c2.0 data, independently reproduces the full plasma assay's BH corrections, checks a ten-test Holm sensitivity analysis, compares uncontrolled and control-adjusted responses, and includes direct exercise-mode tests for plasma and tissue RNA. It generates figures, tables, findings, and provenance in [results/fractalkine](results/fractalkine/).

The notebook confirms the two during-endurance plasma findings. It also identifies a stronger muscle RNA response after endurance than resistance at 15 minutes (direct EE–RE source adjusted p = 4.22 × 10⁻¹²), while direct post-exercise plasma differences remain nonsignificant. This is an analysis of published differential-model results, not a new participant-level model fit.

1. Integrate a versioned human disease gene/variant resource and an appropriate independent disease-expression dataset, with an explicitly defined tested feature universe and multiple-testing plan.
2. Prioritize CX3CL1 for the circulating-protein hypothesis, FGF21/IL-15 for exploratory follow-up, and MSTN/ANGPTL4/VEGFA/SDC4 for tissue-response hypotheses.
3. Revisit individual metabolomics platform annotations and results to resolve BAIBA and assess coverage beyond the package's selected metabolite measurements.
4. Add epigenomic evidence only where it addresses a specific regulatory hypothesis.
5. Define falsifiable validation experiments or independent datasets for shortlisted mechanisms.

No systematic disease-genetic enrichment, colocalization, disease-signature reversal analysis, participant-level correlation, mediation, clinical prediction, or new differential model has been completed. Individual response analyses would require approved participant-level data. The current work establishes exercise-response evidence and candidate priorities; it does not establish causal disease protection.

## CCN1 / CYR61 follow-up — 26 September 2026

Created a separate [CCN1 notebook](ccn1_differential_analysis.ipynb) using the same pinned **human c2.0 acute-exercise data**. The notebook explains CCN1's extracellular-matrix signaling role and reviews diabetes relevance before analyzing its exercise responses. CCN1 is already highlighted in the [MoTrPAC multi-tissue preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/); this analysis is not an independent discovery or replication.

The disease rationale includes exercise-responsive muscle CYR61 expression in people with T2D, circulating CCN1 associations with diabetic retinopathy, and experimental retinal injury mechanisms. These support investigation of diabetic complications and tissue remodeling, not a claim that CCN1 lowers glucose or mediates exercise's protection against diabetes. [Human exercise study](https://pubmed.ncbi.nlm.nih.gov/29924476/), [human retinopathy study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10273100/), [retinal mechanism study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11097549/).

The extraction retains **105 CCN1 source contrasts**. The focused analysis contains **30 exercise-versus-control** and **14 direct EE–RE** tests across muscle/adipose RNA, plasma Olink protein, and muscle/adipose protein. Full-assay BH recalculation reproduces the source values for **426,296 feature/contrast rows in 44 families**; the maximum absolute difference is approximately **3.37 × 10⁻¹³**. No CCN1 result is available in the analyzed blood RNA or phosphoprotein tables. Adipose protein has only a 4-hour recovery result.

At source BH-adjusted p < 0.05:

- **Muscle RNA:** positive control-adjusted effects after endurance at 15 minutes, and after resistance at 15 minutes, 3.5 hours, and 24 hours. Direct EE–RE tests favor a larger resistance response at 3.5 and 24 hours.
- **Adipose RNA:** positive effects after both exercise modes at 45 minutes; later results are nonsignificant.
- **Plasma CCN1:** positive effects during endurance at 40 minutes (adjusted p = 0.0222) and after resistance at 10 minutes (0.00102). None of the direct plasma mode comparisons passes source BH < 0.05.
- **Muscle protein:** positive resistance response at 3.5 hours (adjusted p = 4.67 × 10⁻⁹), with a direct difference favoring resistance. Adipose protein at 4 hours is nonsignificant.

Eight of the nine significant control-adjusted results also pass an exploratory Holm correction across the 30 CCN1 tests. The 24-hour resistance muscle RNA result narrowly misses this additional check (Holm p = 0.0506), although its source BH p is 0.0152. Three direct comparisons pass the smaller CCN1-only Holm family but not source BH; the notebook labels that distinction and retains source BH as the primary rule.

The notebook includes five figures, all available focal comparisons, source model formulas, checksums, and falsifiable diabetes-focused hypotheses. Results and provenance are in [results/ccn1](results/ccn1/), with a concise [findings report](results/ccn1/FINDINGS.md). All 12 code cells were executed successfully. MoTrPAC supplied the model estimates and raw p-values; this is a reanalysis of those published results, without participant-level model fitting, diabetes-group testing, or a secretion/causality claim.
