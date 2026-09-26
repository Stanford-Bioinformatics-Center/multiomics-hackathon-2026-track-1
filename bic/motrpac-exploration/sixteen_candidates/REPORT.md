# Significance of 16 Table 3 candidates in human acute exercise

Analysis date: 2026-09-26. MoTrPAC human pre-suspension sedentary adult collection c2.0; public Analysis package 2.0.8, commit `535b4044e7417413de471104c619120337602b77`.

## Findings

- **Fractalkine/CX3CL1 and lactate** have circulating responses passing the original MoTrPAC FDR threshold.
- **FGF21 and IL-15** have decreasing plasma responses passing a separate exploratory BH correction over these 230 candidate tests, but not the original source correction.
- **Lactate alone** has a source-FDR-significant direct endurance-versus-resistance plasma contrast among these 16 candidates.
- Several candidates with nonsignificant plasma results have tissue RNA responses. These are separate observations, not evidence that the tissue released the protein into circulation.

## What was tested

Fifteen protein candidates comprise 21 Olink assay IDs: IL-6 and IL-8 each have four IDs, retained separately. Lactate has a research metabolomics measurement and a clinical chemistry measurement. The resulting 23 assay features each have six endurance-versus-control and four resistance-versus-control tests: 230 tests. The same 23 features each have four direct endurance-versus-resistance tests: 92 tests. There are no during-resistance samples.

All effects are consortium difference-in-changes estimates: change from pre-exercise in an exercise arm minus the time-matched change in controls. EE = endurance; RE = resistance; CON = non-exercise control. We did not refit individual-level models.

## Strongest source-adjusted plasma result for every candidate

Each row selects the smallest original adjusted p-value across that candidate's tested assay/time/mode conditions. This is a descriptive minimum, not a candidate-level omnibus p-value. Directions for nonsignificant results describe estimates only. Effect and confidence interval use the native source model scale; do not compare effect magnitudes across the clinical lactate and other assays.

| Candidate | Condition | Effect [pointwise 95% CI] | Raw p | Source adjusted p | Candidate-family BH | Interpretation |
|---|---|---|---|---|---|---|
| Lactate | RE-CON, post 10 min | +2.686 [2.448, 2.925] | 4.02e-82 | 3.62e-81 | 9.25e-80 | Source FDR < 0.05 |
| Fractalkine | EE-CON, during 40 min | +0.561 [0.332, 0.789] | 2.89e-06 | 0.000227 | 4.75e-05 | Source FDR < 0.05 |
| FGF21 | RE-CON, post 3.5 h | -0.887 [-1.469, -0.305] | 0.00302 | 0.0714 | 0.0435 | Secondary family only |
| IL-15 | RE-CON, post 30 min | -0.271 [-0.448, -0.093] | 0.00298 | 0.0717 | 0.0435 | Secondary family only |
| VEGF | EE-CON, during 40 min | +0.369 [0.068, 0.671] | 0.0168 | 0.101 | 0.154 | No source-FDR hit |
| IL-8 | EE-CON, during 40 min | +0.344 [0.041, 0.646] | 0.0263 | 0.127 | 0.209 | No source-FDR hit |
| GDF15 | EE-CON, during 40 min | +0.360 [0.017, 0.702] | 0.0397 | 0.157 | 0.285 | No source-FDR hit |
| Angiopoietin-like protein 4 | RE-CON, post 30 min | -0.304 [-0.549, -0.060] | 0.015 | 0.164 | 0.152 | No source-FDR hit |
| IL-6 | RE-CON, post 3.5 h | -0.878 [-1.570, -0.186] | 0.0132 | 0.177 | 0.152 | No source-FDR hit |
| SPARC | RE-CON, post 30 min | -0.820 [-1.500, -0.139] | 0.0186 | 0.184 | 0.165 | No source-FDR hit |
| Follistatin | RE-CON, post 24 h | -0.474 [-0.795, -0.154] | 0.00399 | 0.22 | 0.054 | No source-FDR hit |
| Angiopoietin 1 | RE-CON, post 30 min | -0.899 [-1.749, -0.049] | 0.0384 | 0.263 | 0.285 | No source-FDR hit |
| HSP72 | RE-CON, post 10 min | +0.533 [-0.207, 1.273] | 0.157 | 0.552 | 0.676 | No source-FDR hit |
| IL-7 | EE-CON, during 20 min | +0.355 [-0.234, 0.944] | 0.235 | 0.555 | 0.774 | No source-FDR hit |
| SDC4 | EE-CON, during 40 min | +0.358 [-0.366, 1.082] | 0.33 | 0.588 | 0.808 | No source-FDR hit |
| Myostatin (GDF8) | EE-CON, during 40 min | +0.153 [-0.183, 0.489] | 0.369 | 0.632 | 0.824 | No source-FDR hit |

The lactate row above is the clinical assay. The research metabolomics assay also passes source FDR, with its smallest adjusted p = 2.83e-64. Agreement of two assays in this cohort is not independent-cohort replication.

## Every source-FDR-significant plasma result

| Candidate | Platform | Condition | Effect | Source adjusted p |
|---|---|---|---|---|
| Fractalkine | prot-ol | EE-CON, during 20 min | +0.5681 | 0.000364 |
| Fractalkine | prot-ol | EE-CON, during 40 min | +0.5605 | 0.000227 |
| Lactate | metab-t-clinical | EE-CON, during 20 min | +1.8881 | 1.63e-45 |
| Lactate | metab-t-clinical | EE-CON, during 40 min | +1.9068 | 5.7e-46 |
| Lactate | metab-t-clinical | EE-CON, post 10 min | +1.4916 | 1e-29 |
| Lactate | metab-t-clinical | EE-CON, post 30 min | +0.8550 | 4.04e-11 |
| Lactate | metab-t-clinical | RE-CON, post 10 min | +2.6865 | 3.62e-81 |
| Lactate | metab-t-clinical | RE-CON, post 30 min | +1.9971 | 2.95e-52 |
| Lactate | metab-u-ionpneg | EE-CON, during 20 min | +1.1458 | 1.56e-39 |
| Lactate | metab-u-ionpneg | EE-CON, during 40 min | +1.1024 | 2.32e-36 |
| Lactate | metab-u-ionpneg | EE-CON, post 10 min | +0.8586 | 3.43e-23 |
| Lactate | metab-u-ionpneg | EE-CON, post 30 min | +0.5226 | 1.71e-09 |
| Lactate | metab-u-ionpneg | RE-CON, post 10 min | +1.5265 | 2.83e-64 |
| Lactate | metab-u-ionpneg | RE-CON, post 30 min | +1.1953 | 2.49e-44 |

## Exercise-mode comparisons

Only direct EE-RE contrasts support a difference between exercise protocols. Significance in one exercise-versus-control arm alone is insufficient. A negative EE-RE effect means the endurance change is lower than the resistance change. Timing refers to time after the end of each protocol, not equal elapsed time after starting exercise.

| Candidate | Platform | Time | EE-RE effect | Source adjusted p |
|---|---|---|---|---|
| Lactate | metab-u-ionpneg | post 10 min | -0.6679 | 5.01e-21 |
| Lactate | metab-u-ionpneg | post 30 min | -0.6728 | 2.24e-21 |
| Lactate | metab-t-clinical | post 10 min | -1.1949 | 3.42e-28 |
| Lactate | metab-t-clinical | post 30 min | -1.1421 | 2.33e-26 |

CX3CL1 has no significant direct post-exercise EE-RE contrast (minimum source adjusted p = 0.472). Its during-endurance response cannot be compared with during-resistance, which was not sampled.

## Tissue and blood-cell context

Shown below is the strongest source-FDR result per candidate/tissue/assay layer with at least one hit. The full JSON retains all significant and nonsignificant rows. Blood RNA is a cellular expression measurement, distinct from plasma protein. No new cross-tissue omnibus correction was applied. Phosphosite effects do not directly establish protein activity.

| Candidate | Tissue | Layer / matched gene | Condition | Effect | Source adjusted p |
|---|---|---|---|---|---|
| Angiopoietin 1 | blood | transcript-rna-seq / ANGPT1 | RE-CON, post 30 min | +0.734 | 0.00365 |
| Angiopoietin-like protein 4 | muscle | transcript-rna-seq / ANGPTL4 | EE-CON, post 3.5 h | +1.565 | 0.00135 |
| Fractalkine | adipose | transcript-rna-seq / CX3CL1 | EE-CON, post 45 min | +0.588 | 0.0273 |
| Fractalkine | muscle | transcript-rna-seq / CX3CL1 | EE-CON, post 15 min | +3.626 | 5.5e-53 |
| HSP72 | blood | transcript-rna-seq / HSPA1B | EE-CON, post 10 min | +0.694 | 3.57e-10 |
| HSP72 | muscle | prot-ph / HSPA1A | RE-CON, post 15 min | -0.173 | 0.041 |
| HSP72 | muscle | transcript-rna-seq / HSPA1A | RE-CON, post 3.5 h | +1.581 | 1.04e-12 |
| IL-8 | blood | transcript-rna-seq / CXCL8 | RE-CON, post 30 min | +0.447 | 0.00486 |
| IL-15 | blood | transcript-rna-seq / IL15 | RE-CON, post 10 min | -0.381 | 0.044 |
| IL-15 | muscle | transcript-rna-seq / IL15 | RE-CON, post 3.5 h | -0.391 | 7.66e-05 |
| Lactate | muscle | metab / Lactic acid | EE-CON, post 24 h | -0.470 | 0.0478 |
| Myostatin (GDF8) | muscle | transcript-rna-seq / MSTN | RE-CON, post 3.5 h | -1.449 | 1.43e-10 |
| SPARC | blood | transcript-rna-seq / SPARC | RE-CON, post 30 min | +0.318 | 1.52e-05 |
| SDC4 | muscle | transcript-rna-seq / SDC4 | RE-CON, post 3.5 h | +2.228 | 3.03e-32 |
| VEGF | blood | transcript-rna-seq / VEGFA | EE-CON, post 3.5 h | +0.370 | 0.0101 |
| VEGF | muscle | transcript-rna-seq / VEGFA | EE-CON, post 3.5 h | +1.579 | 6.64e-27 |

HSP72 plasma is mapped to HSPA1A. Tissue context also includes HSPA1B as an explicitly labeled related gene; its RNA response is not a measurement of plasma HSPA1A.

## Multiple testing and uncertainty

The primary threshold is the original `adj_p_value < 0.05`, which uses BH within tissue, assay, platform, and contrast. It does not control an additional omnibus family spanning every timepoint and assay in this report. The secondary candidate-family calculation applies BH to all 230 raw plasma exercise-versus-control p-values. This is the same family used in the earlier Table 3 screen: the 12 unmatched/ambiguous candidates supplied no included tests. It is exploratory because these data have already been examined. FGF21 and IL-15 should be described as secondary findings requiring confirmation, not upgraded to primary discoveries.

A 95% CI excluding zero or a raw p below 0.05 can coexist with a nonsignificant adjusted p because the displayed CI is pointwise, not multiplicity-adjusted. Selecting a candidate's smallest p-value can exaggerate its apparent strength. All assay IDs, timepoints, and nonsignificant results are retained in `results.json` to make that selection visible. A nonsignificant result is not an equivalence test. Neither plasma abundance nor tissue RNA establishes secretion, disease benefit, or causal mediation.

## Interpretation for follow-up

1. **CX3CL1** is the clearest protein candidate for a disease-evidence overlay: circulating endurance-associated increases and strong muscle RNA responses, with source attribution and disease direction still unresolved.
2. **Lactate** provides a strong acute metabolic response and a directly supported difference between exercise protocols. Its known abundance response does not by itself establish a novel disease mechanism.
3. **FGF21 and IL-15** merit targeted validation of the observed post-resistance decreases, retaining their exploratory status and specific sampling times.
4. **MSTN, ANGPTL4, VEGFA, and SDC4** are candidates for tissue-response hypotheses; their RNA responses should not be presented as significant circulating exerkine responses.

## Reproducibility

The source-object verification independently checks extracted effects, confidence intervals, raw p-values, and adjusted p-values against the pinned R data objects. The Python analysis checks unique keys, valid p-value ranges, confidence-interval bounds, candidate counts, contrast coverage, and reproduces the 230-test BH correction. No restricted participant data were used.

Sources: [public MoTrPAC package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77), [contrast and FDR definitions](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/R/load_differential_analysis.R), [portal collection](https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0).
