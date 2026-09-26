# Systematic MoTrPAC screen of Table 3

Source: Chow et al., *Exerkines in health, resilience and disease* (2022), Table 3, printed pages 279–280 (PDF pages 7–8). The supplied PDF was extracted and both complete table pages visually checked. All 28 entries were retained, including absent and unresolved matches.

This is a systematic candidate lookup in the public human MoTrPAC acute-exercise results, not a systematic review of all literature or independent experimental replication.

## Findings

- 16 of 28 entries have an identity-resolved match in the analyzed circulating protein/metabolite tables.
- 2 entries show at least one exercise-versus-control change at the original MoTrPAC adjusted p<0.05: Fractalkine, Lactate.
- 14 entries are tested in plasma but have no source-FDR hit; this is not evidence of equivalence or absence of biological response.
- 11 entries have no match in the analyzed plasma tables; that can reflect assay coverage or analysis/QC eligibility.
- 1 entry (BAIBA) is retained as an unresolved chemical annotation.

## Every Table 3 entry

Tissue columns summarize any source-adjusted p<0.05 across EE-CON or RE-CON: up, down, or both. “NS” means tested with no detected change; “—” means no matching analyzed feature. Tissue results are context and do not establish secretion.

| Table 3 entry | Plasma status | Muscle RNA | Adipose RNA | Muscle protein | Adipose protein |
|---|---|---|---|---|---|
| 12,13-diHOME | Absent from analyzed plasma tables | — | — | — | — |
| Adiponectin | Absent from analyzed plasma tables | NS | NS | NS | NS |
| Angiopoietin 1 | Tested; no source-FDR hit | NS | NS | — | — |
| Angiopoietin-like protein 4 | Tested; no source-FDR hit | up + down | NS | — | NS |
| Apelin | Absent from analyzed plasma tables | up | NS | — | — |
| BAIBA | Chemical identity unresolved | — | — | — | — |
| Catecholamines | Absent from analyzed plasma tables | — | — | — | — |
| Fetuin-A | Absent from analyzed plasma tables | — | — | NS | NS |
| Fractalkine | Detected response | up | up | — | — |
| FGF21 | Tested; no source-FDR hit | — | — | — | — |
| Follistatin | Tested; no source-FDR hit | NS | NS | — | NS |
| GDF15 | Tested; no source-FDR hit | — | NS | — | — |
| HSP72 | Tested; no source-FDR hit | up | NS | NS | NS |
| IL-6 | Tested; no source-FDR hit | — | NS | — | — |
| IL-7 | Tested; no source-FDR hit | — | NS | — | — |
| IL-8 | Tested; no source-FDR hit | — | NS | — | — |
| IL-15 | Tested; no source-FDR hit | down | NS | — | — |
| Irisin | Absent from analyzed plasma tables | NS | NS | — | — |
| Lactate | Detected response | — | — | — | — |
| METRNL | Absent from analyzed plasma tables | up | NS | — | — |
| Myonectin (CTRP15) | Absent from analyzed plasma tables | up | NS | — | — |
| Musclin (osteocrin) | Absent from analyzed plasma tables | — | — | — | — |
| Myostatin (GDF8) | Tested; no source-FDR hit | down | NS | — | — |
| SPARC | Tested; no source-FDR hit | NS | NS | NS | NS |
| SDC4 | Tested; no source-FDR hit | up | NS | NS | — |
| TGF-beta1 | Absent from analyzed plasma tables | up | NS | — | — |
| TGF-beta2 | Absent from analyzed plasma tables | up | NS | — | — |
| VEGF | Tested; no source-FDR hit | up | NS | NS | — |

## Circulating responses passing the original MoTrPAC FDR

These are existing consortium model estimates for exercise-related changes relative to time-matched non-exercise controls. Effects are reported in each source assay’s model scale; they are not assumed to be directly comparable between assay platforms.

| Candidate | Assay/platform | Contrast | Time | Effect | Source adjusted p | Table-3-family BH |
|---|---|---|---|---:|---:|---:|
| Fractalkine | prot-ol | EE-CON | during_20_min | 0.5681 | 0.000364 | 2.27e-05 |
| Fractalkine | prot-ol | EE-CON | during_40_min | 0.5605 | 0.000227 | 4.75e-05 |
| Lactate | metab-u-ionpneg | EE-CON | during_20_min | 1.146 | 1.56e-39 | 1.68e-39 |
| Lactate | metab-u-ionpneg | EE-CON | during_40_min | 1.102 | 2.32e-36 | 2.19e-36 |
| Lactate | metab-u-ionpneg | EE-CON | post_10_min | 0.8586 | 3.43e-23 | 3.88e-23 |
| Lactate | metab-u-ionpneg | EE-CON | post_15_30_45_min | 0.5226 | 1.71e-09 | 1.62e-09 |
| Lactate | metab-u-ionpneg | RE-CON | post_10_min | 1.527 | 2.83e-64 | 1.6e-63 |
| Lactate | metab-u-ionpneg | RE-CON | post_15_30_45_min | 1.195 | 2.49e-44 | 3.13e-44 |
| Lactate | metab-t-clinical | EE-CON | during_20_min | 1.888 | 1.63e-45 | 8.33e-45 |
| Lactate | metab-t-clinical | EE-CON | during_40_min | 1.907 | 5.7e-46 | 3.64e-45 |
| Lactate | metab-t-clinical | EE-CON | post_10_min | 1.492 | 1e-29 | 2.84e-29 |
| Lactate | metab-t-clinical | EE-CON | post_15_30_45_min | 0.855 | 4.04e-11 | 9.37e-11 |
| Lactate | metab-t-clinical | RE-CON | post_10_min | 2.686 | 3.62e-81 | 9.25e-80 |
| Lactate | metab-t-clinical | RE-CON | post_15_30_45_min | 1.997 | 2.95e-52 | 2.51e-51 |

## Multiple-testing definitions

The original MoTrPAC `adj_p_value` adjusts within tissue/assay/platform/contrast. It is preserved and used for the main coverage figure and response classifications.

An additional exploratory BH correction covers all 230 identity-resolved circulating candidate-by-assay-by-time-by-mode tests, using their raw p-values. It includes negative results, all measured Olink assay IDs, and both research and clinical lactate assays. BAIBA is excluded pending chemical-identity resolution. The smaller candidate family has a different null family from the source correction; passing it does not mean passing the source correction.

Candidates with hits only under this additional Table-3-family correction: FGF21, IL-15. These should be labeled exploratory and secondary.

The tissue-context screen does not receive a new cross-tissue/cross-time FDR correction here. Its source-FDR hits should not be treated as a newly controlled omnibus discovery list.

## Identity and interpretation checks

- **BAIBA:** `METABOLOMICS_CVS` includes a raw `beta-Aminoisobutyric-acid` assay label mapped to `Aminoisobutyric acid`. The public gene/feature map associates that name with KEGG C03665, which identifies 2-/alpha-aminoisobutyric acid. The selected blood assay uses the generic name; the muscle raw label uses beta. Those inconsistent annotations require resolution before making beta-BAIBA claims. [KEGG C03665](https://www.genome.jp/dbget-bin/www_bget?cpd%3AC03665=).
- **12,13-diHOME:** no exact analyzed match was found. 9,10-diHOME, which appears in annotation metadata, is a different regioisomer and was excluded.
- **Irisin/FNDC5:** FNDC5 transcript or tissue protein is a precursor measurement, not proof of mature circulating irisin.
- **Myonectin:** mapped to ERFE, retaining aliases CTRP15/FAM132B. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/151176). Musclin maps to OSTN. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/344901).
- **HSP72:** HSPA1A/HSPA1B are grouped as candidate paralogs, but unique physical assay features are tested once per candidate/contrast.
- **VEGF:** operationally mapped to VEGFA; other VEGF-family genes are not silently substituted.
- **Review arrows:** the table footnote defines them as changes in plasma levels. The chronic-training column is retained as metadata but is not tested with this acute-bout dataset.
- **Species and effect:** H/A/C encode human/animal/cell evidence; A/P/E in the source effect column mean autocrine/paracrine/endocrine. A row listing H/A/C does not imply that every biological action was demonstrated in humans.
- **Timing:** `post_15_30_45_min` is 15 min for muscle, 30 min for blood, 45 min for adipose. `post_3.5_4_hr` is 3.5 h for muscle/blood and 4 h for adipose. During-exercise RE sampling is unavailable.
- **Nonsignificance:** not a formal equivalence test, and it does not contradict a heterogeneous review claim across different exercise doses, populations, species, sampling times, or assay sensitivity.
- **Causality:** RNA, tissue abundance, or phosphosite changes do not establish release into blood, source tissue, enzyme activity, target-organ response, or clinical benefit.

## Scope and reproducibility

Public human pre-suspension acute-exercise aggregate release, pinned to commit `535b4044e7417413de471104c619120337602b77`. Thirteen differential-result objects were searched: RNA in blood, muscle and adipose; plasma Olink protein; muscle/adipose protein and phosphoprotein; metabolomics in plasma, muscle and adipose; and the available blood clinical chemistry/protein tables. Clinical protein tables added no Table 3 matches. Epigenomics is outside this abundance/response screen.

Full output contains 1,171 candidate/assay/contrast rows, including the direct EE-RE comparisons. The primary summaries use only exercise-with-controls contrasts. No individual-level data, source assignment, disease-cohort comparison, or new differential model was inferred.

Files: `table3_candidates.csv` (28-row source manifest), `coverage_summary.csv`, `all_candidate_results.csv`, `plasma_results.csv`, `assay_inventory.csv`, `coverage_matrix.png`, and `input_provenance.json`. Source scripts are one directory above: `table3_manifest.py`, `table3_screen.R`, and `table3_report.py`.

```sh
python3 table3_manifest.py
Rscript table3_screen.R .
python3 table3_report.py
```
