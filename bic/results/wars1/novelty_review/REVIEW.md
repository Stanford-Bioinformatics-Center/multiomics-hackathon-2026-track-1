# WARS1: literature review and MoTrPAC novelty audit

Reviewed 27 September 2026. This is a targeted scoping review, not an exhaustive systematic review of all exercise-proteomics supplements.

**Conclusion: remove the unqualified word “novel.” WARS1 is explicitly nominated in the current MoTrPAC exerkine table and appears in its secretome figure. Earlier exercise-related studies also exist. Its vascular function after exercise remains a hypothesis to test.**

The earlier project statement concerned one older supplementary table. Extending that observation to “MoTrPAC has not identified WARS1” was incorrect. This audit corrects that interpretation without changing the differential-analysis results.

## What changed in the MoTrPAC check

| Source actually inspected | Candidate-table result | Interpretation |
|---|---|---|
| Original landscape preprint, Table S8; cached file MD5 `beaf62ed40e4618b6655078fef49d310` | 582 rows, 168 unique genes; no WARS1/WARS match | Our earlier absence statement was true for this file |
| Current public figure repository, Supplementary Table 7, pinned at `ed8d37fd2c9805b47919ec7d38e0aae27aa22ccf` | 341 rows, 115 unique genes; **two WARS1 rows** | WARS1 is already a consortium-nominated candidate |
| Same repository, `Figure7_Secretome.pdf` | **WARS1 label present** in extracted figure text | It is also displayed in the current figure |

Sources: [original preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/), [pinned current Table 7](https://github.com/MoTrPAC/motrpac-human-presuspension-acute/blob/ed8d37fd2c9805b47919ec7d38e0aae27aa22ccf/figures/landscape/assembled/Landscape_Supplementary_Table_7.xlsx), [pinned current Figure 7](https://github.com/MoTrPAC/motrpac-human-presuspension-acute/blob/ed8d37fd2c9805b47919ec7d38e0aae27aa22ccf/figures/landscape/assembled/Figure7_Secretome.pdf).

The current WARS1 rows are muscle RNA responses at `EE_Post_3.5_4hr` and `RE_Post_3.5_4hr`. Both report tissue direction **Up**, plasma direction **Up**, extracellular score **4.371**, directional concordance **Yes**, and temporal concordance **No**. The extracted rows are saved in [motrpac_current_wars1_rows.csv](motrpac_current_wars1_rows.csv).

The table's `Modality = Both` is assigned from tissue responses across exercise arms. It does **not** mean that plasma WARS1 passed BH after both arms. The code matches tissue and plasma results by gene and uses source adjusted p < 0.05. [Candidate-generation code](https://github.com/MoTrPAC/motrpac-human-presuspension-acute/blob/ed8d37fd2c9805b47919ec7d38e0aae27aa22ccf/figures/landscape/figure_7/FIG7_ED9_helpers.R).

The current workbook path appears in a repository reorganization dated 22 September 2026. That establishes availability by that commit; it does **not** establish the first date on which MoTrPAC nominated WARS1. We did not reconstruct the complete history across renamed files. Repository artifacts are distinguished here from the original published preprint supplement.

## Earlier exercise evidence

| Study | Direct observation | What it establishes and what it leaves open |
|---|---|---|
| **Nurbekov et al., 2015**, athlete overtraining pilot | Examined tryptophanyl-tRNA synthase gene expression around training | Earlier exercise-related transcript work. Publisher abstract checked; assay-to-modern-gene mapping and full experimental details were not independently revalidated. It does not establish circulating secretion. [Primary article](https://genescells.ru/2313-1829/article/view/120500) |
| **Aas et al., 2023**, human muscle-cell experiment | Table 2 reports **WARS protein 2.01-fold higher, p = 0.009**, in the exosome fraction after electrical pulse stimulation | Exercise-model evidence for altered released vesicle cargo already exists. Cells came from six women with severe obesity and type 2 diabetes. The differential-protein test used paired t-tests; the reported p is **raw**, not a BH-adjusted q. This is not an in-vivo plasma or vascular-function experiment. [Full paper and Table 2](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2023.1143966/full) |
| **Abbasi et al., 2025**, CHEST conference abstract | WARS1 **decreased** in serum-derived extracellular vesicles during an acute exercise test before training in a 14-person long-COVID cohort | Prior human exercise association, with a different population, compartment and direction. Abstract-only evidence: no WARS1-specific numerical effect/p/q or clearly named multiple-testing correction is supplied. This is not replication of the MoTrPAC plasma increase. [Primary abstract](https://doi.org/10.1016/j.chest.2025.07.3648) |

The 2023 protein-name field includes T1/T2 aliases; that label alone does not establish isoform- or cleavage-specific detection. None of these exercise reports demonstrates a WARS1-dependent endothelial response.

## The vascular rationale also predates this project

- **T2-WARS:** experimental work linked this truncated form's anti-angiogenic activity to VE-cadherin. [Tzima et al., 2005](https://doi.org/10.1074/jbc.C400431200).
- **Mini-WARS:** experimental work showed inhibition of NRP1-associated VE-cadherin turnover and vascular permeability. [Gioelli et al., 2022](https://www.nature.com/articles/s41467-022-31904-1).
- **Extracellular release:** direct and vesicular WARS1 secretion and inflammatory activity have been studied outside exercise. [Nguyen et al., 2023](https://pubmed.ncbi.nlm.nih.gov/36640342/).

These papers support biological plausibility. They do not identify the molecular form detected by the exercise Olink assay or demonstrate vascular mediation of an exercise response. “Vascular” can refer to barrier regulation or vessel growth, and those outcomes should be tested separately.

## What our analysis contributes

We extracted and checked existing MoTrPAC summary statistics from human acute exercise, using analysis package 2.0.8 at commit `535b4044e7417413de471104c619120337602b77`. We did not generate a new cohort, refit participant-level models, or independently discover a signal absent from MoTrPAC's current candidate list. [Pinned source package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77).

The useful contribution is a focused, transparent assessment:

- Plasma WARS1 after resistance exercise, 10 minutes: control-adjusted source effect **+0.567**, original BH-adjusted p **0.0312**.
- No endurance-versus-control plasma WARS1 result passes the original BH threshold of 0.05.
- Muscle RNA increases at 3.5 hours after both exercise modes. The later RNA response does not identify the source of the earlier plasma signal.
- Direct endurance-versus-resistance comparisons do not establish a difference between modes; the post-10-minute adjusted p is **0.227**.
- Olink feature **OID21084**, mapped to **WARS1/P23381**, does not resolve full-length, mini-WARS, T2-WARS, or vesicular carriage in these samples.

See [primary plasma results](../wars1_plasma_primary.csv) and [detailed WARS1 findings](../FINDINGS.md). The significance checks remain valid; the novelty interpretation changes.

## A defensible research claim and next test

Suggested slide wording:

> **WARS1 is a MoTrPAC-nominated candidate exerkine with a plausible vascular role; whether its exercise-associated circulating forms alter endothelial function remains untested in the studies reviewed.**

A testable hypothesis is that a particular WARS1 species contributes to endothelial responses during recovery from exercise. Start by resolving soluble versus vesicular protein and full-length versus truncated forms in paired pre/post-exercise samples using an independent assay. Then test whether WARS1 depletion changes an endothelial response and whether form-specific add-back restores it. Measure permeability and angiogenesis separately. Include sample timing, plasma-volume change, cell injury and assay specificity as alternative explanations.

This is a proposed mechanistic contribution, not a verified claim that no one has previously proposed the idea. The available summary data support candidate prioritization; they cannot demonstrate this mechanism.

## Search scope, reproducibility and limitations

Searches used PubMed, primary publisher pages, full-text tables and MoTrPAC's public repositories. Terms included WARS1, WARS, TrpRS, tryptophanyl-tRNA synthetase/synthase, and tryptophan-tRNA ligase, combined with exercise, training, contraction, myokine, exerkine or electrical pulse stimulation. A further web query checked WRS with exerkine/WARS. Protein identity was checked against cytoplasmic WARS1 rather than mitochondrial WARS2.

The first broad PubMed query produced 5,022 results because of alias ambiguity, including “wars”; those records were **not all screened**. A date-limited title/abstract refinement returned 10 records, all screened and excluded as irrelevant to WARS1 exercise biology. A narrower exerkine/myokine query returned zero. These negative searches do not establish absence: the relevant 2023 result is embedded in a full-text protein table and the 2025 finding is a conference abstract.

Exact PubMed queries are in [search_log.json](search_log.json), with retrieved records and exclusion reasons in [pubmed_screening.csv](pubmed_screening.csv). Downloaded sources are in `reference_data/`; [source_manifest.json](source_manifest.json) records URLs, the pinned repository version and file checksums. [candidate_list_audit.json](candidate_list_audit.json) records workbook counts and matches.

The cached Frontiers HTML has one unrelated publisher Google Maps API key redacted. The source manifest records both the original download checksum and the sanitized file checksum; the scientific article content is unchanged.

No exhaustive scan of every published exercise-proteomics supplement, paywalled full text, preprint version or renamed repository artifact was completed. Positive prior evidence is sufficient to reject the broad novelty claim; it does not establish which study was first.
