# MoTrPAC exercise–disease hypothesis shortlist

Exploratory work, 25 September 2026 (Los Angeles). All quantitative exercise results below were extracted from public human MoTrPAC aggregate tables, not inferred from Exerkine Atlas. This is a screen of existing results, not a newly fitted analysis or a claim of discovery.

## Recommendation

Choose **type 2 diabetes (T2D)** as the single disease, with **CX3CL1/fractalkine as the lead signaling candidate** and **branched-chain amino acids (BCAAs) as a contrasting metabolic-biomarker module**. The core question is: *Which diabetes-related exercise responses plausibly represent inter-organ signaling, and which reflect short-term substrate handling?*

This fits the [hackathon disease track](https://stanford.bioinformatics-center.org/). The suggested deliverable is a small evidence explorer with observed effects, uncertainty, disease links, competing explanations, and one falsifiable hypothesis per candidate. It should label evidence as measured in MoTrPAC, reported in an external experiment, or proposed interpretation.

## Available inputs actually checked

- [Public human analysis package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis), pinned to commit `535b4044e7417413de471104c619120337602b77`.
- Plasma Olink differential results; muscle and adipose RNA differential results; adipose phosphoproteomics; plasma and muscle metabolomics; official feature-to-gene mapping.
- Both exercise-vs-control difference-in-change contrasts and direct endurance-vs-resistance contrasts.
- [Exerkine Atlas](https://exerkineatlas.org/explore.html), used for candidate nomination and references. Manual gene aliases are exported. A gene/precursor measurement is not automatically a measurement of its mature circulating product.
- [Exerkines in health, resilience and disease](https://www.nature.com/articles/s41574-022-00641-2), used as a conceptual framework; disease mechanism claims below are anchored to primary studies.

The cohort comprises healthy sedentary adults undergoing acute exercise. These aggregate inputs do not provide patient disease outcomes, participant-level cross-tissue covariance, or exercise-response genetic associations. The consortium distinguishes its public analysis package from its access-gated individual-level data package in its [repository documentation](https://github.com/MoTrPAC/motrpac-human-presuspension-acute).

## 1. Fractalkine signaling and type 2 diabetes

**Observed.** Plasma CX3CL1 increases during endurance exercise at 40 minutes (logFC +0.5605; adjusted p=0.0002274). Muscle CX3CL1 RNA increases at 15 minutes after endurance (+3.6261; adjusted p=5.50e-53) and resistance (+2.5854; adjusted p=1.55e-33). Adipose RNA increases at 45 minutes after endurance (+0.5878; adjusted p=0.0273).

**Disease anchor.** Experiments implicate CX3CL1/CX3CR1 in islet biology. A [human-islet study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4209359/) found reduced beta-cell apoptosis and reduced glucagon secretion, but no increase in human-islet glucose-stimulated insulin secretion. Protective effects against TNF-alpha on secretion were studied in rat beta cells. This distinction should remain visible. A separate [Cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3717389/) supplies additional mechanistic evidence.

**Hypothesis.** The exercise-associated CX3CL1 pulse may support islet resilience or glucagon regulation rather than simply increase insulin secretion.

**Analysis possible now.** Assemble the measured plasma and tissue-RNA time courses with the externally documented CX3CR1/islet relationship. Test whether the pattern is compatible with new transcription driving the initial circulating pulse. It is already temporally challenging: the plasma signal occurs during exercise, before the first post-exercise muscle biopsy. Pre-existing protein release, shedding, clearance, or an unmeasured source remain alternatives.

**What would weaken the claim.** A proposed *new muscle transcription causes the immediate plasma rise* mechanism is not supported by temporal precedence here. The direct post-exercise EE–RE plasma comparisons are all nonsignificant, and there are no during-RE samples; do not label this an endurance-specific exerkine.

**Later validation, not a hackathon dependency.** A short, physiologically calibrated CX3CL1 exposure in human islets, with CX3CR1 blockade, could test survival and glucagon responses. MoTrPAC alone does not demonstrate islet benefit or identify the source of circulating CX3CL1.

## 2. BCAA reduction: diabetes-related biomarker or increased oxidation?

**Observed.** At 30 minutes after resistance exercise, plasma leucine decreases (logFC -0.2024; adjusted p=9.67e-19), as do isoleucine and valine. The direct EE–RE leucine contrast is significant (+0.1434; adjusted p=1.07e-13), supporting a larger decrease after RE under these protocols. Muscle leucine decreases at 3.5 hours after RE (-0.3239; adjusted p=0.00277). Yet muscle BCKDHB and PPM1K RNA decrease at 24 hours after RE; both have strongly significant EE–RE contrasts.

**Disease anchor.** BCAAs are connected with T2D and insulin resistance, but causal interpretation is unsettled. A [bidirectional genetic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10827349/) cautions against simply treating high BCAA concentrations as upstream causes of diabetes.

**Hypothesis.** The rapid BCAA decrease after RE reflects substrate utilization or redistribution and is not adequately explained by transcriptional induction of BCAA oxidation enzymes.

**Analysis possible now.** Compare plasma and muscle leucine/isoleucine/valine, measured ketoleucine, and the available BCAA-pathway transcripts. Preserve the direct modality contrasts and the time lag. Attach diabetes biomarker and genetic evidence as separate annotations, rather than treating lower BCAA as proof of benefit.

**What would weaken the claim.** A coherent, correctly timed rise in available oxidation-related molecular evidence would favor an oxidation mechanism; the observed later transcript decreases already weaken a simple transcriptional-upregulation explanation. Concentrations and RNA cannot establish metabolic flux or distinguish uptake, protein synthesis, oxidation, and clearance by themselves.

**Later validation.** Stable-isotope BCAA tracing and protein-synthesis measurements could distinguish these mechanisms.

## 3. Kynurenine handling and depression biology — strongest alternative disease

**Observed.** Ten minutes after endurance exercise, plasma kynurenic acid increases (+0.4478; adjusted p=9.74e-8) while kynurenine decreases (-0.1464; adjusted p=0.0291). Direct EE–RE contrasts are significant for both molecules at that time (adjusted p=0.00341 and 1.24e-5 respectively). At 3.5 hours, the kynurenic-acid change is greater after RE (EE–RE adjusted p=0.0369). Muscle PPARGC1A RNA increases at 3.5 hours in both groups, while AADAT RNA decreases at that time in both groups.

**Disease anchor.** [Agudelo et al.](https://pubmed.ncbi.nlm.nih.gov/25259918/) linked muscle PGC-1alpha1, kynurenine conversion, and resistance to stress-induced depression in mouse experiments. That mechanism supplies a testable connection, not evidence that these MoTrPAC participants experienced an antidepressant effect.

**Hypothesis.** Acute human kynurenine handling differs by exercise modality and may use existing enzyme activity or extra-muscular processes before transcriptional adaptation. The simple explanation that acute PPARGC1A induction immediately increases AADAT transcription is inconsistent with the measured timing and direction.

**Analysis possible now.** Build separate plasma kynurenine and kynurenic-acid trajectories and align the measured muscle pathway transcripts. Use the available direct EE–RE tests. Do not construct a tested participant-level kynurenic-acid/kynurenine ratio from marginal summary statistics: its uncertainty requires covariance that is not available here.

**What would weaken the disease interpretation.** Failure to connect these human peripheral changes to a disease-relevant mechanism, or persistence of alternative clearance/source explanations, would leave this as a metabolic response. Brain measurements, depression outcomes, enzyme flux, and PGC-1alpha1-specific isoform activity were not established by this screen.

## 4. Neprilysin and heart failure — a useful counterexample

**Observed.** Plasma MME/neprilysin increases during endurance exercise at 40 minutes (+0.3178; adjusted p=0.0262).

**Disease anchor.** The [PARADIGM-HF trial](https://www.nejm.org/doi/full/10.1056/NEJMoa1409077) demonstrated benefit of combined angiotensin-receptor/neprilysin inhibition relative to enalapril. That trial does not isolate neprilysin inhibition from the combination therapy.

**Hypothesis.** Exercise-associated plasma MME abundance can mark transient protease release or regulation without implying that increasing tissue neprilysin activity is beneficial in heart failure.

**Analysis possible now.** Plot measured MME dynamics and annotate the disease intervention direction. Use this as a transparent example of why exercise direction, disease association, drug action, protein abundance, and enzyme activity must be separate fields in an evidence map.

**Limit.** This is a strong interpretive demonstration but a weaker standalone mechanistic project: the available protein assay does not resolve source, membrane versus soluble pools, or catalytic activity.

## Analysis rules

1. Use exercise-versus-time-matched-control contrasts, not uncontrolled pre/post changes.
2. Claim differences between modalities only from direct EE–RE contrasts. A significant result in one arm and a nonsignificant result in another is insufficient.
3. Keep tissue-specific time labels: early blood 30 min, muscle 15 min, adipose 45 min; middle blood/muscle 3.5 h and adipose 4 h.
4. Source adjusted p-values are BH corrections within tissue/assay/platform/contrast, not a new correction across this exploratory screen. Preserve uncertainty and preregister the smaller follow-up question.
5. Distinguish absent-from-analyzed-table, measured-but-nonsignificant, and detected response. For example, CLU and BDNF are absent from the analyzed plasma Olink table; IL6 is present but no tested exercise-control result passes source adjusted p<0.05 in this screen. Neither observation disproves their biology.
6. Tissue RNA plus plasma protein does not prove secretion or tissue origin. Phosphosite abundance is not automatically enzyme activity. A pathway concentration change is not flux.
7. Do not estimate individual responder classes, mediation, cross-tissue participant correlations, or clinical prediction from aggregate group statistics.
8. If adding genetics, use established variant-to-gene/colocalization evidence with tissue and direction recorded. A nearest-gene annotation or generic target–disease score does not establish causation. No new genetics analysis was performed here.

## Files and reproduction

- `candidate_evidence.csv`: 1,577 selected feature/contrast rows, including nonsignificant results, source model effects, 95% confidence intervals, raw p-values, and source adjusted p-values.
- `atlas_plasma_coverage.csv` and `atlas_protein_gene_mapping.csv`: explicit candidate mapping and coverage checks.
- `exercise_candidate_trajectories.png` / `.svg`: selected plasma results; equal-spaced categorical time points, not a fitted kinetic model.
- `reproduce_screen.R`: downloads the pinned public data and exports the evidence table with base R.
- `plot_results.py`: creates the figure with pandas, NumPy, and Matplotlib.
- `analysis_provenance.txt` and `input_sha256.json`: source version and input checksums.

Run:

```sh
Rscript reproduce_screen.R .
python3 plot_results.py
```

The project proposals are new synthesis of existing observations. Their clinical relevance remains a hypothesis; the plotted effects themselves are existing consortium results, not independent replication.
