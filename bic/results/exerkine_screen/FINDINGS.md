# Additional candidate exerkines in human acute-exercise data

Analysis: 26 September 2026. Open the executed [notebook](../../exerkine_candidate_screen.ipynb) for the code, full tables, and figures.

**Three candidates deserve follow-up: CD300LG for a diabetes-focused project, ANGPT2 for vascular remodeling, and WARS1 as an exploratory immune-signaling candidate.** These are nominations supported by exercise responses and biological evidence, not validated mediators of exercise benefit.

## What we screened

We used the available **human Acute Exercise in Human Sedentary Adults c2.0** differential-result objects, from package 2.0.8 pinned to commit `535b4044e7417413de471104c619120337602b77`. These are single-bout endurance (EE) and resistance (RE) responses, compared with resting controls (CON), not training adaptations. Raw p-values, effects, and confidence intervals came from MoTrPAC; we did not fit new participant-level models. [Pinned source package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77).

The primary nomination rule was a positive plasma effect with **source BH-adjusted p < 0.05**, plus a positive muscle/adipose RNA or total-protein effect passing the same threshold for the same gene and exercise mode. All available sampling times were considered and retained. Here “increase” means a greater baseline-to-timepoint change than in resting controls, not necessarily an absolute increase within the exercise group.

- **1,417 Olink plasma features** were tested in each comparison. Three features could not nominate a uniquely mapped gene, but remained in the BH denominator.
- **142 genes** had a qualifying positive plasma response.
- **15 genes** also had a qualifying positive muscle/adipose response after the same exercise mode.
- **13 of those 15** are additional to the previously analyzed CCN1 and CX3CL1.
- We prioritized **three** of those 13 after reviewing extracellular signaling evidence and research usefulness. The other ten remain visible and require further biological review.

BH was independently reproduced over complete assay families before gene selection: **58 families, 713,660 feature–contrast tests**, maximum absolute difference from the source adjusted values **4.44 × 10⁻¹⁵**. The previous 10 primary plasma results for each of CCN1 and CX3CL1 were reproduced. These checks establish consistency, not independent biological replication.

## Priority candidates

Each row identifies a specific matching comparison. “q” means the original assay-wide BH-adjusted p-value. Effects retain MoTrPAC's native `logFC` scale and are not concentrations or fold changes asserted here.

| Candidate | Exercise comparison | Positive plasma result | Positive muscle RNA result | Reason to prioritize |
|---|---|---|---|---|
| **CD300LG** | EE-CON | During 40 min: effect **0.696**, q **0.00176** | Post 3.5 h: effect **0.439**, q **0.00380** | Strongest external human diabetes rationale; circulating form needs resolution |
| **ANGPT2** | EE-CON | During 40 min: effect **0.401**, q **0.000206** | Post 3.5 h: effect **0.717**, q **0.000985** | Secreted vascular ligand; remodeling hypothesis |
| **WARS1** | RE-CON | Post 10 min: effect **0.567**, q **0.0312** | Post 3.5 h: effect **0.329**, q **3.63 × 10⁻⁶** | Nonclassical secretion and immune-signaling evidence outside exercise |

These matched-mode findings do not themselves establish differences between EE and RE. None of the three has a direct plasma EE-RE comparison passing source BH < 0.05 at the available post-exercise timepoints. The notebook displays every direct test.

**CD300LG.** A 2024 human training study linked circulating CD300LG with clamp-measured insulin sensitivity; population associations and Mendelian-randomization analyses supplied further glucose-related evidence. Our acute result motivates a hypothesis connecting an early circulating response to longer-term vascular/metabolic adaptation. The protein is membrane-associated, so the circulating species and release route require verification. Existing associations do not establish treatment benefit. [Lee-Ødegård et al., 2024](https://pubmed.ncbi.nlm.nih.gov/39190027/).

**ANGPT2.** It is a secreted vascular ligand already highlighted by MoTrPAC. Adipose-specific mouse experiments connected ANGPT2-driven vascularization with metabolic improvement. Our muscle RNA/plasma observations do not identify the releasing cell or establish that raising systemic ANGPT2 benefits diabetes. Vascular bed, dose, and timing matter to the hypothesis. [An et al., 2017](https://elifesciences.org/articles/24071), [MoTrPAC landscape preprint](https://doi.org/10.64898/2026.02.27.702183).

**WARS1.** This cytoplasmic tRNA synthetase can also be released directly or in vesicles and act in innate immune signaling, as shown experimentally outside exercise. That supplies a reason to consider its plasma response as more than an intracellular abundance marker. The exercise-specific release mechanism and biological activity remain untested. It is not currently our strongest diabetes candidate. [Nguyen et al., 2023](https://pubmed.ncbi.nlm.nih.gov/36640342/).

## What is new relative to the published candidate list?

**Update, 27 September 2026:** the comparison below refers to the original preprint's Table S8. WARS1 is explicitly included in the current MoTrPAC Supplementary Table 7 and secretome Figure 7. The historical Table S8 absence does not establish a new nomination beyond the current consortium analysis. Earlier exercise-related literature also exists. See the [WARS1 literature review and version audit](../wars1/novelty_review/REVIEW.md).

We retrieved the original MoTrPAC landscape preprint's **Table S8**, verified its MD5 against the article's XML, and read the `Exerkine_Candidate_List` sheet: **582 rows representing 168 genes**. **Fourteen of our 15 genes are already listed; WARS1 is absent**, including a check for its alias WARS. CD300LG was also independently proposed as an exerkine in 2024. [MoTrPAC Table S8 source article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/).

Thus WARS1 extended this older checked list under our screen, but is already present in the current consortium list. This does **not** establish literature-wide novelty. Our source release and rule differ from the original preprint's broader screen: we require positive same-mode RNA/total-protein support at 0.05 and do not treat phosphosite responses as protein-production evidence. Absence from Exerkine Atlas is also not evidence of novelty.

The full primary overlap is **NADK, SERPINB8, ANGPT2, CX3CL1, CCN1, PXN, CD300LG, ACVRL1, EGLN1, SCLY, FGR, CEACAM8, WARS1, TMSB10, and CAMKK1**. Small p-values alone do not establish active secretion. Only ANGPT2, CX3CL1, and CCN1 have explicit secreted-location flags in our cached UniProt annotations. WARS1 and CD300LG illustrate why published secretion/signaling evidence must supplement that annotation rule. Some other genes have extracellular COMPARTMENTS evidence in Table S8; they have not been ruled out.

## Temporal evidence and statistical limits

**None of the 15 primary genes has a qualifying sampled tissue increase at or before a qualifying plasma increase.** For our three priorities, muscle RNA increases at 3.5 hours after an early plasma signal. Later transcription might replenish released protein or be a parallel response; it cannot establish that new transcription caused the earlier pulse.

Existing protein stores, other tissues, vascular/immune cells within biopsies, shedding, extracellular vesicles, clearance, plasma-volume shifts, and cell injury are alternative explanations. Affinity-assay measurements also need orthogonal verification of the circulating molecular species. These summary tables cannot provide within-person plasma–tissue correlations or establish causal effects on diabetes.

The source BH procedure addresses multiple tests within each assay/contrast family. Our post hoc intersection across times and assays is **not a combined hypothesis test** and does not establish a 5% FDR for the final candidate list. As a sensitivity check, BH across all **14,170 primary plasma tests together** retains the highlighted results for ANGPT2 (**0.000832**), CD300LG (**0.00580**), and WARS1 (**0.0408**). This different correction also does not validate the cross-tissue causal hypothesis.

## Useful secondary leads

| Candidate | Evidence | Why separate from the primary tissue-supported list? |
|---|---|---|
| **FLT3LG** | Multiple plasma hits; best q **9.67 × 10⁻⁸** | No positive same-mode muscle/adipose RNA or total-protein hit at 0.05 |
| **ANGPTL7** | Endurance plasma: during 40 min q **1.12 × 10⁻⁵** | No matching positive tissue hit at 0.05 |
| **HGF** | Plasma q **0.00808**; matching muscle RNA q **0.0514** | Tissue misses the primary cutoff |
| **STC2** | Plasma q **0.0602**; matching muscle RNA q **2.33 × 10⁻¹³** | Plasma misses the primary cutoff |
| **TIMP3** | Plasma q **0.0628**; matching muscle RNA q **0.0226** | Plasma misses the primary cutoff |

FLT3LG has prior human acute-exercise evidence. ANGPTL7 has prior evidence for lower levels after a training program in obesity; the acute rise here concerns a different study and timescale. Those are follow-up hypotheses, not necessarily conflicting findings. [FLT3LG study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9999360/), [ANGPTL7 study](https://doi.org/10.1371/journal.pone.0173024).

The separate q < 0.10 sensitivity screen yields **247 plasma-up genes and 38 same-mode tissue-supported genes**. It is explicitly exploratory; it does not change the primary threshold.

## Proposed tests

| Candidate | Hypothesis to test | Result that would weaken it |
|---|---|---|
| CD300LG | An exercise-responsive circulating CD300LG species participates in vascular or glucose regulation. First resolve full-length versus cleaved/vesicular species with an independent assay; then test a defined species in endothelial/glucose-handling experiments. | Failure to reproduce the plasma signal, or no functional effect after verified perturbation at relevant exposure. |
| ANGPT2 | The acute plasma response modifies endothelial TIE2 signaling or barrier behavior. Compare pre/post-exercise plasma with ANGPT2 depletion and controlled add-back. | Endothelial response is unchanged by effective ANGPT2 depletion/add-back. |
| WARS1 | The RE-associated circulating WARS1 species conveys an immune signal. Separate soluble/vesicular fractions, verify protein identity and injury markers, and test depletion/add-back in an innate-immune reporter assay. | Signal is not reproducible, is explained by nonspecific leakage, or has no WARS1-dependent extracellular activity. |

These are proposed experiments, not completed findings. The available data support prioritization and exact response hypotheses; they do not demonstrate secretion, target-organ activity, or disease protection.
