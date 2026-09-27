# Fractalkine and WARS1 differential analysis

Four-slide presentation based on the available human MoTrPAC c2.0 acute-exercise results.

## Slide 1: Human MoTrPAC differential analysis

- Data: human MoTrPAC c2.0 acute exercise in sedentary adults, analysis package 2.0.8.
- Effect = (exercise at timepoint − exercise baseline) − (control at timepoint − control baseline).
- Null hypothesis: this difference in changes equals zero for the molecule and timepoint.
- MoTrPAC supplied mixed-model effects and raw p-values. We extracted candidates and independently reproduced BH.
- A BH-adjusted p < 0.05 defines a hit. BH accounts for all features in each assay/comparison, including 1,417 plasma Olink features. It targets a 5% expected false-discovery fraction among reported hits under its assumptions.

## Slide 2: Fractalkine rises during endurance exercise

| Endurance vs control | Effect | 95% CI | Raw p | BH-adjusted p |
| --- | --- | --- | --- | --- |
| During 20 min | +0.568 | [0.345, 0.791] | 1.29 × 10⁻⁶ | 0.000364 |
| During 40 min | +0.561 | [0.332, 0.789] | 2.89 × 10⁻⁶ | 0.000227 |

- CX3CL1 plasma protein, Olink feature OID20976. Effects use source logFC and confidence intervals are pointwise.
- 2 of 10 plasma exercise-versus-control tests pass BH < 0.05. Both also pass Holm correction across the ten CX3CL1 plasma tests.
- Muscle RNA rises 15 min after endurance exercise (+3.626, BH p = 5.50 × 10⁻⁵³) and resistance exercise (+2.585, BH p = 1.55 × 10⁻³³).
- No post-exercise plasma test passes BH < 0.05. RNA and protein effects are on different assay scales.

## Slide 3: WARS1 rises after resistance exercise

| Measurement | Exercise vs control | Time after exercise | Effect [95% CI] | BH-adjusted p |
| --- | --- | --- | --- | --- |
| Plasma protein | Resistance | 10 min | +0.567 [0.242, 0.892] | 0.0312 |
| Muscle RNA | Endurance | 3.5 h | +0.374 [0.245, 0.503] | 1.48 × 10⁻⁶ |
| Muscle RNA | Resistance | 3.5 h | +0.329 [0.206, 0.453] | 3.63 × 10⁻⁶ |

- Early plasma protein response and later muscle RNA responses. Plasma Olink feature OID21084.
- These are the only 3 BH hits among 52 primary WARS1 estimates. Muscle total protein does not show a significant exercise response.
- Direct endurance–resistance plasma comparison at 10 min: raw p = 0.0442, BH p = 0.227. Resistance specificity is not established.
- Effects use source logFC. CIs are pointwise. Muscle RNA at 3.5 h does not establish the origin of the earlier plasma response.

## Slide 4: Interpretation and testable next steps

- Fractalkine: reproduce the 20–40 min endurance plasma increase using an independent protein assay and matched resting controls.
- WARS1: confirm the early resistance plasma increase and determine which WARS1 molecular form the assay detects.
- Neither plasma analysis establishes exercise-mode specificity. During-resistance samples are unavailable and the direct post-exercise plasma comparisons do not pass BH.
- Later tissue RNA cannot identify the source of an earlier plasma signal. These exploratory candidates require replication. The current analysis does not test diabetes benefit.

## Sources

- [Pinned MoTrPAC analysis package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77)
- [Fractalkine notebook](../../fractalkine_differential_analysis.ipynb)
- [WARS1 notebook](../../wars1_differential_analysis.ipynb)
- [Fractalkine plasma estimates](../../results/fractalkine/cx3cl1_plasma_differential.csv)
- [Fractalkine RNA estimates](../../results/fractalkine/cx3cl1_tissue_rna_differential.csv)
- [WARS1 primary estimates](../../results/wars1/wars1_primary.csv)

Effects, raw p-values and pointwise CIs originate from MoTrPAC. The notebooks independently verified source BH corrections. They did not refit participant-level models. The deck's speaker notes contain additional statistical detail and result provenance.
