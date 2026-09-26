# Human MoTrPAC fractalkine / CX3CL1 findings

Collection c2.0; public Analysis package 2.0.8; commit `535b4044e7417413de471104c619120337602b77`.

- **Plasma:** 2 of 10 exercise-versus-control tests pass original source FDR < 0.05.
  - EE-CON, During 20 min: effect +0.5681, pointwise 95% CI [0.3452, 0.7910], raw p 1.29e-06, source adjusted p 0.000364.
  - EE-CON, During 40 min: effect +0.5605, pointwise 95% CI [0.3324, 0.7887], raw p 2.89e-06, source adjusted p 0.000227.
- **Temporal sensitivity:** 2 plasma tests pass Holm correction across the ten CX3CL1 exercise-versus-control comparisons.
- **Direct exercise-mode comparison:** 0 of 4 post-exercise plasma tests pass source FDR; minimum source adjusted p = 0.472.
- **Tissue RNA:** 4 of 12 available muscle/adipose exercise-versus-control tests pass source FDR.

## Interpretation

The human results support a transient circulating CX3CL1 response during endurance exercise, alongside post-exercise tissue RNA responses. During-resistance blood was not sampled, and the direct post-exercise plasma tests do not establish an endurance-specific response.

Tissue origin, secretion or cleavage, clearance, and disease benefit remain unresolved. A transcriptional response after exercise is not evidence that it generated the earlier plasma pulse. This cohort does not provide an islet outcome or a mediation analysis.

The next computational step is an independent disease-evidence overlay with explicit gene/receptor mapping and direction. A subsequent experiment could test a physiologically calibrated CX3CL1 pulse in human islets with and without CX3CR1 blockade, measuring survival and glucagon responses; that experiment has not been performed here.

## Statistical limits

These are published mixed-model results, not new fits to individual participants. Source FDR is per assay/contrast; pointwise CIs are not multiplicity-adjusted. CX3CL1 was selected after exploratory screening. Nonsignificance does not establish equivalence or no biological response. Model effects should not be pooled across protein and RNA scales.

## Verification

Nine source-file checksums verified; 75 source contrasts extracted; 32 previously reported rows matched; full plasma BH values independently reproduced for 19838 feature-by-contrast rows. Contrast identities and data-integrity checks passed.

## Direct tissue RNA differences

- muscle, Post 15 min: EE–RE effect +1.0407, source adjusted p 4.22e-12. This difference concerns tissue RNA and does not establish a plasma protein difference.