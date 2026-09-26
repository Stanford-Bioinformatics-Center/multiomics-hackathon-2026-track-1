# CCN1: human MoTrPAC c2.0 findings

Analyzed 105 published CCN1 contrasts. Primary analysis includes 30 exercise-versus-control tests and 14 direct EE–RE tests.

Independent BH recalculation matched 426,296 feature/contrast rows across 44 assay/contrast families.

## Control-adjusted results

| measurement | tests | source_BH_hits |
| --- | --- | --- |
| Adipose RNA | 6 | 2 |
| Adipose protein | 2 | 0 |
| Muscle RNA | 6 | 4 |
| Muscle protein | 6 | 1 |
| Plasma protein (Olink) | 10 | 2 |

All source-BH-significant control-adjusted CCN1 results:

| measurement | contrast_category | time_label | logFC | p_value | adj_p_value | holm_CCN1_control_tests |
| --- | --- | --- | --- | --- | --- | --- |
| Adipose RNA | EE-CON | Post 45 min | 1.6905 | 7.9856e-12 | 1.2404e-07 | 2.0762e-10 |
| Adipose RNA | RE-CON | Post 45 min | 1.0147 | 2.644e-05 | 0.023095 | 0.00063457 |
| Muscle RNA | EE-CON | Post 15 min | 4.194 | 2.3389e-30 | 2.2713e-27 | 6.7828e-29 |
| Muscle RNA | RE-CON | Post 15 min | 4.7509 | 1.2422e-37 | 2.9296e-34 | 3.7265e-36 |
| Muscle RNA | RE-CON | Post 3.5 h | 2.5859 | 1.9962e-15 | 9.9865e-14 | 5.5894e-14 |
| Muscle RNA | RE-CON | Post 24 h | 0.87078 | 0.0022992 | 0.015153 | 0.050583 |
| Muscle protein | RE-CON | Post 3.5 h | 0.92236 | 1.5047e-12 | 4.6729e-09 | 4.0628e-11 |
| Plasma protein (Olink) | EE-CON | During 40 min | 0.7642 | 0.0013578 | 0.022199 | 0.03123 |
| Plasma protein (Olink) | RE-CON | Post 10 min | 1.0172 | 6.493e-06 | 0.0010223 | 0.00016232 |

## Direct modality differences

| measurement | contrast_category | time_label | logFC | p_value | adj_p_value | holm_CCN1_direct_tests |
| --- | --- | --- | --- | --- | --- | --- |
| Muscle RNA | EE-RE | Post 3.5 h | -2.2182 | 2.6611e-17 | 2.4407e-15 | 3.7255e-16 |
| Muscle RNA | EE-RE | Post 24 h | -0.8581 | 0.00039264 | 0.0039817 | 0.0043191 |
| Muscle protein | EE-RE | Post 3.5 h | -0.80867 | 2.4879e-13 | 7.7262e-10 | 3.2343e-12 |

Holm sensitivity families contain 30 control-adjusted and 14 direct CCN1 tests, respectively. This is exploratory and does not account for earlier selection among candidate molecules.

8 of the 9 source-BH-significant control-adjusted results also pass the CCN1-specific Holm check. The muscle RNA response at 24 hours after resistance narrowly misses it (Holm p = 0.0506; source BH p = 0.0152). The 24-hour direct muscle EE–RE contrast passes both corrections.

Three direct comparisons pass the CCN1-only Holm check but not source BH: adipose RNA at 45 minutes, muscle protein at 24 hours, and plasma protein at 10 minutes. These remain exploratory under our primary assay-wide rule. Holm here uses a much smaller hypothesis family, so its values need not be larger than assay-wide BH values.

## Interpretation

CCN1 RNA increases relative to controls in both muscle and adipose during early recovery. Resistance-associated muscle RNA remains significant at 3.5 and 24 hours. Direct muscle RNA contrasts support a larger resistance response at those later times.

Plasma CCN1 protein passes source BH at 40 minutes during endurance and 10 minutes after resistance. None of the available direct plasma EE–RE contrasts passes source BH < 0.05. Thus the plasma data do not establish a modality difference at that threshold.

Muscle CCN1 protein increases after resistance at 3.5 hours, with a significant direct EE–RE contrast favoring resistance. Adipose protein at 4 hours is not significant; other adipose protein times are unavailable.

These are exercise-response findings, not evidence of diabetes prevention. The literature provides a rationale involving diabetic complications, including retinopathy, with potentially different effects from a transient exercise response. RNA induction, plasma abundance, and tissue protein do not establish secretion, source cell, causality, or target-organ activity.

## Reproduction and limits

MoTrPAC supplied raw p-values, effects, and confidence intervals. We checked inputs, mapped CCN1/CYR61, reproduced assay-wide BH from raw p-values, and applied exploratory CCN1-specific Holm corrections. Participant-level models were not refitted. The MoTrPAC preprint already highlights CCN1; this reanalysis is not independent replication.

The full executed notebook contains the evidence review, source links, complete nonsignificant results, figures, and testable hypotheses. See ../../ccn1_differential_analysis.ipynb and provenance.json.
