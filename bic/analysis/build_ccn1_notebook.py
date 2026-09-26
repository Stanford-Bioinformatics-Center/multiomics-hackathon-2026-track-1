"""Build the separate CCN1 biology and human differential-response notebook."""
from pathlib import Path
import textwrap
import nbformat as nbf

root = Path(__file__).resolve().parents[1]
cells = []
def md(s): cells.append(nbf.v4.new_markdown_cell(textwrap.dedent(s).strip()))
def code(s): cells.append(nbf.v4.new_code_cell(textwrap.dedent(s).strip()))

md("""
# CCN1 / CYR61: exerkine biology, diabetes relevance, and human exercise responses

**Question:** Is CCN1 relevant to diabetes, and how do its RNA and protein measurements respond to acute endurance and resistance exercise in humans?

**Yes, there is a research rationale.** Human exercise studies implicate muscle CCN1 expression; human and experimental studies link CCN1 to diabetic vascular complications. These findings justify examining its exercise response, but do not establish that increasing CCN1 prevents or treats diabetes.

We use the available **human MoTrPAC c2.0** release, **Acute Exercise in Human Sedentary Adults** (`human-precovid-sed-adu`), package **2.0.8** dated **2026-09-23**, frozen at commit `535b4044e7417413de471104c619120337602b77`. These are responses to a **single exercise bout**, not adaptations after weeks of training.

This notebook analyzes published differential-model results. MoTrPAC supplied the effect estimates, confidence intervals, raw p-values, and adjusted p-values. We extract CCN1, independently reproduce the assay-wide BH correction, compare time courses, and examine direct exercise-mode contrasts. We **do not fit a new participant-level model**: individual measurements are not in these public result files.

Run all cells to regenerate `results/ccn1/`. Inputs are cached locally; no account or download is required.

[Human collection](https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0) · [Pinned package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77) · [Version history](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/NEWS.md)
""")
md("""
## 1. What CCN1 does

**CCN1** means *cellular communication network factor 1*; **CYR61** is an older name for the same gene. It is distinct from **CCN2/CTGF**. Its protein is secreted and associates with the extracellular matrix, the material surrounding cells. As a **matricellular signaling protein**, it helps cells respond to that surrounding environment through receptors including integrins. Its functions include cell adhesion, migration, blood-vessel formation, and tissue remodeling. [NCBI human CCN1 record](https://www.ncbi.nlm.nih.gov/gene/3491).

CCN1's effects depend on cell type and context. In experimental skin wounds, CCN1 induced fibroblast senescence and limited excessive scarring; a mouse mutant unable to induce that response developed greater fibrosis. That is evidence for a tissue-repair mechanism, not an exercise or diabetes treatment result. [Jun & Lau, 2010](https://pubmed.ncbi.nlm.nih.gov/20526329/).

We call CCN1 a **candidate exerkine** because exercise-responsive tissue expression and circulating protein support a possible signaling role. An increase in muscle RNA alone does not show secretion into blood. A plasma increase does not identify the source tissue or prove an effect on another organ.

The MoTrPAC **2026 preprint** already highlights CCN1 as an exerkine candidate, using tissue RNA, tissue protein, and plasma responses. Our work is a focused reanalysis of the available release with a diabetes-oriented interpretation, not an independent discovery or replication of that MoTrPAC observation. [MoTrPAC multi-tissue preprint, Figure 8](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/).
""")
md("""
## 2. What connects CCN1 to diabetes?

**T2D** means type 2 diabetes. **Diabetic retinopathy** is damage to the retina's blood vessels associated with diabetes. **NETs** are extracellular webs of DNA and proteins released by neutrophils, a type of immune cell.

| Evidence | What was observed | What it supports—and its limit |
|---|---|---|
| Human exercise study: 13 men with T2D and 14 weight-matched glucose-tolerant men; CCN1 RNA available in 13 per group | Muscle CYR61 RNA rose after 60 minutes of cycling and returned toward baseline by 3 hours; no significant between-group difference in the response was detected | Exercise-responsive muscle transcription can occur with T2D. This does not demonstrate equal responses statistically, improved insulin sensitivity, or circulating CCN1 secretion. [Sabaratnam et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29924476/) |
| Human cross-sectional study: 50 healthy controls, 74 T2D without retinopathy, 69 with retinopathy | Plasma CCN1 was higher with diabetic retinopathy; healthy controls and T2D without retinopathy did not differ significantly | A complication-associated biomarker signal, not evidence that CCN1 causes T2D or prospectively predicts retinopathy. [Xiang et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10273100/) |
| Human observations plus diabetic mouse and cell experiments | CCN1 was associated with neutrophil accumulation/NET formation; reducing CCN1 alleviated retinal leakage in diabetic mice | Experimental support for a potentially harmful retinal mechanism. It does not establish that the acute exercise pulse has this effect in humans. [Li et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11097549/) |
| Human kidney tissue | Podocyte CCN1 expression was lower in diabetic nephropathy, as well as some other kidney diseases | Tissue and compartment matter; a direction seen in plasma cannot be assigned to every organ. [Sawai et al., 2007](https://pubmed.ncbi.nlm.nih.gov/17699553/) |

**Decision:** proceed with differential-response analysis. The strongest disease link reviewed here concerns **diabetic complications and tissue remodeling**, rather than demonstrated glucose lowering. A transient exercise response and a chronic disease-associated elevation could have different consequences; that distinction remains a hypothesis to test. This is a focused literature review, not a systematic review or a diabetes genetic-association analysis.
""")
md(r"""
## 3. Precisely what is tested?

For each tissue, assay, and timepoint, the primary effect is:

\[
\beta_{EE-CON,t}=(EE_t-EE_{baseline})-(CON_t-CON_{baseline}).
\]

The corresponding **null hypothesis is \(H_0:\beta_{EE-CON,t}=0\)**. The resistance comparison replaces EE with RE. A positive effect means a greater baseline-to-timepoint change than in resting controls; it does not alone establish an absolute increase within the exercise arm.

The direct EE–RE contrast tests whether the two exercise arms' changes differ. It is needed to claim a modality difference; significance versus control in only one arm is insufficient.

**Primary layer:** CCN1 RNA in the analyzed tissue tables. **Supporting layers:** plasma Olink protein and tissue proteomics. We retain every available timepoint, including nonsignificant comparisons, and keep assay scales separate.

**BH rule:** original assay-wide BH-adjusted p < 0.05. We recompute BH from *all features' raw p-values in the same tissue/assay/platform/contrast*, before selecting CCN1. This checks MoTrPAC's correction; it does not change the null hypothesis or generate new biological p-values. Source BH addresses each contrast separately, not the entire search across times, tissues, and modes.

**Additional exploratory check:** Holm adjustment across all available CCN1 exercise-versus-control tests together, and separately across all CCN1 direct EE–RE tests. This checks selection among the displayed CCN1 comparisons. It does not correct for choosing CCN1 after exploring other molecules or replace assay-wide BH. No p-value here tests diabetes protection.

Intervals are source **pointwise 95% confidence intervals**. They are not corrected for multiple comparisons. Effect magnitudes remain on each source model's `logFC` scale; we do not convert them to concentrations or pool RNA and protein effects.

The study includes 175 sedentary adults overall; the number differs by assay and sampling schedule. Counts of features below are **molecular features, not participants**. No diabetes-versus-control analysis is possible from these contrast tables. Participant-level correlations, mediation, and disease outcomes are not estimated.
""")
code(r"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from IPython.display import display, Markdown, Image

def locate_project():
    here = Path.cwd().resolve()
    for parent in [here, *here.parents]:
        for candidate in (parent, parent / 'BIC'):
            if (candidate / 'analysis/extract_ccn1.R').is_file():
                return candidate
    raise FileNotFoundError('Open this notebook from BIC or a subdirectory.')

ROOT = locate_project()
DATA, OUT = ROOT / 'motrpac-exploration', ROOT / 'results/ccn1'
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05
COMMIT = '535b4044e7417413de471104c619120337602b77'
RSCRIPT = shutil.which('Rscript') or '/usr/local/bin/Rscript'
assert Path(RSCRIPT).is_file(), 'Install base R and make Rscript available.'
pd.set_option('display.max_rows', 60)
pd.set_option('display.max_columns', 18)
pd.set_option('display.float_format', lambda x: f'{x:.5g}')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})
print('Project:', ROOT, '\nTarget: CCN1 / CYR61 | human collection c2.0')
""")
md("""
## 4. Verify inputs and identify available CCN1 measurements

The extractor checks eight cached human assay result objects and matches the source gene map using **CCN1 or CYR61**. It exports all CCN1 contrasts plus full feature backgrounds for assays with a CCN1 result. Epigenomic mappings are saved separately but are not analyzed here.

An assay feature ID identifies a measurement: `ENSG00000142871.18` is the versioned RNA gene identifier, `OID21368` is the plasma Olink assay identifier, and `O00622` is the tissue protein accession. A gene-map entry alone does not guarantee an analyzed result.
""")
code(r"""
SOURCE_OBJECTS = ['HUMAN_FEATURE_TO_GENE', 'BLOOD_PROT_OL_DA', 'BLOOD_TRNSCRPT_DA',
                  'MUSCLE_TRNSCRPT_DA', 'ADIPOSE_TRNSCRPT_DA', 'MUSCLE_PROT_PR_DA',
                  'MUSCLE_PROT_PH_DA', 'ADIPOSE_PROT_PR_DA', 'ADIPOSE_PROT_PH_DA']
expected = json.loads((DATA / 'table3/input_provenance.json').read_text())
assert expected['motrpac_commit'] == COMMIT
input_hashes = {}
for name in SOURCE_OBJECTS:
    path = DATA / f'{name}.rda'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == expected['input_sha256'][path.name], f'Input mismatch: {path.name}'
    input_hashes[path.name] = digest
print(f'Verified {len(input_hashes)} cached files against the pinned release.')
run = subprocess.run([RSCRIPT, str(ROOT / 'analysis/extract_ccn1.R'), str(DATA), str(OUT)],
                     check=True, text=True, capture_output=True)
print(run.stdout)
if run.stderr.strip(): print(run.stderr)
results = pd.read_csv(OUT / 'ccn1_source_contrasts.csv')
coverage = pd.read_csv(OUT / 'ccn1_assay_coverage.csv')
mapping = pd.read_csv(OUT / 'ccn1_feature_mapping.csv')
background = pd.read_csv(OUT / 'assay_background.csv.gz')
display(coverage)
display(mapping)
""")
md("""
**Coverage interpretation:** CCN1 has analyzed results for muscle/adipose RNA, plasma Olink protein, and muscle/adipose total protein. There is no CCN1 result in the analyzed blood RNA or phosphoprotein tables. This means **unavailable in these analyzed results**, not absent from the body. Adipose protein results are available only at the intermediate recovery timepoint; missing times are not filled with zero.

Sampling differs by tissue: the combined source label `post_15_30_45_min` means muscle **15 min**, plasma **30 min**, or adipose **45 min**. `post_3.5_4_hr` means muscle/plasma **3.5 h** or adipose **4 h**. EE and control have during-exercise blood samples; RE does not.
""")
code(r"""
keys = ['source_object', 'feature_id', 'contrast']
assert not results.duplicated(keys).any()
assert not background.duplicated(keys).any()
assert set(results.gene_symbol) <= {'CCN1', 'CYR61'}
for frame in [results, background]:
    assert np.isfinite(frame[['p_value', 'adj_p_value', 'logFC']]).all().all()
    assert frame.p_value.between(0, 1).all() and frame.adj_p_value.between(0, 1).all()
assert results['CI.L_calculated'].le(results.logFC).all()
assert results.logFC.le(results['CI.R_calculated']).all()
assert set(results[results.assay.eq('prot-ol')].feature_id) == {'OID21368'}
assert set(results[results.assay.eq('transcript-rna-seq')].feature_id) == {'ENSG00000142871.18'}
assert set(results[results.assay.eq('prot-pr')].feature_id) == {'O00622'}

TIMES = ['during_20_min', 'during_40_min', 'post_10_min',
         'post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
def time_label(row):
    if row.Timepoint == 'post_15_30_45_min':
        return 'Post ' + {'blood': '30', 'muscle': '15', 'adipose': '45'}[row.tissue] + ' min'
    if row.Timepoint == 'post_3.5_4_hr':
        return 'Post ' + ('4' if row.tissue == 'adipose' else '3.5') + ' h'
    return {'during_20_min': 'During 20 min', 'during_40_min': 'During 40 min',
            'post_10_min': 'Post 10 min', 'post_24_hr': 'Post 24 h', 'pre_exercise': 'Baseline'}[row.Timepoint]
results['time_label'] = results.apply(time_label, axis=1)
results['time_order'] = results.Timepoint.map({t: i for i, t in enumerate(TIMES)}).fillna(-1)
def layer(row):
    if row.assay == 'prot-ol': return 'Plasma protein (Olink)'
    return row.tissue.capitalize() + (' RNA' if row.assay == 'transcript-rna-seq' else ' protein')
results['measurement'] = results.apply(layer, axis=1)
results = results.sort_values(['measurement', 'contrast_category', 'time_order']).reset_index(drop=True)
print(f'Validated {len(results)} CCN1 source contrasts; these include baseline and within-arm contrasts.')
display(results[['measurement', 'full_model']].drop_duplicates().reset_index(drop=True))
""")
md("""
## 5. Reproduce the full-assay BH corrections

Within each comparison, BH sorts all raw p-values, adjusts for the number of tested features, and enforces monotonicity. We verify that the recomputed values agree with the published adjusted p-values. **We do not apply BH to values already adjusted by MoTrPAC.**

The audit table counts all features and all source FDR hits for each comparison. It describes assay-wide context. The later tables contain only CCN1. Holm is a separate sensitivity calculation from the CCN1 **raw** p-values.
""")
code(r"""
def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(adjusted, 1)
    return out

def holm_adjust(values):
    p = np.asarray(values, dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    adjusted = np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1))
    out = np.empty(len(p)); out[order] = np.minimum(adjusted, 1)
    return out

family = ['tissue', 'assay', 'platform', 'contrast']
background['bh_recomputed'] = background.groupby(family, observed=True, dropna=False).p_value.transform(bh_adjust)
np.testing.assert_allclose(background.bh_recomputed, background.adj_p_value, rtol=1e-10, atol=1e-12)
background['bh_abs_difference'] = (background.bh_recomputed - background.adj_p_value).abs()
background['p_rank_in_assay'] = background.groupby(family, observed=True, dropna=False).p_value.rank(method='min')
background['n_features_in_assay'] = background.groupby(family, observed=True, dropna=False).feature_id.transform('size')
audit = background.groupby(family + ['contrast_category', 'Timepoint'], observed=True, dropna=False).agg(
    n_features=('feature_id', 'size'), source_FDR_hits=('adj_p_value', lambda p: int((p < ALPHA).sum())),
    max_BH_difference=('bh_abs_difference', 'max')).reset_index()
audit.to_csv(OUT / 'bh_audit.csv', index=False)
print(f'BH verified for {len(background):,} feature/contrast rows in {len(audit)} comparison families.')
print('Maximum absolute difference:', background.bh_abs_difference.max())
display(audit.drop(columns=['contrast', 'platform', 'max_BH_difference']))

focused = results[results.contrast_type.isin(['exercise_with_controls', 'Endur_vs_Resist'])].copy()
focused = focused.merge(background[keys + ['bh_recomputed', 'p_rank_in_assay', 'n_features_in_assay']],
                        on=keys, how='left', validate='one_to_one')
assert focused.bh_recomputed.notna().all()
focused['source_significant'] = focused.adj_p_value < ALPHA
primary = focused[focused.contrast_type.eq('exercise_with_controls')].copy()
direct = focused[focused.contrast_type.eq('Endur_vs_Resist')].copy()
primary['holm_CCN1_control_tests'] = holm_adjust(primary.p_value)
direct['holm_CCN1_direct_tests'] = holm_adjust(direct.p_value)
print(f'Exploratory Holm families: {len(primary)} control-adjusted CCN1 tests; {len(direct)} direct CCN1 tests.')
""")
md("""
## 6. Primary differential-expression results: CCN1 RNA

These rows focus **only on CCN1**, with separate muscle and adipose results. Each raw p-value tests no control-adjusted expression change at that time. `adj_p_value` is the source assay-wide BH value. `holm_CCN1_control_tests` covers the combined CCN1 comparisons across all available RNA and protein layers.

`n_features_in_assay` is the number of molecular features used in that comparison's BH correction. It is not the number of CCN1 measurements or participants.
""")
code(r"""
TABLE = ['measurement', 'feature_id', 'contrast_category', 'time_label', 'logFC',
         'CI.L_calculated', 'CI.R_calculated', 'p_value', 'adj_p_value', 'source_significant']
rna = primary[primary.assay.eq('transcript-rna-seq')].copy()
plasma = primary[primary.assay.eq('prot-ol')].copy()
tissue_protein = primary[primary.assay.eq('prot-pr')].copy()
display(rna[TABLE + ['holm_CCN1_control_tests', 'n_features_in_assay']].reset_index(drop=True))
for name, frame in [('ccn1_rna_primary', rna), ('ccn1_plasma_primary', plasma),
                    ('ccn1_tissue_protein_primary', tissue_protein), ('ccn1_primary_all_layers', primary),
                    ('ccn1_direct_EE_RE', direct)]:
    frame.to_csv(OUT / f'{name}.csv', index=False)
""")
code(r"""
COLORS = {'EE-CON': '#2267a4', 'RE-CON': '#bd5c29'}
def save_figure(fig, stem):
    fig.savefig(OUT / f'{stem}.png', dpi=170, facecolor='white', bbox_inches='tight')
    fig.savefig(OUT / f'{stem}.svg', facecolor='white', bbox_inches='tight')
    display(Image(filename=str(OUT / f'{stem}.png')))
    plt.close(fig)

def response_panel(ax, frame, times, labels, title):
    ax.axhline(0, color='#777', linewidth=.8)
    for group, offset, label in [('EE-CON', -.10, 'Endurance vs control'), ('RE-CON', .10, 'Resistance vs control')]:
        q = frame[frame.contrast_category.eq(group)].copy()
        q['x'] = q.Timepoint.map({t: i for i, t in enumerate(times)})
        assert q.x.notna().all()
        q = q.sort_values('x')
        if q.empty: continue
        x, y = q.x.to_numpy() + offset, q.logFC.to_numpy()
        error = np.vstack([y - q['CI.L_calculated'].to_numpy(), q['CI.R_calculated'].to_numpy() - y])
        ax.errorbar(x, y, yerr=error, fmt='none', color=COLORS[group], capsize=3)
        ax.scatter(x, y, facecolors='white', edgecolors=COLORS[group], s=58, label=label, zorder=3)
        sig = q.source_significant.to_numpy()
        ax.scatter(x[sig], y[sig], color=COLORS[group], s=58, zorder=4)
    ax.set_xticks(range(len(times)), labels)
    ax.set_title(title, loc='left', fontweight='bold')
    ax.set_ylabel('Control-adjusted effect (source logFC)')
    ax.grid(axis='y', color='#e5e9ed'); ax.set_axisbelow(True)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), sharey=True)
for ax, tissue, labels in zip(axes, ['muscle', 'adipose'],
                            [['Post 15 min', 'Post 3.5 h', 'Post 24 h'], ['Post 45 min', 'Post 4 h', 'Post 24 h']]):
    response_panel(ax, rna[rna.tissue.eq(tissue)], TIMES[3:], labels, f'{tissue.capitalize()} CCN1 RNA')
axes[0].legend(frameon=False, fontsize=9)
fig.text(.02, .01, 'Filled: source BH-adjusted p < 0.05. Open: does not pass that threshold. Bars: pointwise 95% CIs.\nTimes are equally spaced; tissues were sampled at different recovery times.', fontsize=9)
fig.tight_layout(rect=(0, .10, 1, 1)); save_figure(fig, 'ccn1_rna_timecourse')
""")
md("""
**Read the RNA plots:** points are CCN1 effects, not separate participants. Muscle and adipose are separate panels. Filled points pass source BH < 0.05; open points do not. The zero line is no change relative to controls. Bar overlap is not a test of the difference between exercise modes; use the direct tests below. Tissue biopsies contain multiple cell types, so these RNA results cannot identify the cell type producing CCN1.

## 7. Supporting circulating protein results

This section concerns **plasma CCN1 protein**, measured by Olink `OID21368`. These are circulating relative-abundance measurements, not muscle RNA or absolute concentrations in pg/mL. A plasma rise can reflect release, distribution, clearance, or concentration effects; tissue origin is not resolved by these results.
""")
code(r"""
display(plasma[TABLE + ['holm_CCN1_control_tests']].reset_index(drop=True))
fig, ax = plt.subplots(figsize=(10.5, 4.9))
response_panel(ax, plasma, TIMES, ['During\n20 min', 'During\n40 min', 'Post\n10 min',
               'Post\n30 min', 'Post\n3.5 h', 'Post\n24 h'], 'Plasma CCN1 protein — Olink')
ax.legend(frameon=False, fontsize=9)
fig.text(.02, .01, 'Filled: source BH-adjusted p < 0.05. Bars: pointwise 95% CIs.\nNo during-resistance samples. Times are equally spaced for display.', fontsize=9)
fig.tight_layout(rect=(0, .10, 1, 1)); save_figure(fig, 'ccn1_plasma_timecourse')
""")
md("""
## 8. Supporting tissue protein results

Tissue proteomics measures CCN1 protein abundance in the biopsy, not secretion. Adipose has only one analyzed post-exercise time here. Early RNA induction followed by a later protein response is compatible with production, but group-level timing cannot prove that the RNA change caused the protein change.
""")
code(r"""
display(tissue_protein[TABLE + ['holm_CCN1_control_tests']].reset_index(drop=True))
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
response_panel(axes[0], tissue_protein[tissue_protein.tissue.eq('muscle')], TIMES[3:],
               ['Post 15 min', 'Post 3.5 h', 'Post 24 h'], 'Muscle CCN1 protein')
response_panel(axes[1], tissue_protein[tissue_protein.tissue.eq('adipose')], [TIMES[4]],
               ['Post 4 h'], 'Adipose CCN1 protein')
axes[0].legend(frameon=False, fontsize=9)
fig.text(.02, .01, 'Filled: source BH-adjusted p < 0.05. Bars: pointwise 95% CIs. Panels use separate y-axis ranges.\nAdipose protein has no analyzed early or 24-hour contrast in these source tables.', fontsize=9)
fig.tight_layout(rect=(0, .10, 1, 1)); save_figure(fig, 'ccn1_tissue_protein_timecourse')
""")
md("""
## 9. Does endurance differ directly from resistance?

Here the null is **equal baseline-to-timepoint changes in EE and RE**. Positive EE–RE effects favor a larger change with endurance; negative effects favor resistance. All available direct CCN1 contrasts are shown, including nonsignificant ones. The post-exercise clocks refer to recovery from each protocol; during-exercise EE–RE comparisons are unavailable.
""")
code(r"""
display(direct[TABLE + ['holm_CCN1_direct_tests']].reset_index(drop=True))
# A numerical identity check confirms that direct effects match the difference of control-adjusted effects.
wide = primary.pivot(index=['source_object', 'feature_id', 'Timepoint'], columns='contrast_category', values='logFC')
matched = direct.merge((wide['EE-CON'] - wide['RE-CON']).rename('EE_minus_RE').reset_index(),
                       on=['source_object', 'feature_id', 'Timepoint'], validate='one_to_one')
assert len(matched) == len(direct)
np.testing.assert_allclose(matched.logFC, matched.EE_minus_RE, rtol=1e-10, atol=1e-12)

fig, ax = plt.subplots(figsize=(10.5, 7))
q = direct.sort_values(['measurement', 'time_order']).reset_index(drop=True)
y = np.arange(len(q)); x = q.logFC.to_numpy()
err = np.vstack([x-q['CI.L_calculated'].to_numpy(), q['CI.R_calculated'].to_numpy()-x])
ax.errorbar(x, y, xerr=err, fmt='none', ecolor='#64748b', capsize=3)
ax.scatter(x, y, facecolors='white', edgecolors='#7255a5', s=52, zorder=3)
sig = q.source_significant.to_numpy()
ax.scatter(x[sig], y[sig], color='#7255a5', s=52, zorder=4)
ax.axvline(0, color='#777', linewidth=.8)
ax.set_yticks(y, q.measurement + ' | ' + q.time_label); ax.invert_yaxis()
ax.set_xlabel('EE − RE effect (source logFC; assay scales differ)')
ax.set_title('Direct CCN1 exercise-mode comparisons', loc='left', fontweight='bold')
ax.grid(axis='x', color='#e5e9ed'); ax.set_axisbelow(True)
fig.text(.02, .01, 'Negative = greater change with resistance. Filled = source BH p < 0.05; bars = pointwise 95% CIs.\nEach row is its own assay comparison; magnitudes across RNA and protein are not directly comparable.', fontsize=9)
fig.tight_layout(rect=(0, .08, 1, 1)); save_figure(fig, 'ccn1_direct_EE_RE')
""")
md("""
## 10. CCN1 in the full assay background

These four illustrative volcano panels show early responses for muscle RNA and plasma protein after each exercise mode; they do not replace the complete tables above. Every point is a feature in that assay, not a tissue or person. **Blue** means source BH-adjusted p < 0.05; **grey** means it does not pass. The red ring locates CCN1 regardless of significance. Height uses **raw p**, while color uses **BH-adjusted p**. All figures use the same fixed 0.05 threshold.
""")
code(r"""
panels = [('MUSCLE_TRNSCRPT_DA', 'EE-CON', 'post_15_30_45_min', 'Muscle RNA: endurance, post 15 min'),
          ('MUSCLE_TRNSCRPT_DA', 'RE-CON', 'post_15_30_45_min', 'Muscle RNA: resistance, post 15 min'),
          ('BLOOD_PROT_OL_DA', 'EE-CON', 'during_40_min', 'Plasma protein: endurance, during 40 min'),
          ('BLOOD_PROT_OL_DA', 'RE-CON', 'post_10_min', 'Plasma protein: resistance, post 10 min')]
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
for ax, (obj, group, timepoint, title) in zip(axes.flat, panels):
    q = background[background.source_object.eq(obj) & background.contrast_category.eq(group) & background.Timepoint.eq(timepoint)].copy()
    y = -np.log10(q.p_value.clip(lower=np.finfo(float).tiny)); sig = q.adj_p_value < ALPHA
    ax.scatter(q.loc[~sig, 'logFC'], y[~sig], color='#c7cdd4', s=10, alpha=.5)
    ax.scatter(q.loc[sig, 'logFC'], y[sig], color='#2267a4', s=12, alpha=.65)
    fid = results[results.source_object.eq(obj)].feature_id.unique()
    focal = q[q.feature_id.isin(fid)]
    assert len(focal) == 1
    focal = focal.iloc[0]; fy = -np.log10(max(focal.p_value, np.finfo(float).tiny))
    ax.scatter([focal.logFC], [fy], facecolors='none', edgecolors='#b63737', linewidths=2, s=110, zorder=5)
    ax.annotate('CCN1', (focal.logFC, fy), xytext=(7, 7), textcoords='offset points', color='#9e2828')
    ax.axvline(0, color='#aaa', linewidth=.7)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_xlabel('Source model effect (logFC)'); ax.set_ylabel('−log10(raw p-value)')
    ax.text(.03, .97, f'{sig.sum():,}/{len(q):,} features pass BH < 0.05\nCCN1 adjusted p = {focal.adj_p_value:.3g}',
            transform=ax.transAxes, ha='left', va='top', fontsize=9,
            bbox={'facecolor': 'white', 'alpha': .9, 'edgecolor': 'none'})
fig.suptitle('CCN1 within each human assay — exercise vs resting control', x=.02, ha='left', fontweight='bold')
fig.tight_layout(rect=(0, 0, 1, .95)); save_figure(fig, 'ccn1_volcano_context')
""")
md("""
## 11. Results and diabetes hypotheses

The following report is generated from the computed tables, using source BH < 0.05 as the main decision rule. “Significant” describes an exercise contrast only. It does not establish a diabetes mechanism, tissue source, or long-term benefit.
""")
code(r"""
def markdown_table(frame, columns):
    def val(x):
        if isinstance(x, (float, np.floating)): return f'{x:.5g}'
        return str(x).replace('|', '/')
    return '\n'.join(['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---']*len(columns)) + ' |'] +
                     ['| ' + ' | '.join(val(v) for v in row) + ' |' for row in frame[columns].itertuples(index=False, name=None)])

hits = primary[primary.source_significant]
direct_hits = direct[direct.source_significant]
counts = primary.groupby('measurement').agg(tests=('p_value', 'size'), source_BH_hits=('source_significant', 'sum')).reset_index()
lines = ['# CCN1: human MoTrPAC c2.0 findings', '',
         f'Analyzed {len(results)} published CCN1 contrasts. Primary analysis includes {len(primary)} exercise-versus-control tests and {len(direct)} direct EE–RE tests.', '',
         f'Independent BH recalculation matched {len(background):,} feature/contrast rows across {len(audit)} assay/contrast families.', '',
         '## Control-adjusted results', '', markdown_table(counts, list(counts.columns)), '',
         'All source-BH-significant control-adjusted CCN1 results:', '',
         markdown_table(hits, ['measurement', 'contrast_category', 'time_label', 'logFC', 'p_value', 'adj_p_value', 'holm_CCN1_control_tests']), '',
         '## Direct modality differences', '',
         markdown_table(direct_hits, ['measurement', 'contrast_category', 'time_label', 'logFC', 'p_value', 'adj_p_value', 'holm_CCN1_direct_tests']), '',
         f'Holm sensitivity families contain {len(primary)} control-adjusted and {len(direct)} direct CCN1 tests, respectively. This is exploratory and does not account for earlier selection among candidate molecules.', '',
         f'{int((hits.holm_CCN1_control_tests < ALPHA).sum())} of the {len(hits)} source-BH-significant control-adjusted results also pass the CCN1-specific Holm check. The muscle RNA response at 24 hours after resistance narrowly misses it (Holm p = 0.0506; source BH p = 0.0152). The 24-hour direct muscle EE–RE contrast passes both corrections.', '',
         'Three direct comparisons pass the CCN1-only Holm check but not source BH: adipose RNA at 45 minutes, muscle protein at 24 hours, and plasma protein at 10 minutes. These remain exploratory under our primary assay-wide rule. Holm here uses a much smaller hypothesis family, so its values need not be larger than assay-wide BH values.', '',
         '## Interpretation', '',
         'CCN1 RNA increases relative to controls in both muscle and adipose during early recovery. Resistance-associated muscle RNA remains significant at 3.5 and 24 hours. Direct muscle RNA contrasts support a larger resistance response at those later times.', '',
         'Plasma CCN1 protein passes source BH at 40 minutes during endurance and 10 minutes after resistance. None of the available direct plasma EE–RE contrasts passes source BH < 0.05. Thus the plasma data do not establish a modality difference at that threshold.', '',
         'Muscle CCN1 protein increases after resistance at 3.5 hours, with a significant direct EE–RE contrast favoring resistance. Adipose protein at 4 hours is not significant; other adipose protein times are unavailable.', '',
         'These are exercise-response findings, not evidence of diabetes prevention. The literature provides a rationale involving diabetic complications, including retinopathy, with potentially different effects from a transient exercise response. RNA induction, plasma abundance, and tissue protein do not establish secretion, source cell, causality, or target-organ activity.', '',
         '## Reproduction and limits', '',
         'MoTrPAC supplied raw p-values, effects, and confidence intervals. We checked inputs, mapped CCN1/CYR61, reproduced assay-wide BH from raw p-values, and applied exploratory CCN1-specific Holm corrections. Participant-level models were not refitted. The MoTrPAC preprint already highlights CCN1; this reanalysis is not independent replication.', '',
         'The full executed notebook contains the evidence review, source links, complete nonsignificant results, figures, and testable hypotheses. See ../../ccn1_differential_analysis.ipynb and provenance.json.', '']
report = '\n'.join(lines)
(OUT / 'FINDINGS.md').write_text(report)
display(Markdown(report))
""")
md("""
### Testable hypotheses supported by the available results

| Hypothesis | What the current data can establish | What would test the disease mechanism |
|---|---|---|
| CCN1 participates in the more sustained muscle remodeling response to resistance exercise | Late RNA/protein and direct EE–RE contrasts support a modality-associated response | Perturb CCN1 in a muscle/endothelial model under contraction-related stimuli; measure remodeling, capillary function, and insulin-stimulated glucose uptake. Loss of the proposed benefit with CCN1 depletion, restored by add-back, would support mediation. |
| Transient exercise-related CCN1 and persistent diabetes-associated CCN1 exposure have different vascular effects | MoTrPAC establishes a sampled acute response; it does not measure a chronic diabetes state | Compare matched-dose pulses with sustained CCN1 exposure in retinal endothelial/neutrophil systems under normal and high-glucose conditions. Measure barrier permeability and NET formation, with CCN1 blockade and appropriate vehicle/osmotic controls. Similar effects at matched exposure would weaken a duration-specific hypothesis. |
| Diabetes alters the coupling between tissue CCN1 RNA and circulating protein | Current public summaries cannot estimate within-person coupling or diabetes interactions | Use independent participants with and without T2D, matched exercise protocols and serial RNA/protein samples, with resting controls. Test the group × exercise × time interaction and participant-level associations; collect glucose-disposal outcomes. The 2018 RNA result does not answer protein coupling. |

These are proposals, not experiments completed here. The tissue and plasma sampling clocks differ, and the early plasma response precedes the first post-exercise biopsy. We cannot infer that newly induced muscle RNA caused the initial plasma rise. Bulk RNA may also reflect stromal or vascular cells. No tissue-to-plasma correlation, glucose association, disease-signature reversal, or genetic causal analysis is claimed.

The appropriate follow-up is to separate **exercise responsiveness**, **tissue source**, and **disease function** experimentally. A secretion study would need direct evidence such as tissue release measurements or appropriate cell-specific perturbation; a glucose-benefit claim requires metabolic outcomes.
""")
md("""
## 12. Save provenance

The record below includes input checksums, exact source model formulas, full hypothesis-family sizes, software versions, and output checksums. The extractor and builder are stored in `analysis/`; execution can be repeated with `python3 analysis/run_notebook.py ccn1_differential_analysis.ipynb` from BIC. The existing fractalkine notebook is separate.
""")
code(r"""
source_base = f'https://raw.githubusercontent.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/{COMMIT}/data/'
provenance = {
    'generated_utc': datetime.now(timezone.utc).isoformat(),
    'gene': 'CCN1', 'aliases_searched': ['CCN1', 'CYR61'],
    'study': 'Acute Exercise in Human Sedentary Adults', 'dataset_id': 'human-precovid-sed-adu',
    'collection': 'c2.0', 'package_version': '2.0.8', 'package_date': '2026-09-23', 'commit': COMMIT,
    'analysis_type': 'Reanalysis of published human differential-model summary statistics; no participant-level refit',
    'alpha': ALPHA, 'primary_rule': 'Source BH-adjusted p < 0.05 in each tissue/assay/platform/contrast',
    'BH_family_columns': family, 'BH_family_count': len(audit), 'BH_feature_contrast_rows': len(background),
    'BH_max_abs_difference': float(background.bh_abs_difference.max()),
    'Holm_control_CCN1_tests': len(primary), 'Holm_direct_CCN1_tests': len(direct),
    'CCN1_source_contrasts': len(results),
    'inputs': {name: {'sha256': digest, 'url': source_base + name} for name, digest in input_hashes.items()},
    'source_models': results[['source_object', 'full_model']].drop_duplicates().to_dict('records'),
    'software': {'python': platform.python_version(),
                 **{name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'matplotlib', 'nbformat', 'nbclient']},
                 'R': subprocess.check_output([RSCRIPT, '--version'], text=True).strip()},
    'literature_reviewed': ['https://www.ncbi.nlm.nih.gov/gene/3491', 'https://pubmed.ncbi.nlm.nih.gov/20526329/',
       'https://pubmed.ncbi.nlm.nih.gov/29924476/', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10273100/',
       'https://pmc.ncbi.nlm.nih.gov/articles/PMC11097549/', 'https://pubmed.ncbi.nlm.nih.gov/17699553/',
       'https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/'],
    'scripts': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [ROOT / 'analysis/extract_ccn1.R', ROOT / 'analysis/build_ccn1_notebook.py', ROOT / 'analysis/run_notebook.py']},
    'outputs': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(OUT.iterdir()) if path.is_file() and path.name != 'provenance.json'},
}
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print('Saved results to', OUT)
print('Primary CCN1 source BH hits:', int(primary.source_significant.sum()), '/', len(primary))
print('Direct CCN1 source BH hits:', int(direct.source_significant.sum()), '/', len(direct))
""")

notebook = nbf.v4.new_notebook(cells=cells)
notebook.metadata = {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                     'language_info': {'name': 'python', 'version': '3.9'}}
for i, cell in enumerate(cells):
    if cell.cell_type == 'code': compile(cell.source, f'ccn1-cell-{i}', 'exec')
nbf.validate(notebook)
path = root / 'ccn1_differential_analysis.ipynb'
nbf.write(notebook, path)
print(f'Built {len(cells)} cells: {path}')
