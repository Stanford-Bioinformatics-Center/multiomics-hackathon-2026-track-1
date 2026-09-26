"""Build the readable CX3CL1 notebook. Execution is handled by run_notebook.py."""
from pathlib import Path
import textwrap
import nbformat as nbf

root = Path(__file__).resolve().parents[1]
cells = []
def md(s): cells.append(nbf.v4.new_markdown_cell(textwrap.dedent(s).strip()))
def code(s): cells.append(nbf.v4.new_code_cell(textwrap.dedent(s).strip()))

md("""
# Fractalkine / CX3CL1: differential responses to human acute exercise

**Question:** Does circulating fractalkine respond to endurance or resistance exercise, when does it respond, and do muscle/adipose CX3CL1 RNA results provide compatible tissue context?

**Dataset:** MoTrPAC **Acute Exercise in Human Sedentary Adults**, `human-precovid-sed-adu`, **collection c2.0**. Public analysis package **2.0.8**, dated **23 September 2026**, pinned to commit `535b4044e7417413de471104c619120337602b77`.

This notebook completes an analysis of the **published human differential-model results**. It extracts original contrasts, independently recalculates assay-wide plasma BH adjustments, evaluates effects and uncertainty, and checks direct endurance–resistance differences. It does **not** refit a participant-level model: individual measurements and phenotypes require approved access.

All inputs are cached locally. Run all cells to regenerate the tables, figures, findings, and provenance in `results/fractalkine/`.

[Human data collection](https://motrpac-data.org/data-download/file-browser/analysis/human-precovid-sed-adu/c2.0) · [Pinned public package](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77) · [Collection version history](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/NEWS.md)
""")
md("""
## 1. Study design and analysis decisions

The study includes 175 sedentary adults across all assays, assigned to endurance exercise (**EE**), resistance exercise (**RE**), or resting control (**CON**). This release focuses on one **acute bout**, before long-term training adaptation. Sample counts vary by assay and sampling schedule; 175 is not the sample size for every contrast.

| Analysis | Estimand | Role |
|---|---|---|
| EE–CON / RE–CON | Exercise-arm change from baseline minus time-matched control change | Primary response analysis |
| EE–RE | Difference between changes in the two exercise arms | Direct modality comparison |
| EE–EE / RE–RE / CON–CON | Within-arm change from baseline | Descriptive check of the role of controls |

**Primary significance:** original MoTrPAC BH-adjusted p < 0.05, calculated across the assay's features within each contrast. The notebook checks those values independently against all plasma Olink features before selecting CX3CL1.

**Sensitivity analysis:** Holm correction across the ten CX3CL1 plasma exercise-versus-control tests. This addresses the temporal/mode search within the selected molecule, but is exploratory because CX3CL1 was selected after an earlier screen. It does not replace the assay-wide source FDR or establish independent replication.

**Timing:** blood/plasma early recovery = 30 min, muscle = 15 min, adipose = 45 min; intermediate recovery = 3.5 h in plasma/muscle and 4 h in adipose. During-exercise samples exist for EE and controls, not RE. Post-exercise times are measured from the end of each protocol.

Confidence intervals are **pointwise 95% intervals**, not adjusted for multiple testing. Effects remain on each source model's scale; their magnitudes are not pooled across protein and RNA assays.
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
            if (candidate / 'analysis/extract_fractalkine.R').is_file() and (candidate / 'motrpac-exploration').is_dir():
                return candidate
    raise FileNotFoundError('Open the notebook from the BIC directory or one of its subdirectories.')

ROOT = locate_project()
DATA = ROOT / 'motrpac-exploration'
OUT = ROOT / 'results/fractalkine'
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05
GENE = 'CX3CL1'
COMMIT = '535b4044e7417413de471104c619120337602b77'
RSCRIPT = shutil.which('Rscript')
if RSCRIPT is None and Path('/usr/local/bin/Rscript').is_file():
    RSCRIPT = '/usr/local/bin/Rscript'
if RSCRIPT is None:
    raise RuntimeError('Base R is required: install R and make Rscript available on PATH.')

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.float_format', lambda x: f'{x:.5g}')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.spines.top': False, 'axes.spines.right': False})
print('Project:', ROOT)
print('Human collection: c2.0 | gene:', GENE)
print('Python:', sys.version.split()[0], '| Rscript:', RSCRIPT)
""")
md("""
## 2. Verify the frozen human inputs and extract CX3CL1

Input SHA-256 checksums must match the earlier pinned download. The base-R extractor reads eight human RNA/protein/phosphoprotein result objects and the feature-to-gene map. It retains every available CX3CL1 contrast, including nonsignificant results. It also exports the full plasma Olink universe for independent BH recalculation and volcano plots.
""")
code(r"""
SOURCE_OBJECTS = [
    'HUMAN_FEATURE_TO_GENE', 'BLOOD_PROT_OL_DA', 'BLOOD_TRNSCRPT_DA',
    'MUSCLE_TRNSCRPT_DA', 'ADIPOSE_TRNSCRPT_DA',
    'MUSCLE_PROT_PR_DA', 'MUSCLE_PROT_PH_DA',
    'ADIPOSE_PROT_PR_DA', 'ADIPOSE_PROT_PH_DA',
]
expected = json.loads((DATA / 'table3/input_provenance.json').read_text())
assert expected['motrpac_commit'] == COMMIT
input_hashes = {}
for name in SOURCE_OBJECTS:
    path = DATA / f'{name}.rda'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == expected['input_sha256'][path.name], f'Input mismatch: {path.name}'
    input_hashes[path.name] = digest
print(f'Verified {len(input_hashes)} source files against the pinned human release.')

run = subprocess.run([RSCRIPT, str(ROOT / 'analysis/extract_fractalkine.R'), str(DATA), str(OUT)],
                     check=True, text=True, capture_output=True)
print(run.stdout)
if run.stderr.strip():
    print(run.stderr)

results = pd.read_csv(OUT / 'cx3cl1_source_contrasts.csv')
coverage = pd.read_csv(OUT / 'cx3cl1_assay_coverage.csv')
mapping = pd.read_csv(OUT / 'cx3cl1_feature_mapping.csv')
background = pd.read_csv(OUT / 'plasma_olink_background.csv')
display(coverage)
display(mapping)
""")
md("""
An analyzed-feature count of zero means there is **no CX3CL1 result in that analyzed assay table**. It does not establish biological absence. The plasma protein should map to Olink assay **OID20976**; the tissue RNA should map to CX3CL1's Ensembl identifier. Blood-cell RNA and plasma protein are separate measurement layers.
""")
code(r"""
keys = ['source_object', 'feature_id', 'contrast_type', 'contrast_category', 'Timepoint']
assert not results.duplicated(keys).any()
assert not background.duplicated(keys).any()
assert set(results.gene_symbol) == {GENE}
assert results[['p_value', 'adj_p_value', 'logFC', 'CI.L_calculated', 'CI.R_calculated']].notna().all().all()
assert results.p_value.between(0, 1).all() and results.adj_p_value.between(0, 1).all()
assert results['CI.L_calculated'].le(results.logFC).all()
assert results.logFC.le(results['CI.R_calculated']).all()

TIMES = ['during_20_min', 'during_40_min', 'post_10_min',
         'post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
def time_label(row):
    if row.Timepoint == 'post_15_30_45_min':
        return 'Post ' + {'blood': '30', 'muscle': '15', 'adipose': '45'}[row.tissue] + ' min'
    if row.Timepoint == 'post_3.5_4_hr':
        return 'Post ' + ('4' if row.tissue == 'adipose' else '3.5') + ' h'
    return {'during_20_min': 'During 20 min', 'during_40_min': 'During 40 min',
            'post_10_min': 'Post 10 min', 'post_24_hr': 'Post 24 h',
            'pre_exercise': 'Baseline'}[row.Timepoint]

results['time_label'] = results.apply(time_label, axis=1)
results['source_significant'] = results.adj_p_value < ALPHA
results['time_order'] = results.Timepoint.map({t: i for i, t in enumerate(TIMES)})
primary = results[results.contrast_type.eq('exercise_with_controls')].copy()
plasma = primary[primary.assay.eq('prot-ol')].sort_values(['contrast_category', 'time_order']).copy()
direct = results[results.contrast_type.eq('Endur_vs_Resist')].copy()
plasma_direct = direct[direct.assay.eq('prot-ol')].sort_values('time_order').copy()
rna = primary[primary.assay.eq('transcript-rna-seq')].copy()
assert len(plasma) == 10 and len(plasma_direct) == 4
assert set(plasma.feature_id) == {'OID20976'}
assert set(plasma.loc[plasma.contrast_category.eq('EE-CON'), 'Timepoint']) == set(TIMES)
assert set(plasma.loc[plasma.contrast_category.eq('RE-CON'), 'Timepoint']) == set(TIMES[2:])

# Match the earlier screen by feature, source object, contrast, and time.
earlier = pd.read_csv(DATA / 'table3/all_candidate_results.csv')
earlier = earlier[earlier.candidate_name.eq('Fractalkine')]
audit = results.merge(earlier, on=keys, suffixes=('_new', '_earlier'), validate='one_to_one')
assert len(audit) == len(earlier)
for column in ['logFC', 'p_value', 'adj_p_value', 'CI.L_calculated', 'CI.R_calculated']:
    np.testing.assert_allclose(audit[column + '_new'], audit[column + '_earlier'], rtol=1e-12, atol=1e-15)
print(f'Validated {len(results)} source contrasts; {len(audit)} match the earlier screen.')
display(results.groupby(['tissue', 'assay', 'contrast_type']).size().rename('rows').reset_index())
""")
md("""
## 3. Independently recalculate the plasma multiple-testing adjustment

BH is recalculated using **all analyzed Olink features for each comparison**, before focusing on CX3CL1. Recalculating BH from only one chosen protein would answer a different statistical question. This check must reproduce the stored source adjusted p-values.
""")
code(r"""
def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1)
    return out

def holm_adjust(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind='stable')
    adjusted = np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1))
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1)
    return out

family = ['tissue', 'assay', 'platform', 'contrast']
background['bh_recomputed'] = background.groupby(family, observed=True, dropna=False).p_value.transform(bh_adjust)
np.testing.assert_allclose(background.bh_recomputed, background.adj_p_value, rtol=1e-10, atol=1e-12)
background['p_rank_in_assay'] = background.groupby(family, observed=True, dropna=False).p_value.rank(method='min')
background['n_features_in_assay'] = background.groupby(family, observed=True, dropna=False).feature_id.transform('size')
checks = background.groupby(['contrast_category', 'Timepoint']).agg(
    n_features=('feature_id', 'size'), source_FDR_hits=('adj_p_value', lambda p: int((p < ALPHA).sum())))
print('Maximum absolute BH difference:', float((background.bh_recomputed - background.adj_p_value).abs().max()))
display(checks.reset_index())

plasma = plasma.merge(background[['feature_id', 'contrast', 'bh_recomputed', 'p_rank_in_assay', 'n_features_in_assay']],
                      on=['feature_id', 'contrast'], validate='one_to_one')
plasma['holm_10_plasma_tests'] = holm_adjust(plasma.p_value)
plasma_direct['holm_4_direct_tests'] = holm_adjust(plasma_direct.p_value)
""")
md("""
## 4. Primary plasma differential results

Positive effects mean a greater change than in resting controls; negative effects mean a lower change. The source field is named `logFC`; it is a model-scale estimate, not an absolute plasma concentration. Pointwise intervals and nominal p-values describe uncertainty before multiplicity correction.
""")
code(r"""
columns = ['contrast_category', 'time_label', 'logFC', 'CI.L_calculated', 'CI.R_calculated',
           'p_value', 'adj_p_value', 'holm_10_plasma_tests', 'p_rank_in_assay', 'n_features_in_assay']
display(plasma[columns])
plasma.to_csv(OUT / 'cx3cl1_plasma_differential.csv', index=False)
hits = plasma[plasma.source_significant]
print('Primary source-FDR hits:', len(hits), 'of', len(plasma), 'tests')
print('Holm sensitivity hits:', int((plasma.holm_10_plasma_tests < ALPHA).sum()))
display(hits[columns])
""")
code(r"""
COLORS = {'EE-CON': '#2267a4', 'RE-CON': '#bd5c29'}
def save_figure(fig, stem):
    fig.savefig(OUT / f'{stem}.png', dpi=180, facecolor='white', bbox_inches='tight')
    fig.savefig(OUT / f'{stem}.svg', facecolor='white', bbox_inches='tight')
    display(Image(filename=str(OUT / f'{stem}.png')))
    plt.close(fig)

def response_panel(ax, frame, timepoints, labels, title):
    ax.axhline(0, color='#777', linewidth=.9)
    for group, offset, label in [('EE-CON', -.11, 'Endurance vs control'), ('RE-CON', .11, 'Resistance vs control')]:
        q = frame[frame.contrast_category.eq(group)].copy()
        q['x'] = q.Timepoint.map({t: i for i, t in enumerate(timepoints)})
        q = q.sort_values('x')
        if q.empty:
            continue
        x, y = q.x.to_numpy() + offset, q.logFC.to_numpy()
        error = np.vstack([y - q['CI.L_calculated'].to_numpy(), q['CI.R_calculated'].to_numpy() - y])
        ax.errorbar(x, y, yerr=error, fmt='none', color=COLORS[group], capsize=3)
        ax.scatter(x, y, facecolors='white', edgecolors=COLORS[group], s=62, label=label, zorder=3)
        sig = q.source_significant.to_numpy()
        ax.scatter(x[sig], y[sig], color=COLORS[group], s=62, zorder=4)
    ax.set_xticks(range(len(timepoints)), labels)
    ax.set_title(title, loc='left', fontweight='bold')
    ax.set_ylabel('Control-adjusted effect (source logFC)')
    ax.grid(axis='y', color='#e5e9ed')
    ax.set_axisbelow(True)

fig, ax = plt.subplots(figsize=(11, 5.2))
response_panel(ax, plasma, TIMES, ['During\n20 min', 'During\n40 min', 'Post\n10 min', 'Post\n30 min', 'Post\n3.5 h', 'Post\n24 h'],
               'Human plasma fractalkine / CX3CL1')
ax.legend(frameon=False)
fig.text(.01, .01, 'Filled points: original BH-adjusted p < 0.05. Bars: pointwise 95% CIs.\nTimepoints are equally spaced for display. No during-resistance samples.', fontsize=9, color='#555')
fig.tight_layout(rect=(0, .10, 1, 1))
save_figure(fig, 'cx3cl1_plasma_timecourse')
""")
md("""
## 5. Where CX3CL1 sits in the full plasma assay

These volcano plots retain every tested Olink feature. The x-axis is the source effect and the y-axis is −log10(raw p). Blue points pass the original assay-wide FDR; the red outlined point is CX3CL1. The panels illustrate the during-endurance result and matched post-exercise comparisons. They are descriptive views, not additional independent tests.
""")
code(r"""
panels = [('EE-CON', 'during_40_min', 'Endurance vs control: during 40 min'),
          ('EE-CON', 'post_10_min', 'Endurance vs control: post 10 min'),
          ('RE-CON', 'post_10_min', 'Resistance vs control: post 10 min'),
          ('EE-RE', 'post_10_min', 'Endurance vs resistance: post 10 min')]
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
for ax, (group, timepoint, title) in zip(axes.flat, panels):
    q = background[background.contrast_category.eq(group) & background.Timepoint.eq(timepoint)].copy()
    y = -np.log10(q.p_value.clip(lower=np.finfo(float).tiny))
    sig = q.adj_p_value < ALPHA
    ax.scatter(q.loc[~sig, 'logFC'], y[~sig], color='#c7cdd4', s=13, alpha=.7)
    ax.scatter(q.loc[sig, 'logFC'], y[sig], color='#2267a4', s=16, alpha=.8)
    focal = q[q.feature_id.eq('OID20976')].iloc[0]
    fy = -np.log10(focal.p_value)
    ax.scatter([focal.logFC], [fy], facecolors='none', edgecolors='#b63737', linewidths=2, s=105, zorder=5)
    ax.annotate('CX3CL1', (focal.logFC, fy), xytext=(7, 10), textcoords='offset points', color='#9e2828', fontsize=10)
    ax.axvline(0, color='#aaa', linewidth=.7)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_xlabel('Source model effect (logFC)')
    ax.set_ylabel('−log10(raw p-value)')
    ax.text(.03, .97, f'{int(sig.sum())}/{len(q)} assay features pass FDR\nCX3CL1 adjusted p = {focal.adj_p_value:.3g}',
            transform=ax.transAxes, ha='left', va='top', fontsize=9,
            bbox={'facecolor': 'white', 'alpha': .85, 'edgecolor': 'none'})
fig.suptitle('Fractalkine in the human plasma Olink differential analysis', x=.02, ha='left', fontweight='bold', fontsize=15)
fig.tight_layout(rect=(0, 0, 1, .95))
save_figure(fig, 'cx3cl1_plasma_volcano_context')
""")
md("""
## 6. Direct endurance-versus-resistance test

An EE–RE effect tests whether the changes differ between protocols. A significant EE–CON result alongside a nonsignificant RE–CON result does not itself establish that difference. No during-exercise EE–RE comparison is possible in this release.
""")
code(r"""
display(plasma_direct[['time_label', 'logFC', 'CI.L_calculated', 'CI.R_calculated', 'p_value', 'adj_p_value', 'holm_4_direct_tests']])
plasma_direct.to_csv(OUT / 'cx3cl1_plasma_endurance_vs_resistance.csv', index=False)
fig, ax = plt.subplots(figsize=(9, 4.7))
q = plasma_direct
x, y = np.arange(len(q)), q.logFC.to_numpy()
err = np.vstack([y - q['CI.L_calculated'].to_numpy(), q['CI.R_calculated'].to_numpy() - y])
ax.errorbar(x, y, yerr=err, fmt='o', markerfacecolor='white', color='#6d4a8e', capsize=4)
ax.axhline(0, color='#777', linewidth=.9)
ax.set_xticks(x, q.time_label)
ax.set_ylabel('Difference in changes: EE − RE')
ax.set_title('Plasma CX3CL1: direct comparison of exercise protocols', loc='left', fontweight='bold')
ax.grid(axis='y', color='#e5e9ed')
fig.text(.01, .01, 'Pointwise 95% CIs; post-exercise times are relative to the end of each protocol.', fontsize=9, color='#555')
fig.tight_layout(rect=(0, .07, 1, 1))
save_figure(fig, 'cx3cl1_direct_mode_comparison')
""")
md("""
## 7. Why control-adjusted and simple pre/post results differ

This table compares the source within-arm tests with the primary exercise-versus-control tests. They estimate different quantities and have different uncertainty and multiple-testing families. A pre/post p-value cannot replace the control-adjusted result. The effect identity below checks that the source contrast estimates are internally consistent; no independence or covariance is inferred from it.
""")
code(r"""
within = results[results.assay.eq('prot-ol') & results.contrast_type.eq('exercise_no_controls')].copy()
control = results[results.assay.eq('prot-ol') & results.contrast_type.eq('control_only')].copy()
rows = []
for _, r in plasma.iterrows():
    arm = 'EE-EE' if r.contrast_category == 'EE-CON' else 'RE-RE'
    w = within[within.contrast_category.eq(arm) & within.Timepoint.eq(r.Timepoint)].iloc[0]
    c = control[control.Timepoint.eq(r.Timepoint)].iloc[0]
    np.testing.assert_allclose(w.logFC - c.logFC, r.logFC, rtol=1e-10, atol=1e-12)
    rows.append({'exercise': r.contrast_category[:2], 'time': r.time_label,
                 'within_arm_effect': w.logFC, 'within_arm_source_q': w.adj_p_value,
                 'control_change': c.logFC, 'control_adjusted_effect': r.logFC,
                 'control_adjusted_source_q': r.adj_p_value})
control_check = pd.DataFrame(rows)
display(control_check)
control_check.to_csv(OUT / 'cx3cl1_within_arm_vs_control_adjusted.csv', index=False)
""")
md("""
## 8. Muscle and adipose RNA context

These tests measure **CX3CL1 gene expression within sampled tissue**, not release of protein into plasma. Source-FDR results are shown without adding a new cross-tissue significance claim. Tissue biopsies begin after exercise, while the significant plasma measurements occur during endurance exercise; this timing cannot establish new tissue transcription as the cause of the initial circulating pulse.
""")
code(r"""
rna = rna.sort_values(['tissue', 'contrast_category', 'time_order'])
display(rna[['tissue', 'feature_id', 'contrast_category', 'time_label', 'logFC', 'CI.L_calculated', 'CI.R_calculated', 'p_value', 'adj_p_value']])
rna.to_csv(OUT / 'cx3cl1_tissue_rna_differential.csv', index=False)
rna_direct = direct[direct.assay.eq('transcript-rna-seq')].sort_values(['tissue', 'time_order']).copy()
display(Markdown('**Direct EE–RE tissue RNA comparisons.** Positive effects mean a larger increase after endurance. These tests concern tissue RNA, not plasma protein.'))
display(rna_direct[['tissue', 'time_label', 'logFC', 'CI.L_calculated', 'CI.R_calculated', 'p_value', 'adj_p_value']])
rna_direct.to_csv(OUT / 'cx3cl1_tissue_rna_endurance_vs_resistance.csv', index=False)
direct.to_csv(OUT / 'cx3cl1_all_direct_mode_contrasts.csv', index=False)
rna_times = ['post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
fig, axes = plt.subplots(1, 2, figsize=(12, 5.3))
for ax, tissue, labels in zip(axes, ['muscle', 'adipose'],
                               [['Post\n15 min', 'Post\n3.5 h', 'Post\n24 h'], ['Post\n45 min', 'Post\n4 h', 'Post\n24 h']]):
    response_panel(ax, rna[rna.tissue.eq(tissue)], rna_times, labels, f'CX3CL1 RNA: {tissue}')
axes[0].legend(frameon=False, fontsize=9)
fig.text(.01, .01, 'Filled points: source BH-adjusted p < 0.05. Bars: pointwise 95% CIs.\nRNA is tissue context and does not establish circulating secretion. Panels have separate y-axis ranges.', fontsize=9, color='#555')
fig.tight_layout(rect=(0, .12, 1, 1))
save_figure(fig, 'cx3cl1_tissue_rna_timecourse')
""")
md("""
## 9. Findings, limits, and next experiment

The following summary is calculated from the results above, rather than copied from the earlier candidate screen. The source-FDR threshold remains primary. The focused Holm calculation is a sensitivity check, not a new confirmatory study.
""")
code(r"""
source_hits = plasma[plasma.source_significant].sort_values('time_order')
direct_hits = plasma_direct[plasma_direct.source_significant]
rna_hits = rna[rna.source_significant]
lines = [
    '# Human MoTrPAC fractalkine / CX3CL1 findings', '',
    f'Collection c2.0; public Analysis package 2.0.8; commit `{COMMIT}`.', '',
    f'- **Plasma:** {len(source_hits)} of {len(plasma)} exercise-versus-control tests pass original source FDR < {ALPHA}.',
]
for _, r in source_hits.iterrows():
    lines.append(f'  - {r.contrast_category}, {r.time_label}: effect {r.logFC:+.4f}, pointwise 95% CI [{r["CI.L_calculated"]:.4f}, {r["CI.R_calculated"]:.4f}], raw p {r.p_value:.3g}, source adjusted p {r.adj_p_value:.3g}.')
lines += [
    f'- **Temporal sensitivity:** {int((plasma.holm_10_plasma_tests < ALPHA).sum())} plasma tests pass Holm correction across the ten CX3CL1 exercise-versus-control comparisons.',
    f'- **Direct exercise-mode comparison:** {len(direct_hits)} of {len(plasma_direct)} post-exercise plasma tests pass source FDR; minimum source adjusted p = {plasma_direct.adj_p_value.min():.3g}.',
    f'- **Tissue RNA:** {len(rna_hits)} of {len(rna)} available muscle/adipose exercise-versus-control tests pass source FDR.', '',
    '## Interpretation', '',
    'The human results support a transient circulating CX3CL1 response during endurance exercise, alongside post-exercise tissue RNA responses. During-resistance blood was not sampled, and the direct post-exercise plasma tests do not establish an endurance-specific response.', '',
    'Tissue origin, secretion or cleavage, clearance, and disease benefit remain unresolved. A transcriptional response after exercise is not evidence that it generated the earlier plasma pulse. This cohort does not provide an islet outcome or a mediation analysis.', '',
    'The next computational step is an independent disease-evidence overlay with explicit gene/receptor mapping and direction. A subsequent experiment could test a physiologically calibrated CX3CL1 pulse in human islets with and without CX3CR1 blockade, measuring survival and glucagon responses; that experiment has not been performed here.', '',
    '## Statistical limits', '',
    'These are published mixed-model results, not new fits to individual participants. Source FDR is per assay/contrast; pointwise CIs are not multiplicity-adjusted. CX3CL1 was selected after exploratory screening. Nonsignificance does not establish equivalence or no biological response. Model effects should not be pooled across protein and RNA scales.', '',
    '## Verification', '',
    f'Nine source-file checksums verified; {len(results)} source contrasts extracted; {len(audit)} previously reported rows matched; full plasma BH values independently reproduced for {len(background)} feature-by-contrast rows. Contrast identities and data-integrity checks passed.',
]
rna_mode_hits = rna_direct[rna_direct.source_significant]
lines += ['', '## Direct tissue RNA differences', '']
for _, r in rna_mode_hits.iterrows():
    lines.append(f'- {r.tissue}, {r.time_label}: EE–RE effect {r.logFC:+.4f}, source adjusted p {r.adj_p_value:.3g}. This difference concerns tissue RNA and does not establish a plasma protein difference.')
findings = '\n'.join(lines)
(OUT / 'FINDINGS.md').write_text(findings)
display(Markdown(findings))
""")
md("""
## 10. Provenance and reproducibility

Outputs retain all estimates, including nonsignificant findings. The original cached `.rda` files are never modified. The run records checksums, model formulas, package versions, and the exact test families. A future participant-level analysis would require approved normalized matrices and phenotype/sample metadata; it should model repeated measures and appropriate covariates rather than treat timepoints as independent participants.

Sources: [MoTrPAC study design](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/), [source contrast definitions](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/blob/535b4044e7417413de471104c619120337602b77/R/load_differential_analysis.R), [CX3CL1 human-islet study for future disease interpretation](https://pmc.ncbi.nlm.nih.gov/articles/PMC4209359/). The exercise analysis does not independently validate the disease mechanism.
""")
code(r"""
results.to_csv(OUT / 'cx3cl1_all_contrasts_annotated.csv', index=False)
background.to_csv(OUT / 'plasma_olink_background_checked.csv', index=False)
provenance = {
    'study': 'human-precovid-sed-adu', 'species': 'Homo sapiens',
    'collection': 'c2.0', 'package': 'MotrpacHumanPreSuspensionAnalysis',
    'package_version': '2.0.8', 'commit': COMMIT,
    'run_utc': datetime.now(timezone.utc).isoformat(),
    'analysis': 'Published human differential-model results; no participant-level model refit',
    'input_sha256': input_hashes,
    'source_model_formulas': sorted(results.full_model.unique().tolist()),
    'primary_family': 'BH across all plasma Olink features within each source contrast',
    'sensitivity_family': 'Holm across ten CX3CL1 plasma EE-CON / RE-CON tests; exploratory',
    'direct_sensitivity_family': 'Holm across four CX3CL1 plasma EE-RE tests; separate exploratory family',
    'source_contrast_rows': len(results), 'plasma_background_rows': len(background),
    'primary_plasma_tests': len(plasma), 'direct_plasma_tests': len(plasma_direct),
    'tissue_rna_tests': len(rna), 'verified_previous_rows': len(audit),
    'maximum_BH_absolute_difference': float((background.bh_recomputed - background.adj_p_value).abs().max()),
    'software': {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'matplotlib', 'nbformat', 'nbclient', 'ipykernel']},
    'python': sys.version, 'platform': platform.platform(),
    'R': subprocess.run([RSCRIPT, '--version'], text=True, capture_output=True, check=True).stdout.strip(),
    'extractor_sha256': hashlib.sha256((ROOT / 'analysis/extract_fractalkine.R').read_bytes()).hexdigest(),
}
provenance['output_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'provenance.json'}
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2))
display(pd.DataFrame([{'file': p.name, 'bytes': p.stat().st_size} for p in sorted(OUT.iterdir()) if p.is_file()]))
print('Analysis complete. Outputs:', OUT)
""")

notebook = nbf.v4.new_notebook(cells=cells)
notebook.metadata = {
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python', 'version': '3.11'},
    'bic': {'study': 'human-precovid-sed-adu', 'collection': 'c2.0', 'gene': 'CX3CL1',
            'commit': '535b4044e7417413de471104c619120337602b77'},
}
nbf.validate(notebook)
for index, cell in enumerate(notebook.cells):
    if cell.cell_type == 'code':
        compile(cell.source, f'notebook cell {index}', 'exec')
target = root / 'fractalkine_differential_analysis.ipynb'
nbf.write(notebook, target)
print(f'Created {target} with {len(cells)} cells.')
