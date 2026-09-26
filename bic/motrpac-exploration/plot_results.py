from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

p = Path(__file__).parent
d = pd.read_csv(p / 'candidate_evidence.csv')
d = d[(d.tissue == 'blood') & (d.contrast_type == 'exercise_with_controls')]
times = ['during_20_min', 'during_40_min', 'post_10_min',
         'post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
labels = ['During\n20 min', 'During\n40 min', 'Post\n10 min',
          'Post\n30 min', 'Post\n3.5 h', 'Post\n24 h']
panels = [('CX3CL1', 'Fractalkine / CX3CL1 · plasma protein', 'prot-ol'),
          ('Kynurenic acid', 'Kynurenic acid · plasma metabolite', 'metab'),
          ('Kynurenine', 'Kynurenine · plasma metabolite', 'metab'),
          ('Leucine', 'Leucine · plasma metabolite', 'metab')]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})
fig, axs = plt.subplots(2, 2, figsize=(13, 8.4))
for ax, (name, title, assay) in zip(axs.flat, panels):
    s = d[(d.assay == assay) & ((d.gene_symbol == name) | (d.feature_id == name))]
    ax.axhline(0, color='#8a9199', linewidth=0.8)
    for group, color, offset, legend in [('EE-CON', '#2867b2', -.12, 'Endurance'),
                                        ('RE-CON', '#c76528', .12, 'Resistance')]:
        q = s[s.contrast_category == group].copy()
        q['x'] = q.Timepoint.map({t: i for i, t in enumerate(times)})
        q = q.sort_values('x')
        x = q.x.to_numpy() + offset
        y = q.logFC.to_numpy()
        err = np.vstack([y-q['CI.L_calculated'], q['CI.R_calculated']-y])
        ax.errorbar(x, y, yerr=err, fmt='none', ecolor=color, capsize=3, alpha=.9)
        ax.scatter(x, y, facecolors='white', edgecolors=color, s=48, label=legend, zorder=3)
        sig = (q.adj_p_value < .05).to_numpy()
        ax.scatter(x[sig], y[sig], color=color, s=48, zorder=4)
    ax.set_title(title, loc='left', fontweight='bold', pad=12)
    ax.set_xticks(range(len(times)), labels)
    ax.set_ylabel('Control-adjusted model effect (logFC)')
    ax.grid(axis='y', color='#e8ebef')
    ax.set_axisbelow(True)
axs[0,0].legend(frameon=False, loc='upper right')
fig.suptitle('Observed human exercise responses: candidates for disease hypotheses',
             x=.06, ha='left', fontsize=16, fontweight='bold', y=.98)
fig.text(.06, .927, 'Public MoTrPAC aggregate results · endurance / resistance vs time-matched controls', color='#555')
fig.text(.06, .035,
         'Points: published model estimates; bars: published 95% CIs. Filled: source BH-adjusted p < 0.05.\n'
         'Time points are equally spaced for display. Resistance has no during-exercise samples. '
         'No disease outcomes are shown.\n'
         'Exploratory selection; source FDR is per contrast, not across this screen. Commit 535b4044e741.',
         fontsize=9, color='#555')
fig.tight_layout(rect=(.015,.11,1,.91), h_pad=2.3, w_pad=2)
fig.savefig(p/'exercise_candidate_trajectories.png', dpi=180, facecolor='white')
fig.savefig(p/'exercise_candidate_trajectories.svg', facecolor='white')
