from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'presentations/output/figures'
OUT.mkdir(parents=True, exist_ok=True)
COLORS = {'EE-CON': '#2267A4', 'RE-CON': '#AC5226'}
PLASMA_TIMES = ['during_20_min', 'during_40_min', 'post_10_min', 'post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
PLASMA_LABELS = ['During\n20 min', 'During\n40 min', 'Post\n10 min', 'Post\n30 min', 'Post\n3.5 h', 'Post\n24 h']
RNA_TIMES = ['post_15_30_45_min', 'post_3.5_4_hr', 'post_24_hr']
RNA_LABELS = ['Post\n15 min', 'Post\n3.5 h', 'Post\n24 h']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.labelsize': 14, 'axes.titlesize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'svg.fonttype': 'none'})

def panel(ax, df, times, labels, title, ylim, q_labels):
    ax.set_title(title, loc='left', fontweight='bold', color='#172F40', pad=12)
    ax.axhline(0, color='#77838B', lw=1.2)
    for group, offset in [('EE-CON', -.11), ('RE-CON', .11)]:
        rows = df[df.contrast_category.eq(group)].copy()
        assert len(rows) == len(set(rows.Timepoint))
        assert set(rows.Timepoint).issubset(times)
        x = rows.Timepoint.map({t:i for i,t in enumerate(times)}).to_numpy() + offset
        y = rows.logFC.to_numpy()
        lo, hi = rows['CI.L_calculated'].to_numpy(), rows['CI.R_calculated'].to_numpy()
        assert np.all(lo <= y) and np.all(hi >= y)
        ax.errorbar(x, y, yerr=np.vstack([y-lo, hi-y]), fmt='none', ecolor=COLORS[group], capsize=4, elinewidth=1.7, zorder=2)
        sig = rows.adj_p_value.lt(.05).to_numpy()
        ax.scatter(x, y, s=70, facecolors='white', edgecolors=COLORS[group], linewidths=1.6, zorder=3)
        ax.scatter(x[sig], y[sig], s=70, c=COLORS[group], edgecolors=COLORS[group], zorder=4)
        for xx, (_,row) in zip(x,rows.iterrows()):
            key=(group,row.Timepoint)
            if key in q_labels:
                text, dy = q_labels[key]
                ax.annotate(text, (xx,row['CI.R_calculated']), xytext=(0,dy), textcoords='offset points', ha='center', fontsize=12, color=COLORS[group], fontweight='bold')
    ax.set_xticks(range(len(times)),labels)
    ax.set_xlim(-.48,len(times)-.52)
    ax.set_ylim(ylim)
    ax.set_ylabel('Effect vs resting control\n(source logFC)', color='#172F40')
    ax.grid(axis='y', color='#E1E6EA', linewidth=.8)
    ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#A7B2BA')

def save(gene, plasma, rna, limits, q_labels):
    fig, axes = plt.subplots(1,2,figsize=(13.6,5.0),gridspec_kw={'width_ratios':[1.55,1]})
    fig.subplots_adjust(left=.073,right=.985,bottom=.18,top=.79,wspace=.3)
    panel(axes[0],plasma,PLASMA_TIMES,PLASMA_LABELS,'Plasma protein (Olink)',limits[0],q_labels[0])
    panel(axes[1],rna,RNA_TIMES,RNA_LABELS,'Muscle RNA',limits[1],q_labels[1])
    handles = [Line2D([0],[0],marker='o',ls='none',markerfacecolor='white',markeredgecolor=color,markeredgewidth=1.6,markersize=8,label=label) for label,color in [('Endurance vs control',COLORS['EE-CON']),('Resistance vs control',COLORS['RE-CON'])]]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False,fontsize=14)
    for ext in ['png','svg']:
        fig.savefig(OUT/f'{gene}_plasma_muscle_timecourses.{ext}',dpi=220,facecolor='white')
    plt.close(fig)

cxp=pd.read_csv(ROOT/'results/fractalkine/cx3cl1_plasma_differential.csv')
cxr=pd.read_csv(ROOT/'results/fractalkine/cx3cl1_tissue_rna_differential.csv').query("tissue == 'muscle'")
w=pd.read_csv(ROOT/'results/wars1/wars1_primary.csv')
wp=w.query("tissue == 'blood' and assay == 'prot-ol'")
wr=w.query("tissue == 'muscle' and assay == 'transcript-rna-seq'")
assert len(cxp)==len(wp)==10 and len(cxr)==len(wr)==6
assert cxp.adj_p_value.lt(.05).sum()==2 and cxr.adj_p_value.lt(.05).sum()==3
assert wp.adj_p_value.lt(.05).sum()==1 and wr.adj_p_value.lt(.05).sum()==2
save('cx3cl1',cxp,cxr,[(-.38,.98),(-.7,4.6)],[{('EE-CON','during_20_min'):('q =\n0.000364',10),('EE-CON','during_40_min'):('q =\n0.000227',10)},{}])
save('wars1',wp,wr,[(-.55,1.10),(-.22,.65)],[{('RE-CON','post_10_min'):('q = 0.0312',10)},{}])
print('Created two figures from 32 source estimates. q denotes original BH-adjusted p.')
