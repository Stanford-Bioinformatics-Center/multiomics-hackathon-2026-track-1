from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'presentations/output/figures'
COLORS={'EE-CON':'#2267A4','RE-CON':'#AC5226'}
PT=['during_20_min','during_40_min','post_10_min','post_15_30_45_min','post_3.5_4_hr','post_24_hr']
PL=['During\n20 min','During\n40 min','Post\n10 min','Post\n30 min','Post\n3.5 h','Post\n24 h']
RT=PT[3:]
RL=['Post\n15 min','Post\n3.5 h','Post\n24 h']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.labelsize':14,'axes.titlesize':16,'xtick.labelsize':12,'ytick.labelsize':12,'svg.fonttype':'none'})

def draw(ax,df,times,labels,title,ylim):
    ax.set_title(title,loc='left',weight='bold',color='#172F40',pad=12)
    ax.axhline(0,color='#77838B',lw=1.2)
    for group,offset in [('EE-CON',-.11),('RE-CON',.11)]:
        d=df[df.contrast_category.eq(group)]
        x=d.Timepoint.map({t:i for i,t in enumerate(times)}).to_numpy()+offset
        y=d.logFC.to_numpy()
        lo=d['CI.L_calculated'].to_numpy();hi=d['CI.R_calculated'].to_numpy()
        assert np.all(lo<=y) and np.all(hi>=y)
        ax.errorbar(x,y,yerr=np.vstack([y-lo,hi-y]),fmt='none',ecolor=COLORS[group],capsize=4,elinewidth=1.7,zorder=2)
        hit=d.adj_p_value.lt(.05).to_numpy()
        ax.scatter(x,y,s=70,facecolors='white',edgecolors=COLORS[group],linewidths=1.6,zorder=3)
        ax.scatter(x[hit],y[hit],s=70,c=COLORS[group],zorder=4)
    ax.set_xticks(range(len(times)),labels)
    ax.set_xlim(-.48,len(times)-.52);ax.set_ylim(ylim)
    ax.set_ylabel('Effect vs resting control\n(source logFC)',color='#172F40')
    ax.grid(axis='y',color='#E1E6EA',lw=.8);ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#A7B2BA')

def get(df,group,t):
    d=df[df.contrast_category.eq(group)&df.Timepoint.eq(t)]
    assert len(d)==1
    return d.iloc[0]

def label(ax,df,group,t,times,text,dx=0,dy=10,below=False,ha='center'):
    r=get(df,group,t)
    x=times.index(t)+(-.11 if group=='EE-CON' else .11)
    y=r['CI.L_calculated' if below else 'CI.R_calculated']
    ax.annotate(text,(x,y),xytext=(dx,dy),textcoords='offset points',ha=ha,
                va='top' if below else 'bottom',fontsize=11.5,weight='bold',color=COLORS[group])

cxp=pd.read_csv(ROOT/'results/fractalkine/cx3cl1_plasma_differential.csv')
cxr=pd.read_csv(ROOT/'results/fractalkine/cx3cl1_tissue_rna_differential.csv').query("tissue == 'muscle'")
w=pd.read_csv(ROOT/'results/wars1/wars1_primary.csv')
wp=w.query("tissue == 'blood' and assay == 'prot-ol'")
wr=w.query("tissue == 'muscle' and assay == 'transcript-rna-seq'")
assert len(cxp)==len(wp)==10 and len(cxr)==len(wr)==6
assert get(wp,'EE-CON','during_20_min').adj_p_value>.05
assert get(wp,'EE-CON','during_40_min').p_value<.05<get(wp,'EE-CON','during_40_min').adj_p_value
assert cxp[cxp.contrast_category.eq('RE-CON')].adj_p_value.ge(.05).all()
assert get(cxr,'RE-CON',RT[0]).adj_p_value<.05
assert get(wr,'EE-CON',RT[1]).adj_p_value<.05

for gene,plasma,rna,limits in [('cx3cl1',cxp,cxr,[(-.58,1.00),(-.7,4.65)]),('wars1',wp,wr,[(-.57,1.15),(-.22,.76)])]:
    fig,ax=plt.subplots(1,2,figsize=(13.6,5.0),gridspec_kw={'width_ratios':[1.55,1]})
    fig.subplots_adjust(left=.073,right=.985,bottom=.18,top=.79,wspace=.3)
    draw(ax[0],plasma,PT,PL,'Plasma protein (Olink)',limits[0])
    draw(ax[1],rna,RT,RL,'Muscle RNA',limits[1])
    if gene=='cx3cl1':
        for t in PT[:2]:
            q=get(plasma,'EE-CON',t).adj_p_value
            label(ax[0],plasma,'EE-CON',t,PT,f'q =\n{q:.6f}')
        label(ax[0],plasma,'RE-CON',PT[2],PT,f'q = {get(plasma,"RE-CON",PT[2]).adj_p_value:.3f}',dy=-14,below=True)
        ax[0].text(.01,.035,'No blood samples\nduring resistance',transform=ax[0].transAxes,color=COLORS['RE-CON'],fontsize=11.5,va='bottom')
        # Separate labels vertically because the two early RNA estimates share a timepoint.
        label(ax[1],rna,'EE-CON',RT[0],RT,'q = 5.50e−53',dx=-6,dy=10)
        label(ax[1],rna,'RE-CON',RT[0],RT,'q = 1.55e−33',dx=6,dy=-16,below=True)
    else:
        label(ax[0],plasma,'EE-CON',PT[0],PT,f'q = {get(plasma,"EE-CON",PT[0]).adj_p_value:.3f}')
        label(ax[0],plasma,'EE-CON',PT[1],PT,f'raw p = {get(plasma,"EE-CON",PT[1]).p_value:.3f}\nq = {get(plasma,"EE-CON",PT[1]).adj_p_value:.3f}',dy=10)
        label(ax[0],plasma,'RE-CON',PT[2],PT,f'q = {get(plasma,"RE-CON",PT[2]).adj_p_value:.4f}',dy=10)
        label(ax[1],rna,'EE-CON',RT[1],RT,'q = 1.48e−6',dx=-18,dy=21)
        label(ax[1],rna,'RE-CON',RT[1],RT,'q = 3.63e−6',dx=22,dy=-16,below=True)
    handles=[Line2D([0],[0],marker='o',ls='none',markerfacecolor='white',markeredgecolor=color,markeredgewidth=1.6,markersize=8,label=label) for label,color in [('Endurance vs control',COLORS['EE-CON']),('Resistance vs control',COLORS['RE-CON'])]]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False,fontsize=14)
    for ext in ('png','svg'):
        fig.savefig(OUT/f'{gene}_both_modes.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
print('Plotted all 32 plasma/muscle RNA estimates with explicit labels for both exercise modes.')
