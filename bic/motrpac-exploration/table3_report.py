"""Create a compact report and coverage figure from the systematic screen."""
from pathlib import Path
import hashlib,json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

base=Path(__file__).parent;p=base/'table3'
s=pd.read_csv(p/'coverage_summary.csv').fillna('')
r=pd.read_csv(p/'all_candidate_results.csv')
b=pd.read_csv(p/'plasma_results.csv')
keys=['candidate_id','source_object','feature_id','contrast_category','Timepoint']
assert not r.duplicated(keys).any(), 'Duplicate hypothesis rows'
assert len(s)==28 and s.candidate_id.nunique()==28
assert ((r['CI.L_calculated']<=r.logFC)&(r.logFC<=r['CI.R_calculated'])).all()
source_hits=b[b.identity_resolved & (b.adj_p_value<.05)]
family_hits=b[b.identity_resolved & (b.table3_screen_bh<.05)]
counts=s.plasma_status.value_counts().to_dict()
family_extra=sorted(set(family_hits.candidate_name)-set(source_hits.candidate_name))
labels={'response_detected':'Detected response',
 'present_no_source_FDR_hit':'Tested; no source-FDR hit',
 'absent_from_analyzed_tables':'Absent from analyzed plasma tables',
 'identity_unresolved':'Chemical identity unresolved'}
md=['# Systematic MoTrPAC screen of Table 3','',
'Source: Chow et al., *Exerkines in health, resilience and disease* (2022), Table 3, printed pages 279–280 (PDF pages 7–8). The supplied PDF was extracted and both complete table pages visually checked. All 28 entries were retained, including absent and unresolved matches.','',
'This is a systematic candidate lookup in the public human MoTrPAC acute-exercise results, not a systematic review of all literature or independent experimental replication.','',
'## Findings','',
f'- {sum(s.plasma_status.isin(["response_detected","present_no_source_FDR_hit"]))} of 28 entries have an identity-resolved match in the analyzed circulating protein/metabolite tables.',
f'- {counts.get("response_detected",0)} entries show at least one exercise-versus-control change at the original MoTrPAC adjusted p<0.05: '+', '.join(sorted(source_hits.candidate_name.unique()))+'.',
f'- {counts.get("present_no_source_FDR_hit",0)} entries are tested in plasma but have no source-FDR hit; this is not evidence of equivalence or absence of biological response.',
f'- {counts.get("absent_from_analyzed_tables",0)} entries have no match in the analyzed plasma tables; that can reflect assay coverage or analysis/QC eligibility.',
f'- {counts.get("identity_unresolved",0)} entry (BAIBA) is retained as an unresolved chemical annotation.',
'',
'## Every Table 3 entry','',
'Tissue columns summarize any source-adjusted p<0.05 across EE-CON or RE-CON: up, down, or both. “NS” means tested with no detected change; “—” means no matching analyzed feature. Tissue results are context and do not establish secretion.','',
'| Table 3 entry | Plasma status | Muscle RNA | Adipose RNA | Muscle protein | Adipose protein |',
'|---|---|---|---|---|---|']
def fmt(v):return {'absent':'—','tested_NS':'NS','identity_unresolved':'Identity?'}.get(v,str(v).replace('down;up','up + down'))
for _,z in s.iterrows():
 md.append('| '+' | '.join([z['name'],labels[z.plasma_status]]+[fmt(z[k]) for k in ['muscle|transcript-rna-seq','adipose|transcript-rna-seq','muscle|prot-pr','adipose|prot-pr']])+' |')
md += ['', '## Circulating responses passing the original MoTrPAC FDR','',
'These are existing consortium model estimates for exercise-related changes relative to time-matched non-exercise controls. Effects are reported in each source assay’s model scale; they are not assumed to be directly comparable between assay platforms.','',
'| Candidate | Assay/platform | Contrast | Time | Effect | Source adjusted p | Table-3-family BH |',
'|---|---|---|---|---:|---:|---:|']
for _,z in source_hits.sort_values(['candidate_name','source_object','contrast_category','Timepoint']).iterrows():
 md.append(f'| {z.candidate_name} | {z.platform} | {z.contrast_category} | {z.Timepoint} | {z.logFC:.4g} | {z.adj_p_value:.3g} | {z.table3_screen_bh:.3g} |')
md += ['', '## Multiple-testing definitions','',
'The original MoTrPAC `adj_p_value` adjusts within tissue/assay/platform/contrast. It is preserved and used for the main coverage figure and response classifications.','',
f'An additional exploratory BH correction covers all {int(b.identity_resolved.sum())} identity-resolved circulating candidate-by-assay-by-time-by-mode tests, using their raw p-values. It includes negative results, all measured Olink assay IDs, and both research and clinical lactate assays. BAIBA is excluded pending chemical-identity resolution. The smaller candidate family has a different null family from the source correction; passing it does not mean passing the source correction.','',
'Candidates with hits only under this additional Table-3-family correction: '+(', '.join(family_extra) if family_extra else 'none')+'. These should be labeled exploratory and secondary.','',
'The tissue-context screen does not receive a new cross-tissue/cross-time FDR correction here. Its source-FDR hits should not be treated as a newly controlled omnibus discovery list.','',
'## Identity and interpretation checks','',
'- **BAIBA:** `METABOLOMICS_CVS` includes a raw `beta-Aminoisobutyric-acid` assay label mapped to `Aminoisobutyric acid`. The public gene/feature map associates that name with KEGG C03665, which identifies 2-/alpha-aminoisobutyric acid. The selected blood assay uses the generic name; the muscle raw label uses beta. Those inconsistent annotations require resolution before making beta-BAIBA claims. [KEGG C03665](https://www.genome.jp/dbget-bin/www_bget?cpd%3AC03665=).',
'- **12,13-diHOME:** no exact analyzed match was found. 9,10-diHOME, which appears in annotation metadata, is a different regioisomer and was excluded.',
'- **Irisin/FNDC5:** FNDC5 transcript or tissue protein is a precursor measurement, not proof of mature circulating irisin.',
'- **Myonectin:** mapped to ERFE, retaining aliases CTRP15/FAM132B. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/151176). Musclin maps to OSTN. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/344901).',
'- **HSP72:** HSPA1A/HSPA1B are grouped as candidate paralogs, but unique physical assay features are tested once per candidate/contrast.',
'- **VEGF:** operationally mapped to VEGFA; other VEGF-family genes are not silently substituted.',
'- **Review arrows:** the table footnote defines them as changes in plasma levels. The chronic-training column is retained as metadata but is not tested with this acute-bout dataset.',
'- **Species and effect:** H/A/C encode human/animal/cell evidence; A/P/E in the source effect column mean autocrine/paracrine/endocrine. A row listing H/A/C does not imply that every biological action was demonstrated in humans.',
'- **Timing:** `post_15_30_45_min` is 15 min for muscle, 30 min for blood, 45 min for adipose. `post_3.5_4_hr` is 3.5 h for muscle/blood and 4 h for adipose. During-exercise RE sampling is unavailable.',
'- **Nonsignificance:** not a formal equivalence test, and it does not contradict a heterogeneous review claim across different exercise doses, populations, species, sampling times, or assay sensitivity.',
'- **Causality:** RNA, tissue abundance, or phosphosite changes do not establish release into blood, source tissue, enzyme activity, target-organ response, or clinical benefit.','',
'## Scope and reproducibility','',
'Public human pre-suspension acute-exercise aggregate release, pinned to commit `535b4044e7417413de471104c619120337602b77`. Thirteen differential-result objects were searched: RNA in blood, muscle and adipose; plasma Olink protein; muscle/adipose protein and phosphoprotein; metabolomics in plasma, muscle and adipose; and the available blood clinical chemistry/protein tables. Clinical protein tables added no Table 3 matches. Epigenomics is outside this abundance/response screen.','',
f'Full output contains {len(r):,} candidate/assay/contrast rows, including the direct EE-RE comparisons. The primary summaries use only exercise-with-controls contrasts. No individual-level data, source assignment, disease-cohort comparison, or new differential model was inferred.','',
'Files: `table3_candidates.csv` (28-row source manifest), `coverage_summary.csv`, `all_candidate_results.csv`, `plasma_results.csv`, `assay_inventory.csv`, `coverage_matrix.png`, and `input_provenance.json`. Source scripts are one directory above: `table3_manifest.py`, `table3_screen.R`, and `table3_report.py`.','',
'```sh','python3 table3_manifest.py','Rscript table3_screen.R .','python3 table3_report.py','```','']
(p/'REPORT.md').write_text('\n'.join(md))

# Coverage figure: source-FDR detections, not a quantitative effect heatmap.
cols=['blood|prot-ol','blood|metab','blood|transcript-rna-seq',
 'muscle|transcript-rna-seq','muscle|prot-pr','muscle|prot-ph','muscle|metab',
 'adipose|transcript-rna-seq','adipose|prot-pr','adipose|prot-ph','adipose|metab']
clabels=['Plasma\nprotein','Plasma\nmetabolite','Blood\nRNA','Muscle\nRNA','Muscle\nprotein',
 'Muscle\nphosphosite','Muscle\nmetabolite','Adipose\nRNA','Adipose\nprotein','Adipose\nphosphosite','Adipose\nmetabolite']
codes={'absent':0,'tested_NS':1,'up':2,'down':3,'down;up':4,'up;down':4,'identity_unresolved':5}
colors=['#e5e7eb','#ffffff','#df766c','#76a4d0','#a38bbb','#e8c66d']
symbols={0:'—',1:'·',2:'↑',3:'↓',4:'↕',5:'?'}
v=np.array([[codes[z[c]] for c in cols] for _,z in s.iterrows()])
fig,ax=plt.subplots(figsize=(14.7,13.2))
ax.imshow(v,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(-.5,6.5),6),aspect='auto')
for i in range(v.shape[0]):
 for j in range(v.shape[1]):ax.text(j,i,symbols[v[i,j]],ha='center',va='center',fontsize=12,color='#333')
ax.set_xticks(range(len(cols)),clabels,fontsize=9);ax.xaxis.tick_top()
ax.set_yticks(range(28),s['name'],fontsize=10)
ax.set_xticks(np.arange(-.5,len(cols),1),minor=True);ax.set_yticks(np.arange(-.5,28,1),minor=True)
ax.grid(which='minor',color='#ccd1d7',linewidth=.5);ax.tick_params(which='both',length=0)
ax.axvline(1.5,color='#505761',linewidth=2)
for side in ax.spines.values():side.set_visible(False)
fig.suptitle('Table 3: what is observable in human MoTrPAC?',x=.035,ha='left',fontsize=19,fontweight='bold',y=.988)
fig.text(.035,.952,'All 28 cardiometabolic exerkine entries · acute exercise versus time-matched controls',fontsize=11,color='#555')
patches=[Patch(facecolor=colors[i],edgecolor='#888',label=l) for i,l in enumerate(['No analyzed match','Tested; no FDR hit','Increase','Decrease','Both directions','Identity unresolved'])]
fig.legend(handles=patches,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.54,.041),fontsize=10)
fig.text(.035,.018,'Source BH-adjusted p < 0.05 at any sampled time/mode. Tissue signals do not establish circulating secretion.\nPlasma metabolite column includes research and clinical lactate assays. “No FDR hit” does not establish no response.',fontsize=9,color='#555')
fig.tight_layout(rect=(.02,.10,.995,.923))
fig.savefig(p/'coverage_matrix.png',dpi=170,facecolor='white')
fig.savefig(p/'coverage_matrix.svg',facecolor='white')

files=list(base.glob('*.rda'))
review_pdf=base.parent/'references'/'s41574-022-00641-2.pdf'
provenance={'motrpac_commit':'535b4044e7417413de471104c619120337602b77',
 'review_pdf':str(review_pdf),
 'review_pdf_sha256':hashlib.sha256(review_pdf.read_bytes()).hexdigest(),
 'input_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
 'counts':counts,'circulating_family_tests':int(b.identity_resolved.sum()),
 'source_fdr_candidates':sorted(source_hits.candidate_name.unique()),
 'secondary_family_only_candidates':family_extra}
(p/'input_provenance.json').write_text(json.dumps(provenance,indent=2))
print(json.dumps({k:v for k,v in provenance.items() if k not in ['input_sha256','review_pdf_sha256']},indent=2))
