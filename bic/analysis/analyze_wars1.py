"""Assemble WARS1 source statistics and published disease evidence without refitting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ['HUMAN_FEATURE_TO_GENE','BLOOD_PROT_OL_DA','BLOOD_TRNSCRPT_DA',
          'MUSCLE_TRNSCRPT_DA','ADIPOSE_TRNSCRPT_DA','MUSCLE_PROT_PR_DA',
          'ADIPOSE_PROT_PR_DA','MUSCLE_PROT_PH_DA','ADIPOSE_PROT_PH_DA',
          'BLOOD_METAB_DA','MUSCLE_METAB_DA','ADIPOSE_METAB_DA','METABOLOMICS_CVS']
SOURCES = {
 'secretion': 'https://pubmed.ncbi.nlm.nih.gov/36640342/',
 'immune': 'https://doi.org/10.3390/biom10091283',
 'mini_wars': 'https://doi.org/10.1038/s41467-022-31904-1',
 't2_wars': 'https://doi.org/10.1074/jbc.C400431200',
 'insulin': 'https://doi.org/10.1007/s00018-023-05082-2',
 'hypertension': 'https://doi.org/10.1038/s41588-025-02096-3',
 'ogtt_genetics': 'https://doi.org/10.1007/s00125-026-06800-8',
 'ogtt_tables': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs00125-026-06800-8/MediaObjects/125_2026_6800_MOESM2_ESM.xlsx',
 'prior_exercise_transcript': 'https://genescells.ru/2313-1829/issue/view/6122',
 'olink': 'https://olink.com/assay/explore/neurology/tryptophan-trna-ligase-cytoplasmic',
 'uniprot': 'https://www.uniprot.org/uniprotkb/P23381/entry',
 'motrpac': 'https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77',
}

def checksum(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify_inputs(root=ROOT):
 source = root / 'motrpac-exploration'
 expected = json.loads((source/'table3/input_provenance.json').read_text())['input_sha256']
 hashes = {n+'.rda': checksum(source/(n+'.rda')) for n in INPUTS}
 assert all(h == expected[n] for n,h in hashes.items())
 return hashes

def add_time(frame):
 out=frame.copy()
 simple={'during_20_min':'During 20 min','during_40_min':'During 40 min',
         'post_10_min':'Post 10 min','post_24_hr':'Post 24 h','pre_exercise':'Baseline'}
 out['sample_time']=out.Timepoint.map(simple)
 for raw,labels in [('post_15_30_45_min',{'blood':'Post 30 min','muscle':'Post 15 min','adipose':'Post 45 min'}),
                    ('post_3.5_4_hr',{'blood':'Post 3.5 h','muscle':'Post 3.5 h','adipose':'Post 4 h'})]:
  ix=out.Timepoint.eq(raw);out.loc[ix,'sample_time']=out.loc[ix,'tissue'].map(labels)
 assert out.sample_time.notna().all()
 out['source_significant']=out.adj_p_value < .05
 return out

def build(root=ROOT):
 root=Path(root);out=root/'results/wars1';ref=out/'reference_data'
 hashes=verify_inputs(root)
 data=add_time(pd.read_csv(out/'wars1_all_source_contrasts.csv'))
 audit=pd.read_csv(out/'bh_audit.csv')
 np.testing.assert_allclose(data.bh_recomputed,data.adj_p_value,rtol=1e-10,atol=1e-12)
 assert audit.max_BH_difference.lt(1e-10).all()
 key=['tissue','assay','platform','feature_id','contrast']
 assert not data.duplicated(key).any()
 # Baseline contrasts also use labels EE-CON / RE-CON. Filter on contrast_type!
 primary=data[data.contrast_type.eq('exercise_with_controls')].copy()
 direct=data[data.contrast_type.eq('Endur_vs_Resist')].copy()
 baseline=data[data.contrast_type.eq('baseline')].copy()
 assert not primary.Timepoint.eq('pre_exercise').any()
 plasma=primary[primary.assay.eq('prot-ol')].copy()
 assert len(plasma)==10 and set(plasma.feature_id)=={'OID21084'}
 previous=pd.read_csv(root/'results/exerkine_screen/focal_genes_all_source_contrasts.csv')
 previous=previous[previous.gene_symbol.eq('WARS1')]
 checked=data.merge(previous,on=key,suffixes=('_new','_prior'),validate='one_to_one')
 assert len(checked)==58
 for column in ['logFC','p_value','adj_p_value']:
  np.testing.assert_allclose(checked[column+'_new'],checked[column+'_prior'],rtol=1e-10,atol=1e-12)
 bg=pd.read_csv(root/'results/exerkine_screen/plasma_background_validated.csv')
 plasma=plasma.merge(bg[['feature_id','contrast','pooled_primary_plasma_bh']],on=['feature_id','contrast'],validate='one_to_one')
 arms=data[data.assay.eq('prot-ol') & data.Timepoint.eq('post_10_min') &
           data.contrast_category.isin(['RE-RE','CON-CON','RE-CON']) & ~data.contrast_type.eq('baseline')].copy()
 assert len(arms)==3
 effects=arms.set_index('contrast_category').logFC
 np.testing.assert_allclose(effects['RE-RE']-effects['CON-CON'],effects['RE-CON'],atol=1e-12)
 for name,frame in [('wars1_primary',primary),('wars1_plasma_primary',plasma),('wars1_direct_EE_RE',direct),
                    ('wars1_baseline_comparisons',baseline),('wars1_post10_plasma_decomposition',arms)]:
  frame.to_csv(out/(name+'.csv'),index=False)
 m=add_time(pd.read_csv(out/'tryptophan_kynurenine_source_context.csv'))
 m['correction_status']='Original source q; attempted independent BH groupings did not reproduce it'
 m.to_csv(out/'metabolite_context_annotated.csv',index=False)
 # Original published workbook is read only. Preserve its columns and allele directions.
 workbook=ref/'uluvar_2026_esm_tables.xlsx'
 diseases=pd.read_excel(workbook,sheet_name='ESM Table 16',header=3)
 diseases=diseases[diseases.prot.isin(['WARS','WARS1'])].copy()
 tissues=pd.read_excel(workbook,sheet_name='ESM Table 17',header=4)
 tissues=tissues[tissues.tissue.notna()].copy()
 assert len(diseases)==2 and len(tissues)==49
 assert tissues.Allele_ref.eq('C').all() and tissues.Allele_alt.eq('G').all()
 assert tissues['Effect.aligned'].lt(0).all()
 # The table gives effects for G. Flip ONLY this single-variant tissue table
 # to the WARS-increasing C allele, as in the authors' text and figure.
 tissues['effect_allele_for_display']='C'
 tissues['eqtl_beta_C']=-tissues.slope
 tissues['eqtl_CI_low_C']=tissues.eqtl_beta_C-1.96*tissues.slope_se
 tissues['eqtl_CI_high_C']=tissues.eqtl_beta_C+1.96*tissues.slope_se
 tissues['plasma_pqtl_beta_C']=-tissues['Effect.aligned']
 # ESM16 uses different lead SNPs for protein and outcome. Do not harmonize
 # disease effect direction by sign-flipping unrelated lead-SNP coefficients.
 diseases['interpretation']='Reported colocalization; outcome beta belongs to topsnp_rsid_other, not necessarily rs2273804'
 diseases.to_csv(out/'published_wars1_disease_colocalization.csv',index=False)
 tissues.to_csv(out/'published_wars1_tissue_eqtl.csv',index=False)
 raw15=pd.read_excel(workbook,sheet_name='ESM Table 15',header=None)
 row15=raw15[raw15[0].eq('OID21084')]
 assert len(row15)==1
 r=row15.iloc[0]
 ogtt=pd.DataFrame([{'olink_id':r[0],'gene':r[1],'interaction_p':float(r[16]),
                     'interaction_q':float(r[17]),'authors_interaction_threshold':.20,
                     'passes_005':bool(float(r[17])<.05),'n_initial_participants':11,
                     'source_sheet':'ESM Table 15','source_excel_row':int(row15.index[0])+1}])
 assert np.isclose(ogtt.interaction_q.iloc[0],.086227,atol=1e-6)
 ogtt.to_csv(out/'published_ogtt_wars1_interaction.csv',index=False)
 counts={'source_WARS1_contrasts':len(data),'primary_WARS1_contrasts':len(primary),
         'primary_BH_hits':int(primary.source_significant.sum()),'BH_families_verified':len(audit),
         'background_feature_tests_verified':int(audit.n_features.sum()),
         'maximum_BH_difference':float(audit.max_BH_difference.max()),'prior_screen_rows_reproduced':len(checked)}
 versions={'python':__import__('sys').version.split()[0],'numpy':np.__version__,'pandas':pd.__version__}
 provenance={'date':'2026-09-26','dataset':'human-precovid-sed-adu c2.0','package_version':'2.0.8',
             'commit':'535b4044e7417413de471104c619120337602b77','input_sha256':hashes,
             'reference_sha256':{p.name:checksum(p) for p in sorted(ref.iterdir()) if p.is_file()},
             'source_urls':SOURCES,'versions':versions,'counts':counts,
             'model_fitting':'Source p-values, estimates and CIs; no participant-level refit.',
             'primary_definition':'contrast_type==exercise_with_controls, not contrast_category alone',
             'metabolite_limitation':'Source q-values retained; correction not independently reproduced with the two tested groupings. See saved diagnostics.',
             'genetic_analysis':'Extracted author-reported results; no new GWAS, MR or colocalization. ESM17 G-allele slopes flipped to C; ESM16 original lead-SNP effects preserved.',
             'selection_limitation':'WARS1 was selected in the preceding assay-wide screen; this is not an independent validation cohort.'}
 (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 return dict(data=data,primary=primary,direct=direct,baseline=baseline,plasma=plasma,arms=arms,
             metabolites=m,diseases=diseases,tissue_eqtl=tissues,ogtt=ogtt,audit=audit,counts=counts)

if __name__=='__main__':
 result=build();print(json.dumps(result['counts'],indent=2))
