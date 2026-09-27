import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';

const ROOT='/Users/acheron/multiomics-hackathon-2026-track-1/bic/presentations';
const TMP=path.join(ROOT,'.build_diabetes_extension');
const OUT=path.join(ROOT,'output');
const SKILL='/Users/acheron/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME='/Users/acheron/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=path.join(RUNTIME,'node/node_modules');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const sourcePath=path.join(TMP,'source.pptx');
const sourceHash=(await fs.readFile(path.join(TMP,'source.sha256'),'utf8')).trim();
const activeSource=path.join(OUT,'fractalkine_wars1_with_both_modes.pptx');
if(createHash('sha256').update(await fs.readFile(activeSource)).digest('hex')!==sourceHash) throw new Error('Source changed after inspection. Reconcile the latest user edits before exporting.');
const presentation=await PresentationFile.importPptx(await FileBlob.load(sourcePath));
if(presentation.slides.items.length!==3) throw new Error('Unexpected source slide count.');
const font='Helvetica Neue';
const C={ink:'#172F40',blue:'#2267A4',orange:'#AC5226',muted:'#526571',paper:'#FAFBFC',white:'#FFFFFF',light:'#EEF3F6'};
function txt(s,text,x,y,w,h,size=25,color=C.ink,bold=false){
  const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  sh.text=text;sh.text.style={typeface:font,fontSize:size,color,bold,autoFit:'none'};return sh;
}
function base(title){
  const s=presentation.slides.add();s.background.fill=C.paper;
  txt(s,title,60,42,1160,100,44,C.ink,true);return s;
}
const refs={
  islet:'[1] Rutti S et al. Fractalkine (CX3CL1), a new factor protecting beta-cells against TNF-alpha. Molecular Metabolism 2014;3:731–741. https://doi.org/10.1016/j.molmet.2014.07.007 ; https://pubmed.ncbi.nlm.nih.gov/25353001/',
  secretion:'[2] Nguyen TTT et al. Tryptophan-dependent and -independent secretions of tryptophanyl-tRNA synthetase mediate innate inflammatory responses. Cell Reports 2023;42:111905. https://doi.org/10.1016/j.celrep.2022.111905 ; https://pubmed.ncbi.nlm.nih.gov/36640342/',
  insulin:'[3] Sun WX et al. Tryptophanylation of insulin receptor by WARS attenuates insulin signaling. Cellular and Molecular Life Sciences 2024. https://doi.org/10.1007/s00018-023-05082-2 ; https://pubmed.ncbi.nlm.nih.gov/38212570/',
  vascular:'[4] Gioelli N et al. Neuropilin 1 and its inhibitory ligand mini-tryptophanyl-tRNA synthetase inversely regulate VE-cadherin turnover and vascular permeability. Nature Communications 2022. https://doi.org/10.1038/s41467-022-31904-1',
  exerkine:'[5] Chow LS et al. Exerkines in health, resilience and disease. Nature Reviews Endocrinology 2022;18:273–289. https://doi.org/10.1038/s41574-022-00641-2',
};
{
 const s=base('Biological significance and diabetes relevance');
 txt(s,'External biological evidence, separate from the MoTrPAC associations',60,143,1160,42,24,C.muted);
 txt(s,'Fractalkine / CX3CL1',60,201,540,42,30,C.blue,true);
 txt(s,'WARS1',674,201,546,42,30,C.orange,true);
 txt(s,'CX3CL1 signals through CX3CR1 in immune and pancreatic islet biology. [1]',60,255,540,85,25);
 txt(s,'In human islets, CX3CL1 reduced β-cell death and glucagon release, without increasing insulin secretion in that study. [1]',60,353,540,135,25);
 txt(s,'This motivates testing islet survival and glucagon control after exercise.',60,501,540,75,25);
 txt(s,'WARS1 normally charges tRNA. Secreted WARS1 can activate innate immune responses. [2]',674,255,546,85,25);
 txt(s,'With excess tryptophan, intracellular WARS1 modifies insulin receptor K1209 and weakens insulin signaling in experimental cells. [3]',674,353,546,135,25);
 txt(s,'Mini-WARS, a shorter form, can stabilize endothelial junctions. [4]',674,501,546,75,25);
 txt(s,'Protein form and compartment matter. These results do not establish a diabetes benefit from the exercise-associated plasma increases.',60,595,1160,65,24,C.ink,true);
 txt(s,'[1] Rutti 2014   [2] Nguyen 2023   [3] Sun 2024   [4] Gioelli 2022. Full references in notes.',60,675,1160,30,17.5,C.muted);
 s.speakerNotes.textFrame.setText([
  'Purpose: relate the exercise-responsive candidates to independent biological evidence without interpreting the MoTrPAC response as a disease outcome.',
  refs.islet,refs.secretion,refs.insulin,refs.vascular,
  'Evidence details and limits: Rutti et al. exposed human islets to CX3CL1 for 24 hours at 1–50 ng/mL. They reported lower basal human beta-cell apoptosis and lower low-glucose-stimulated glucagon secretion, without increased glucose-stimulated insulin secretion in that study. The TNF-alpha rescue experiments on secretion largely used rat beta cells. These dose and duration conditions cannot be equated with the acute Olink response, which does not establish absolute concentration. Other work reported different insulin-secretory effects (Lee et al., Cell 2013, https://pmc.ncbi.nlm.nih.gov/articles/PMC3717389/), so the slide specifically qualifies the finding as belonging to the Rutti study.',
  'Nguyen et al. demonstrated direct tryptophan-dependent WARS1 secretion and tryptophan-independent release through plasma-membrane-derived vesicles, including extracellular innate immune effects. These were not exercise experiments. WARS1 is the cytoplasmic tryptophanyl-tRNA synthetase, distinct from mitochondrial WARS2.',
  'Sun et al. studied WARS-dependent tryptophanylation of insulin receptor K1209 under excess tryptophan. This reduced insulin-stimulated signaling through IR, AKT and AS160 and reduced glucose uptake. SIRT1 reversed the modification. Substrate availability and the intracellular compartment are essential to this model. The exercise plasma rise does not demonstrate insulin-receptor modification or insulin resistance.',
  'Gioelli et al. studied the mini-WARS splice form, which inhibits NRP1-dependent VE-cadherin turnover and can reduce endothelial permeability. This supports resolving WARS1 forms before assigning vascular function. The Olink measurement in our analysis does not identify full-length, spliced, cleaved or vesicular species.',
  'MoTrPAC results are pinned human c2.0 analysis package 2.0.8, commit 535b4044e7417413de471104c619120337602b77. No glucose-control, islet, endothelial or immune-function endpoint was tested in our differential analysis.'
 ].join('\n\n'));
}
{
 const s=base('Testable hypotheses and WARS1 exerkine status');
 txt(s,'Working explanations for the observed timing, with tests that could challenge them',60,143,1160,42,24,C.muted);
 const values=[
  ['Observed pattern','Working hypothesis','Discriminating test'],
  ['CX3CL1\nPlasma rise during exercise',
   'Rapid shedding or release of pre-existing CX3CL1 could generate the early plasma pulse.',
   'Measure soluble forms across time. Test islet survival and glucagon responses to a calibrated pulse, with and without CX3CR1 blockade.'],
  ['WARS1\nPlasma rise at 10 min, muscle RNA at 3.5 h',
   'Release of existing protein could precede RNA induction. Later RNA may replenish a pool or reflect a parallel stress response.',
   'Identify plasma forms and the release source. Deplete WARS1, then add back the identified form and test immune or endothelial responses.'],
 ];
 const t=s.tables.add({rows:3,columns:3,left:60,top:203,width:1160,height:306,columnWidths:[230,395,535],values});
 t.styleOptions={headerRow:true,bandedRows:false};
 t.borders.assign({fill:'#D8E1E6',width:.7,style:'solid'});
 t.cells.block({row:0,column:0,rowCount:3,columnCount:3}).assign({fill:C.white,textStyle:{typeface:font,fontSize:22,color:C.ink},margins:{left:14,right:14,top:11,bottom:10},anchor:'center'});
 t.cells.block({row:0,column:0,rowCount:1,columnCount:3}).assign({fill:C.ink,textStyle:{typeface:font,fontSize:23,color:C.white,bold:true}});
 t.cells.block({row:2,column:0,rowCount:1,columnCount:3}).fill=C.light;
 t.rows[0].height=48;t.rows[1].height=118;t.rows[2].height=140;
 t.getCell(1,0).text.style={typeface:font,fontSize:22,color:C.blue,bold:true};
 t.getCell(2,0).text.style={typeface:font,fontSize:22,color:C.orange,bold:true};
 txt(s,'WARS1 is a candidate exerkine',60,534,1160,38,28,C.orange,true);
 txt(s,'Its plasma response plus known release and extracellular signaling support candidacy. [2,5]\nMoTrPAC has not established exercise-triggered release, tissue source or target-cell action.',60,582,1160,68,24);
 txt(s,'Alternatives include altered clearance, plasma-volume changes or cell injury. Proposed tests require additional measurements.',60,669,1160,35,18.5,C.muted);
 s.speakerNotes.textFrame.setText([
  'These hypotheses and experiments are proposed, not performed. They do not establish muscle origin, exercise-mode specificity or a protective effect on diabetes.',
  'Observed evidence: CX3CL1 plasma EE-CON at during 20 and 40 min has source BH-adjusted p=0.000364 and 0.000227. WARS1 plasma RE-CON at post 10 min has effect +0.566828, pointwise 95% CI [0.241716,0.891941], adjusted p=0.031163. Muscle WARS1 RNA at post 3.5 h rises in both arms. Our primary source is MoTrPAC c2.0, analysis package 2.0.8, https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77 .',
  'Hypothesis 1: an early plasma pulse is compatible with release or cleavage of pre-existing CX3CL1. It cannot be attributed to later measured muscle RNA. Other sources or altered clearance remain possible. An orthogonal quantitative assay, sampling during exercise and soluble-versus-cell-associated measurements would distinguish mechanisms. For the disease hypothesis, an exercise-calibrated human-islet exposure plus CX3CR1 blockade could test a receptor-dependent effect on survival and glucagon. This requires new experiments and measured absolute concentrations.',
  'Hypothesis 2: stored WARS1 may be released rapidly, while the later muscle RNA response could reflect replenishment or an independent stress/translation program. Current summary contrasts cannot measure protein flux or establish that muscle produced the plasma species. Muscle total protein has no primary BH hit. Resolve the circulating form with orthogonal assays and obtain direct cell/tissue release evidence before assigning a myokine source.',
  'Functional test: compare matched pre/post-exercise plasma under sham depletion, WARS1 depletion and add-back of a verified molecular form at an appropriate measured concentration. Endothelial permeability or innate immune readouts are motivated by external studies. A response that persists after effective specific depletion would weaken a WARS1-dependent mechanism. Include plasma-volume, cell-injury, assay-specificity and endotoxin controls.',
  'Available-data check: use the existing tissue/timepoint summary statistics for pathway-level tests of translation/stress and immune programs. This can support prioritization but cannot establish secretion, between-person coupling, receptor activation or disease mediation. Participant-level correlations require data that were not used here.',
  'Exerkine assessment: following Chow et al., an exerkine is a factor released in response to exercise that acts through endocrine, paracrine or autocrine signaling. The exercise-associated circulating WARS1 signal and experimentally established release/signaling outside exercise make it a candidate. This dataset alone does not complete the release-plus-function evidence chain. Its source could be muscle, vascular or immune cells, or another tissue. An exerkine need not be a myokine and its effects need not be uniformly beneficial.',
  refs.secretion,refs.vascular,refs.exerkine,refs.islet,
  'Neither gene has a significant direct post-exercise endurance-versus-resistance plasma comparison in this analysis. During-resistance samples are unavailable. Nonsignificance does not establish absence or equivalence.'
 ].join('\n\n'));
}
const candidate=path.join(TMP,'candidate.pptx');
await (await PresentationFile.exportPptx(presentation)).save(candidate);
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,
 finalPath:path.join(OUT,'fractalkine_wars1_with_both_modes_biology_hypotheses.pptx'),
 explicitTotalSlideCount:5,
 pythonExecutable:path.join(RUNTIME,'python/bin/python3'),
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','5'],
 requiredNativeTableOwnerSlides:[5],requiredNativeChartOwnerSlides:[],
 fontPolicy:{basis:'reference',families:[font],referencePath:sourcePath,referenceSha256:sourceHash},
 verifyArtifactToolImport:true,receiptPath:path.join(TMP,'validation.json')
});
console.log(JSON.stringify({finalPath:result.finalPath,slideCount:presentation.slides.items.length}));
