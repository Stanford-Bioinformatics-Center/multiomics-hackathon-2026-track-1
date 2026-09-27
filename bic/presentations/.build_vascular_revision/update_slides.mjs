import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob, PresentationFile} from '@oai/artifact-tool';

const ROOT='/Users/acheron/multiomics-hackathon-2026-track-1/bic/presentations';
const TMP=path.join(ROOT,'.build_vascular_revision');
const OUT=path.join(ROOT,'output');
const SKILL='/Users/acheron/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME='/Users/acheron/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=path.join(RUNTIME,'node/node_modules');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const sourcePath=path.join(TMP,'source.pptx');
const sourceHash=(await fs.readFile(path.join(TMP,'source.sha256'),'utf8')).trim();
const activeSource=path.join(OUT,'fractalkine_wars1_with_both_modes_biology_hypotheses.pptx');
function hash(b){return createHash('sha256').update(b).digest('hex');}
if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed after inspection. Reconcile before exporting.');
const p=await PresentationFile.importPptx(await FileBlob.load(sourcePath));
if(p.slides.items.length!==5) throw new Error('Unexpected source slide count');
const notes=JSON.parse(await fs.readFile(path.join(TMP,'source-notes.json'),'utf8'));
const set=(id,text)=>{p.resolve(id).text=text;};

set('sh/bip8jmho','WARS1: vascular signaling');
p.resolve('sh/bip8jmho').position={left:674,top:179.87,width:546,height:75.13};
set('sh/r65knqtk','Mini-WARS stabilizes endothelial junctions through NRP1 / VE-cadherin and reduces permeability. [4]');
p.resolve('sh/r65knqtk').position={left:674,top:255,width:546,height:98};
set('sh/q5wjelsz','T2-WARS binds VE-cadherin and inhibits endothelial migration and new blood-vessel formation. [6]');
p.resolve('sh/q5wjelsz').position={left:674,top:365,width:546,height:110};
set('sh/ove9o7yd','Olink does not identify which WARS1 form rises after exercise. Vascular activity remains a testable hypothesis.');
p.resolve('sh/ove9o7yd').position={left:674,top:510,width:546,height:110};
set('sh/mtwrmxg7','[1] Rutti 2014   [4] Gioelli 2022   [6] Tzima 2005. Full references in notes.');

const t=p.resolve('tb/adonexsj');
t.cells.set(2,1,'An extracellular WARS1 form could reduce vascular leakage or restrain vessel sprouting after exercise.');
t.cells.set(2,2,'Identify plasma WARS1 forms. Deplete and add back the detected form. Test endothelial permeability, VE-cadherin junctions and VEGF-driven sprouting.');
set('sh/z2tcnm5s','WARS1: candidate vascular exerkine');
set('sh/yhkbe1o7','The plasma rise and known extracellular activity support a vascular exerkine hypothesis. [2,4–6]\nExercise-triggered release, the active form and vascular effects remain unproven.');

const t2ref='[6] Tzima E et al. VE-cadherin links tRNA synthetase cytokine to anti-angiogenic function. Journal of Biological Chemistry 2005;280:2405–2408. https://doi.org/10.1074/jbc.C400431200 ; https://pubmed.ncbi.nlm.nih.gov/15579907/';
const replaceChecked=(s,oldText,newText)=>{if(!s.includes(oldText)) throw new Error('Missing note text: '+oldText.slice(0,50));return s.replace(oldText,newText);};
notes['2']=replaceChecked(notes['2'],
 'Gioelli et al. studied the mini-WARS splice form, which inhibits NRP1-dependent VE-cadherin turnover and can reduce endothelial permeability. This supports resolving WARS1 forms before assigning vascular function. The Olink measurement in our analysis does not identify full-length, spliced, cleaved or vesicular species.',
 'Vascular evidence: Gioelli et al. found that mini-WARS binds NRP1 and inhibits NRP1/VE-cadherin internalization, stabilizing endothelial junctions and reducing histamine-induced permeability in experimental endothelial cells. Tzima et al. found that the proteolytic T2-WARS fragment binds VE-cadherin, inhibits VEGF-induced ERK activation and endothelial migration, and inhibits angiogenesis. Mini-WARS and T2-WARS are distinct molecular forms. The NRP1 mechanism should not be assigned to T2-WARS. These studies did not test exercise plasma. The Olink result does not identify full-length, spliced, cleaved or vesicular WARS1 species.');
notes['2']+='\n\n'+t2ref+'\n\nInterpretation: vascular signaling is the prioritized follow-up hypothesis because extracellular WARS1 forms have direct experimental vascular actions. MoTrPAC does not establish that this explanation is more likely than other functions. Reduced vessel leakage and inhibition of new vessel growth are different outcomes, and neither implies a uniformly beneficial exercise adaptation. The intracellular insulin-receptor mechanism remains background biology and does not establish a consequence of the plasma response.';
p.resolve('sl/hwbqtkby').speakerNotes.textFrame.setText(notes['2']);

notes['5']=replaceChecked(notes['5'],
 'These hypotheses and experiments are proposed, not performed. They do not establish muscle origin, exercise-mode specificity or a protective effect on diabetes.',
 'These hypotheses and experiments are proposed, not performed. Vascular signaling is the prioritized WARS1 follow-up. Current results do not establish muscle origin, exercise-mode specificity or a beneficial vascular effect.');
notes['5']=replaceChecked(notes['5'],
 'Functional test: compare matched pre/post-exercise plasma under sham depletion, WARS1 depletion and add-back of a verified molecular form at an appropriate measured concentration. Endothelial permeability or innate immune readouts are motivated by external studies. A response that persists after effective specific depletion would weaken a WARS1-dependent mechanism. Include plasma-volume, cell-injury, assay-specificity and endotoxin controls.',
 'Vascular hypothesis and test: compare matched pre/post-exercise plasma under sham depletion, WARS1 depletion and add-back of the molecular form actually detected, at its measured concentration. Mini-WARS motivates testing reduced permeability and preservation of junctional VE-cadherin, with NRP1 dependence. T2-WARS motivates testing reduced VEGF-driven migration or sprouting. These predictions depend on detecting the relevant form. Test both barrier function and sprouting because they are distinct vascular outcomes. If the effect persists after effective specific depletion, or physiological add-back fails to restore it, that would weaken a WARS1-dependent explanation. Include matched resting controls, plasma-volume and cell-injury checks, and assay-specificity and endotoxin controls. Resolving forms and testing function require additional measurements.');
notes['5']=replaceChecked(notes['5'],
 'Available-data check: use the existing tissue/timepoint summary statistics for pathway-level tests of translation/stress and immune programs. This can support prioritization but cannot establish secretion, between-person coupling, receptor activation or disease mediation. Participant-level correlations require data that were not used here.',
 'Available-data follow-up, proposed: test whether endothelial junction and angiogenesis gene sets change in the available tissue/timepoint differential summaries, alongside translation/stress programs. Use each assay’s measured features as its background and correct across the tested gene sets. Such enrichment is contextual evidence from bulk tissues. It cannot identify the active plasma WARS1 form, prove endothelial cell specificity or establish causality. The current summary data cannot test between-person coupling or mediation.');
notes['5']+='\n\n'+t2ref;
p.resolve('sl/i107q5of').speakerNotes.textFrame.setText(notes['5']);

if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed during editing. Reconcile before exporting.');
const candidate=path.join(TMP,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,
 finalPath:path.join(OUT,'fractalkine_wars1_vascular_hypothesis.pptx'),
 explicitTotalSlideCount:5,
 pythonExecutable:path.join(RUNTIME,'python/bin/python3'),
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','5'],
 requiredNativeTableOwnerSlides:[5],requiredNativeChartOwnerSlides:[],
 fontPolicy:{basis:'reference',families:['Helvetica Neue'],referencePath:sourcePath,referenceSha256:sourceHash},
 verifyArtifactToolImport:true,receiptPath:path.join(TMP,'validation.json')
});
console.log(JSON.stringify({finalPath:result.finalPath,slideCount:p.slides.items.length}));
