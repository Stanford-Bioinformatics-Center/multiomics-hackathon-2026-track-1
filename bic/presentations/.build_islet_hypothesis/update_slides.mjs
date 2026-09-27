import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const ROOT='/Users/acheron/multiomics-hackathon-2026-track-1/bic/presentations';
const TMP=path.join(ROOT,'.build_islet_hypothesis');
const OUT=path.join(ROOT,'output');
const SKILL='/Users/acheron/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME='/Users/acheron/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=path.join(RUNTIME,'node/node_modules');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const sourcePath=path.join(TMP,'source.pptx');
const activeSource=path.join(OUT,'fractalkine_wars1_conceptual_hypotheses.pptx');
const sourceHash=(await fs.readFile(path.join(TMP,'source.sha256'),'utf8')).trim();
const hash=b=>createHash('sha256').update(b).digest('hex');
if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed. Reconcile edits before continuing.');
const p=await PresentationFile.importPptx(await FileBlob.load(sourcePath));
if(p.slides.items.length!==5) throw new Error('Unexpected slide count');
const set=(id,text)=>{p.resolve(id).text=text;};
set('sh/a10jqpsj','Human islets: less β-cell death and glucagon secretion, without increased insulin secretion in this study. [1]');
p.resolve('sh/a10jqpsj').position={left:60,top:255,width:540,height:110};
set('sh/x4r21kru','Rat β-cells: CX3CL1 protects insulin secretion against impairment by the inflammatory cytokine TNFα. [1]');
p.resolve('sh/x4r21kru').position={left:60,top:380,width:540,height:100};
set('sh/w3i1sfa9','An additional possible role after exercise: immune-cell communication and vascular remodeling. [7]');
p.resolve('sh/w3i1sfa9').position={left:60,top:510,width:540,height:110};
set('sh/mtwrmxg7','[1] Rutti 2014   [4] Gioelli 2022   [6] Tzima 2005   [7] Strömberg 2016. Full references in notes.');
const t=p.resolve('tb/adonexsj');
t.cells.set(1,1,'Type 2 diabetes hypothesis: the CX3CL1 pulse during endurance could signal to pancreatic islets, supporting β-cell survival and modulating glucagon during recovery. Repeated pulses might help preserve glucose control. [1]');
set('sh/yhkbe1o7','Direct EE–RE plasma tests were not significant. No during-resistance samples were collected.\nPancreatic and vascular effects remain untested. These data do not establish a diabetes benefit.');
const notes=JSON.parse(await fs.readFile(path.join(TMP,'source-notes.json'),'utf8'));
const immuneRef='[7] Strömberg A et al. CX3CL1—a macrophage chemoattractant induced by a single bout of exercise in human skeletal muscle. American Journal of Physiology-Regulatory, Integrative and Comparative Physiology 2016;310:R297–R304. https://doi.org/10.1152/ajpregu.00236.2015 ; https://pubmed.ncbi.nlm.nih.gov/26632602/';
notes['2']+='\n\n'+immuneRef+'\n\nComplementary biology: the human exercise study found increased muscle CX3CL1 with predominantly endothelial localization. Cell experiments supported immune communication and remodeling as possible functions. This is an alternative to, or could coexist with, pancreatic signaling. Neither study establishes which role dominates after exercise. The slide distinguishes the human-islet observations from the rat TNFα experiments.';
p.resolve('sl/hwbqtkby').speakerNotes.textFrame.setText(notes['2']);
const oldConcept='CX3CL1 hypothesis: sustained metabolic activity could prompt local muscle/endothelial communication, with CX3CL1 helping coordinate immune cells and vascular remodeling around active tissue. Strömberg et al. found increased muscle CX3CL1 after cycling, mainly endothelial localization, and induction in endothelial cells exposed to exercised-muscle tissue fluid. In cell experiments CX3CL1 altered inflammatory, chemotactic and angiogenic programs. These observations motivate the conceptual role. They do not identify the source or function of the early MoTrPAC plasma signal. Rapid shedding is one possible release mechanism. Local signaling with spillover into plasma is an alternative to a systemic endocrine role.';
if(!notes['5'].includes(oldConcept)) throw new Error('Expected conceptual note missing.');
const newConcept=[
 'Primary CX3CL1 conceptual hypothesis for type 2 diabetes: sustained metabolic demand during endurance exercise could be accompanied by a circulating signal to pancreatic islets. The proposed downstream role is support of β-cell survival and modulation of glucagon during recovery. If repeated exercise pulses have relevant activity, they might contribute to maintaining glucose regulation over time. This is an inference combining the acute MoTrPAC association with separate islet experiments. It is not a demonstrated explanation for the plasma rise or evidence that CX3CL1 prevents diabetes.',
 'Evidence boundary: the human-islet study used prolonged exposure in culture, so it does not establish the activity of an acute exercise pulse. It did not find increased insulin secretion under those conditions. Its TNFα-related insulin-secretion rescue primarily concerns rat β-cells. The hypothesis therefore emphasizes islet survival and glucagon regulation rather than assuming more insulin release. Dose, timing and glucose context need testing. Suppression of glucagon should not be assumed beneficial in every physiological context.',
 'Alternative CX3CL1 interpretation: immune-cell communication and vascular remodeling around exercised muscle remain plausible, based on Strömberg et al. These local roles could coexist with pancreatic signaling, or the plasma rise could reflect spillover without a systemic action. The tissue source, release mechanism and target organ remain unresolved.',
 'Discriminating diabetes experiment, proposed: quantify the exercise-associated CX3CL1 exposure using an orthogonal assay. Apply the measured pulse to human islets under relevant glucose conditions and compare β-cell survival, glucagon and glucose-stimulated insulin secretion with CX3CL1 depletion or CX3CR1 blockade and appropriate controls. Repeated-pulse and longer-term studies would be needed to test the extended hypothesis about glucose regulation. These experiments require new measurements and were not performed here.'
].join('\n\n');
notes['5']=notes['5'].replace(oldConcept,newConcept);
p.resolve('sl/i107q5of').speakerNotes.textFrame.setText(notes['5']);
if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed during editing. Reconcile before export.');
const candidate=path.join(TMP,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,finalPath:path.join(OUT,'fractalkine_wars1_islet_vascular_hypotheses.pptx'),
 explicitTotalSlideCount:5,pythonExecutable:path.join(RUNTIME,'python/bin/python3'),
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','5'],
 requiredNativeTableOwnerSlides:[5],requiredNativeChartOwnerSlides:[],
 fontPolicy:{basis:'reference',families:['Helvetica Neue'],referencePath:sourcePath,referenceSha256:sourceHash},
 verifyArtifactToolImport:true,receiptPath:path.join(TMP,'validation.json')
});
console.log(JSON.stringify({finalPath:result.finalPath,slideCount:p.slides.items.length}));
