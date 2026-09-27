import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {FileBlob,PresentationFile} from '@oai/artifact-tool';
const ROOT='/Users/acheron/multiomics-hackathon-2026-track-1/bic/presentations';
const TMP=path.join(ROOT,'.build_conceptual_hypotheses');
const OUT=path.join(ROOT,'output');
const SKILL='/Users/acheron/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME='/Users/acheron/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=path.join(RUNTIME,'node/node_modules');
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const sourcePath=path.join(TMP,'source.pptx');
const activeSource=path.join(OUT,'fractalkine_wars1_vascular_hypothesis.pptx');
const sourceHash=(await fs.readFile(path.join(TMP,'source.sha256'),'utf8')).trim();
const hash=b=>createHash('sha256').update(b).digest('hex');
if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed. Reconcile edits before continuing.');
const p=await PresentationFile.importPptx(await FileBlob.load(sourcePath));
if(p.slides.items.length!==5) throw new Error('Unexpected slide count');
const set=(id,text)=>{p.resolve(id).text=text;};
set('sh/032tgr6d','Conceptual hypotheses');
set('sh/dgbulwnm','Why might release occur during endurance or soon after resistance exercise?');
const t=p.resolve('tb/adonexsj');
t.columns.get(0).width=298;
t.columns.get(1).width=862;
t.cells.set(0,1,'Conceptual hypotheses');
t.cells.set(1,0,'CX3CL1\nPlasma rises during endurance exercise\n(20–40 min)');
t.cells.set(2,0,'WARS1\nPlasma rises after resistance exercise\n(10 min)');
t.cells.set(1,1,'Sustained metabolic demand may trigger CX3CL1 release from muscle or its endothelium. The signal could recruit immune cells and coordinate vascular remodeling to help active tissue recover and adapt. [7]');
t.cells.set(2,1,'Mechanical stress from heavy contractions could prompt release of stored WARS1. During recovery, mini-WARS could limit vascular leakage, while T2-WARS could restrain new vessel growth. The active form is unknown. [4,6,8]');
set('sh/yhkbe1o7','Direct EE–RE plasma tests were not significant. No during-resistance samples were collected.\nA shared response with different timing remains possible. Release and function need confirmation.');
const prior=JSON.parse(await fs.readFile(path.join(TMP,'source-notes.json'),'utf8'))['5'];
const conceptual=[
 'Conceptual interpretation: the table proposes possible triggers and adaptive roles. The exercise-mode labels describe where this analysis detected a plasma response. They do not establish mode-specific release, a beneficial effect, or an evolutionary purpose.',
 'CX3CL1 hypothesis: sustained metabolic activity could prompt local muscle/endothelial communication, with CX3CL1 helping coordinate immune cells and vascular remodeling around active tissue. Strömberg et al. found increased muscle CX3CL1 after cycling, mainly endothelial localization, and induction in endothelial cells exposed to exercised-muscle tissue fluid. In cell experiments CX3CL1 altered inflammatory, chemotactic and angiogenic programs. These observations motivate the conceptual role. They do not identify the source or function of the early MoTrPAC plasma signal. Rapid shedding is one possible release mechanism. Local signaling with spillover into plasma is an alternative to a systemic endocrine role.',
 'WARS1 hypothesis: high-force contractions and associated vascular compression could provide a mechanical-stress context for release of existing WARS1 from stressed cells during early recovery. The proposed link from that stress to WARS1 release is an inference, not an experimentally demonstrated mechanism. Mini-WARS motivates a barrier-stabilization hypothesis, whereas T2-WARS motivates restraint of vessel sprouting. The latter is not equivalent to better perfusion or a beneficial training response. Detecting and identifying the extracellular form is essential. Later muscle RNA induction cannot explain the earlier plasma rise by itself or identify its source.',
 'Mode and timing: neither candidate passed BH in a direct post-exercise EE–RE plasma comparison. A significant exercise-versus-control result in one arm and a nonsignificant result in the other is not evidence of a significant difference between arms. There are no during-resistance samples. A shared response with different timing, clearance changes, concentration changes or cell injury remains possible. WARS1 is a candidate vascular exerkine rather than an established exercise-mediated vascular signal. The claim of novelty requires a systematic literature assessment.',
 '[7] Strömberg A et al. CX3CL1—a macrophage chemoattractant induced by a single bout of exercise in human skeletal muscle. American Journal of Physiology-Regulatory, Integrative and Comparative Physiology 2016;310:R297–R304. https://doi.org/10.1152/ajpregu.00236.2015 ; https://pubmed.ncbi.nlm.nih.gov/26632602/',
 '[8] Muscle contraction-blood flow interactions during upright knee extension exercise in humans. Journal of Applied Physiology 2004. https://doi.org/10.1152/japplphysiol.00219.2004 . Physiological context only: contraction-related intramuscular pressure can restrict flow. This study does not establish WARS1 release or compare its response between exercise modes.',
 'Further rationale, proposed tests and sources follow. These tests have not been performed in this analysis.'
].join('\n\n');
p.resolve('sl/i107q5of').speakerNotes.textFrame.setText(conceptual+'\n\n'+prior);
if(hash(await fs.readFile(activeSource))!==sourceHash) throw new Error('Source changed during editing. Reconcile before export.');
const candidate=path.join(TMP,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
const result=await finalizePresentation({
 workspaceDir:ROOT,candidatePath:candidate,finalPath:path.join(OUT,'fractalkine_wars1_conceptual_hypotheses.pptx'),
 explicitTotalSlideCount:5,pythonExecutable:path.join(RUNTIME,'python/bin/python3'),
 integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','5'],
 requiredNativeTableOwnerSlides:[5],requiredNativeChartOwnerSlides:[],
 fontPolicy:{basis:'reference',families:['Helvetica Neue'],referencePath:sourcePath,referenceSha256:sourceHash},
 verifyArtifactToolImport:true,receiptPath:path.join(TMP,'validation.json')
});
console.log(JSON.stringify({finalPath:result.finalPath,slideCount:p.slides.items.length}));
