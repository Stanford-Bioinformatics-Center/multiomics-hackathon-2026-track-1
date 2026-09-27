import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { Presentation, PresentationFile } from '@oai/artifact-tool';

const ROOT = '/Users/acheron/multiomics-hackathon-2026-track-1/bic';
const workspaceDir = path.join(ROOT, 'presentations');
const TMP = path.join(workspaceDir, '.build_cx3cl1_wars1');
const OUT = path.join(workspaceDir, 'output');
const SKILL = '/Users/acheron/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME = '/Users/acheron/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES = path.join(RUNTIME, 'node/node_modules');
const { resolvePresentationFont, finalizePresentation } = await import(pathToFileURL(path.join(SKILL, 'container_tools/artifact_tool_utils.mjs')).href);
const font = resolvePresentationFont();
const presentation = Presentation.create({ slideSize: { width: 1280, height: 720 } });
const C = { ink: '#172F40', blue: '#2267A4', orange: '#AC5226', muted: '#526571', paper: '#FAFBFC', light: '#EEF3F6', white: '#FFFFFF' };
await fs.mkdir(OUT, { recursive: true });
const data = JSON.parse(await fs.readFile(path.join(TMP, 'data.json'), 'utf8'));
const sourceURL = 'https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77';
const commonSource = `Source: MoTrPAC human acute exercise, collection c2.0, analysis package 2.0.8, commit 535b4044e7417413de471104c619120337602b77. ${sourceURL}\nThe saved notebooks reuse MoTrPAC mixed-model estimates, raw p-values and confidence intervals. No new participant-level model was fitted. BH-adjusted p-values were independently reproduced. Primary comparisons have contrast_type=exercise_with_controls. Effects are source logFC, not absolute concentrations or asserted percentage changes.`;
const slideCopy = [];

function txt(slide, text, x, y, w, h, size=25, color=C.ink, bold=false) {
  const s = slide.shapes.add({ geometry: 'textbox', position: { left:x, top:y, width:w, height:h }, fill:'none', line:{ fill:'none', width:0 } });
  s.text = text;
  s.text.style = { typeface:font, fontSize:size, color, bold, autoFit:'none' };
  return s;
}
function base(title, n) {
  const s = presentation.slides.add();
  s.background.fill = C.paper;
  txt(s, title, 60, 42, 1160, 106, 44, C.ink, true);
  txt(s, `${n} / 4`, 1164, 665, 68, 30, 16, C.muted);
  return s;
}
function note(slide, detail) { slide.speakerNotes.textFrame.setText(`${commonSource}\n\n${detail}`); }
function table(slide, values, widths, y, height) {
  const t = slide.tables.add({ rows:values.length, columns:values[0].length, left:60, top:y, width:1160, height, columnWidths:widths, values });
  t.styleOptions = { headerRow:true, bandedRows:false };
  t.borders.assign({ fill:'#D8E1E6', width:.7, style:'solid' });
  const all = t.cells.block({ row:0, column:0, rowCount:values.length, columnCount:values[0].length });
  all.assign({ fill:C.white, textStyle:{ typeface:font, fontSize:23, color:C.ink }, margins:{ left:14, right:12, top:13, bottom:10 }, anchor:'center' });
  t.cells.block({ row:0, column:0, rowCount:1, columnCount:values[0].length }).assign({ fill:C.ink, textStyle:{ typeface:font, fontSize:23, color:C.white, bold:true } });
  for (let r=1; r<values.length; r++) if(r%2===0) t.cells.block({row:r,column:0,rowCount:1,columnCount:values[0].length}).fill = C.light;
  return t;
}
const fmt = x => (x>=0?'+':'')+x.toFixed(3);
const ci = r => `[${r['CI.L_calculated'].toFixed(3)}, ${r['CI.R_calculated'].toFixed(3)}]`;

// 1. The estimand and correction, in plain language.
{
  const s=base('Human MoTrPAC differential analysis',1);
  txt(s,'Acute endurance or resistance exercise in sedentary adults',60,154,1160,38,28,C.blue,true);
  txt(s,'Effect = change from baseline in exercise − change from baseline in resting controls',60,218,1160,84,32,C.ink,true);
  txt(s,'Null hypothesis: this difference in changes equals zero for the molecule and timepoint.',60,311,1160,58,26);
  txt(s,'MoTrPAC supplied mixed-model effects and raw p-values. We extracted each candidate and reproduced Benjamini–Hochberg (BH) correction.',60,394,1160,82,26);
  txt(s,'A BH-adjusted p < 0.05 defines a hit. BH accounts for all features within each assay and comparison, including 1,417 Olink features for each plasma comparison.',60,492,1160,88,26);
  txt(s,'BH targets a 5% expected false-discovery fraction among reported hits, under its assumptions. It does not change the null hypothesis.',60,605,1080,64,21,C.muted);
  note(s,`Analysis notebooks: ${ROOT}/fractalkine_differential_analysis.ipynb and ${ROOT}/wars1_differential_analysis.ipynb.\nPlasma model: ~ 0 + group_timepoint + BMI + calculatedAge + Sex + (1 | pid). RNA models add assay-specific covariates. BH family: tissue, assay, platform, contrast. Source correction applies separately to each comparison, not jointly across all timepoints.\nThe study tests responses to a single bout, not adaptation to a training program. Positive effects mean a larger baseline-to-timepoint change than control. The raw p-value assesses H0: the model contrast is zero. BH is a multiple-testing procedure, not another differential model. FDR is an expectation over repetitions, not a posterior probability that an individual hit is false. Pointwise confidence intervals do not include multiplicity correction.`);
  slideCopy.push({title:'Human MoTrPAC differential analysis', paragraphs:[
    'Data: human MoTrPAC c2.0 acute exercise in sedentary adults, analysis package 2.0.8.',
    'Effect = (exercise at timepoint − exercise baseline) − (control at timepoint − control baseline).',
    'Null hypothesis: this difference in changes equals zero for the molecule and timepoint.',
    'MoTrPAC supplied mixed-model effects and raw p-values. We extracted candidates and independently reproduced BH.',
    'A BH-adjusted p < 0.05 defines a hit. BH accounts for all features in each assay/comparison, including 1,417 plasma Olink features. It targets a 5% expected false-discovery fraction among reported hits under its assumptions.'
  ]});
}
// 2. CX3CL1.
{
  const s=base('Fractalkine: endurance and resistance',2);
  s.images.add({blob:new Uint8Array(await fs.readFile(path.join(OUT,'figures/cx3cl1_both_modes.png'))),contentType:'image/png',fit:'contain',alt:'CX3CL1 plasma and muscle RNA control-adjusted effects across time, with pointwise 95% confidence intervals and BH significance.',position:{left:60,top:148,width:1160,height:426.5}});
  txt(s,'Resistance: no plasma BH hits after exercise. Muscle RNA rises at 15 min: +2.585, q = 1.55 × 10⁻³³.\nDuring-resistance plasma was not sampled.',60,583,1160,60,23);
  txt(s,'Filled: q < 0.05. Hollow: q ≥ 0.05. Bars: pointwise 95% CIs. q = source BH-adjusted p.\nTimepoints are equally spaced. No during-resistance samples. RNA and protein scales differ.',60,648,1080,48,18.5,C.muted);
  const rows=data.cxPlasma.filter(r=>r.adj_p_value<.05);
  const vals=[['Endurance vs control','Effect','95% CI','Raw p','BH-adjusted p'],...rows.map((r,i)=>[`During ${i===0?'20':'40'} min`,fmt(r.logFC),ci(r),i===0?'1.29 × 10⁻⁶':'2.89 × 10⁻⁶',i===0?'0.000364':'0.000227'])];
  note(s,`Sources: ${ROOT}/results/fractalkine/cx3cl1_plasma_differential.csv and cx3cl1_tissue_rna_differential.csv. All numerical values are rounded from the CSVs.\nResistance plasma: during-exercise samples are unavailable. Post 10 min effect +0.156882, raw p=0.146200, BH-adjusted p=0.548753. The remaining resistance plasma adjusted p-values are 0.595239 (30 min), 0.962028 (3.5 h), and 0.791810 (24 h). None passes BH. Exactly two of ten primary plasma comparisons pass source BH. The remaining post-exercise plasma contrasts do not. Holm-adjusted p-values across the ten focal plasma comparisons are 0.0000128585 and 0.0000259997. Selection followed exploratory screening, so these sensitivity checks are not independent validation.\nMuscle RNA at post 15 min: EE-CON effect +3.626050, q=5.500778e-53. RE-CON effect +2.585393, q=1.554822e-33. There are four of twelve significant muscle/adipose primary RNA comparisons in total. The other two are endurance muscle RNA at 3.5 h (effect +0.536484, q=0.041645) and endurance adipose RNA at 45 min (effect +0.587781, q=0.027316). RNA and protein effect scales must not be pooled.\nDuring-resistance blood was not sampled. Zero of four direct post-exercise EE-RE plasma comparisons pass BH (minimum adjusted p approximately 0.472). An earlier circulating response plus later tissue RNA does not identify the tissue source.`);
  slideCopy.push({title:'Fractalkine: endurance and resistance',table:vals,paragraphs:[
    'CX3CL1 plasma protein, Olink feature OID20976. Effects use source logFC and confidence intervals are pointwise.',
    'Resistance plasma: no during-exercise samples. Post-exercise BH-adjusted p-values are 0.549 (10 min), 0.595 (30 min), 0.962 (3.5 h), and 0.792 (24 h). None passes BH.',
    '2 of 10 plasma exercise-versus-control tests pass BH < 0.05. Both also pass Holm correction across the ten CX3CL1 plasma tests.',
    'Muscle RNA rises 15 min after endurance exercise (+3.626, BH p = 5.50 × 10⁻⁵³) and resistance exercise (+2.585, BH p = 1.55 × 10⁻³³).',
    'No post-exercise plasma test passes BH < 0.05. RNA and protein effects are on different assay scales.'
  ]});
}
// 3. WARS1.
{
  const s=base('WARS1: endurance and resistance',3);
  s.images.add({blob:new Uint8Array(await fs.readFile(path.join(OUT,'figures/wars1_both_modes.png'))),contentType:'image/png',fit:'contain',alt:'WARS1 plasma and muscle RNA control-adjusted effects across time. Plasma protein is significant at 10 minutes after resistance exercise; muscle RNA is significant at 3.5 hours after both exercise modes.',position:{left:60,top:148,width:1160,height:426.5}});
  txt(s,'Endurance plasma: q = 0.482 at 20 min and 0.136 at 40 min. Neither passes BH.\nEndurance muscle RNA rises at 3.5 h: +0.374, q = 1.48 × 10⁻⁶.',60,583,1160,60,23);
  txt(s,'Filled: q < 0.05. Hollow: q ≥ 0.05. Bars: pointwise 95% CIs. q = source BH-adjusted p.\nTimepoints are equally spaced. No during-resistance samples. RNA and protein scales differ.',60,648,1080,48,18.5,C.muted);
  const hits=data.warsPrimary.filter(r=>r.adj_p_value<.05);
  const vals=[['Measurement','Exercise vs control','Time after exercise','Effect [95% CI]','BH-adjusted p'],...hits.map(r=>[r.tissue==='blood'?'Plasma protein':'Muscle RNA',r.contrast_category==='EE-CON'?'Endurance':'Resistance',r.tissue==='blood'?'10 min':'3.5 h',`${fmt(r.logFC)}\n${ci(r)}`,r.tissue==='blood'?'0.0312':r.contrast_category==='EE-CON'?'1.48 × 10⁻⁶':'3.63 × 10⁻⁶'])];
  note(s,`Sources: ${ROOT}/results/wars1/wars1_primary.csv, wars1_plasma_primary.csv, wars1_direct_EE_RE.csv, wars1_post10_plasma_decomposition.csv and FINDINGS.md.\nEndurance plasma at 20 min: effect +0.247576, raw p=0.154284, BH-adjusted p=0.482228. At 40 min: effect +0.387097, raw p=0.030154, BH-adjusted p=0.135725. Neither passes BH. Endurance muscle RNA at 3.5 h: effect +0.374024, q=1.479715e-6. Resistance muscle RNA at 3.5 h: effect +0.329233, q=3.628728e-6. Plasma Olink feature OID21084 maps to WARS1 (aliases WARS/WRS/IFI53), UniProt P23381. The plasma raw p-value is 0.0007337048817773. Its pointwise CI is [0.2417157128, 0.8919410064]. Within resistance, post 10 min vs baseline is +0.370578 (q=0.008118). Control post 10 min vs baseline is -0.196250 (q=0.751245). Their difference gives the primary effect +0.566828.\nThe plasma response also passes the earlier exploratory pooled BH across 14,170 primary plasma tests (adjusted p=0.040771). This is a sensitivity analysis on the same data. No direct WARS1 contrast passes source BH. No during-resistance samples exist.\nFifty-two primary estimates span available assays, including two muscle phosphosites. Neither the other compartments nor muscle total protein pass primary source BH. A nonsignificant total protein result does not prove absent release. Muscle RNA at 3.5 h cannot establish the origin of a 10 min plasma response. Absolute effects across RNA and protein scales are not directly comparable.`);
  slideCopy.push({title:'WARS1: endurance and resistance',table:vals,paragraphs:[
    'Early plasma protein response and later muscle RNA responses. Plasma Olink feature OID21084.',
    'Endurance plasma: during 20 min effect +0.248, q=0.482. During 40 min effect +0.387, raw p=0.0302, q=0.136. Neither passes BH. Endurance muscle RNA rises at 3.5 h (+0.374, q=1.48e-6).',
    'These are the only 3 BH hits among 52 primary WARS1 estimates. Muscle total protein does not show a significant exercise response.',
    'Direct endurance–resistance plasma comparison at 10 min: raw p = 0.0442, BH p = 0.227. Resistance specificity is not established.',
    'Effects use source logFC. CIs are pointwise. Muscle RNA at 3.5 h does not establish the origin of the earlier plasma response.'
  ]});
}
// 4. Inference and falsifiable follow-up.
{
  const s=base('Interpretation and testable next steps',4);
  txt(s,'Fractalkine / CX3CL1',60,166,540,42,30,C.blue,true);
  txt(s,'The plasma signal is clearest during endurance exercise.',60,224,540,76,28);
  txt(s,'Test: reproduce the 20–40 min plasma increase using an independent protein assay and matched resting controls.',60,328,540,116,26);
  txt(s,'WARS1',674,166,546,42,30,C.orange,true);
  txt(s,'The plasma signal occurs 10 min after resistance exercise.',674,224,546,76,28);
  txt(s,'Test: confirm the early plasma increase and determine which WARS1 molecular form the assay detects.',674,328,546,116,26);
  txt(s,'Shared limits',60,485,1160,42,28,C.ink,true);
  txt(s,'Neither plasma analysis establishes exercise-mode specificity. Later tissue RNA cannot identify the source of the earlier plasma signal.',60,541,1160,76,26);
  txt(s,'These exploratory candidates require replication. The current analysis does not test diabetes benefit.',60,640,1080,47,23,C.muted);
  note(s,`Sources: ${ROOT}/results/fractalkine/FINDINGS.md and ${ROOT}/results/wars1/FINDINGS.md. Proposed experiments are hypotheses, not completed results.\nA targeted orthogonal plasma assay should measure both exercise arms and time-matched controls, with baseline and during/post-exercise sampling, and prospectively specified contrasts and correction. Distinguishing secretion, cleavage, altered clearance and plasma-volume effects requires further measurements. The available relative Olink assay results do not establish absolute concentrations or the molecular forms recognized.\nDuring-resistance samples are unavailable in the source study, and direct post-exercise plasma EE-RE comparisons do not pass BH for either gene. A result significant in one arm and nonsignificant in another is not itself evidence that the arms differ.\nThe current analysis supplies candidate prioritization from exercise-responsive measurements. It does not measure tissue-to-plasma secretion flux, mediation of clinical outcomes, or prevention/treatment of diabetes.`);
  slideCopy.push({title:'Interpretation and testable next steps',paragraphs:[
    'Fractalkine: reproduce the 20–40 min endurance plasma increase using an independent protein assay and matched resting controls.',
    'WARS1: confirm the early resistance plasma increase and determine which WARS1 molecular form the assay detects.',
    'Neither plasma analysis establishes exercise-mode specificity. During-resistance samples are unavailable and the direct post-exercise plasma comparisons do not pass BH.',
    'Later tissue RNA cannot identify the source of an earlier plasma signal. These exploratory candidates require replication. The current analysis does not test diabetes benefit.'
  ]});
}

let md='# Fractalkine and WARS1 differential analysis\n\nFour-slide presentation with scientific time-course figures, based on the available human MoTrPAC c2.0 acute-exercise results. Figures are available as PNG and SVG in the adjacent figures directory. Values and detailed statistical results follow below.\n\n';
for(let i=0;i<slideCopy.length;i++){
  const s=slideCopy[i]; md+=`## Slide ${i+1}: ${s.title}\n\n`;
  if(s.table){md+='| '+s.table[0].join(' | ')+' |\n| '+s.table[0].map(()=>'---').join(' | ')+' |\n'; for(const r of s.table.slice(1)) md+='| '+r.map(v=>String(v).replaceAll('\n',' ')).join(' | ')+' |\n'; md+='\n';}
  md+=s.paragraphs.map(p=>'- '+p).join('\n')+'\n\n';
}
md+='## Sources\n\n';
md+=`- [Pinned MoTrPAC analysis package](${sourceURL})\n- [Fractalkine notebook](../../fractalkine_differential_analysis.ipynb)\n- [WARS1 notebook](../../wars1_differential_analysis.ipynb)\n- [Fractalkine plasma estimates](../../results/fractalkine/cx3cl1_plasma_differential.csv)\n- [Fractalkine RNA estimates](../../results/fractalkine/cx3cl1_tissue_rna_differential.csv)\n- [WARS1 primary estimates](../../results/wars1/wars1_primary.csv)\n\nEffects, raw p-values and pointwise CIs originate from MoTrPAC. The notebooks independently verified source BH corrections. They did not refit participant-level models. The deck's speaker notes contain additional statistical detail and result provenance.\n`;
await fs.writeFile(path.join(OUT,'fractalkine_wars1_both_modes_slide_text.md'),md);
const candidatePath=path.join(TMP,'candidate_with_both_modes.pptx');
await (await PresentationFile.exportPptx(presentation)).save(candidatePath);
const result=await finalizePresentation({
  workspaceDir,candidatePath,finalPath:path.join(OUT,'fractalkine_wars1_with_both_modes.pptx'),
  explicitTotalSlideCount:4,
  pythonExecutable:path.join(RUNTIME,'python/bin/python3'),
  integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit'],
  requiredNativeTableOwnerSlides:[],requiredNativeChartOwnerSlides:[],fontPolicy:{basis:'design',families:[font]},
  verifyArtifactToolImport:true,receiptPath:path.join(TMP,'validation_both_modes.json')
});
console.log(JSON.stringify({font,result}));
