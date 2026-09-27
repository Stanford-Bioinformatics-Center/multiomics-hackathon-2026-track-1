# Documentation index: who documented what, and where

Team 2-PAC's documentation is spread over 5 branches and many folders. This page lists every README and
write-up in the repository, grouped by author, with a short summary and a link to the original.
It holds no copies: the linked files are the source of truth.

**As of 2026-09-27, `main` at `378e865`.** 34 documents (10 READMEs, 24 other write-ups) by 4 authors, not
counting this page. Authors are taken from `git blame` and shown as git records them. Blame shows who
committed the lines.

## Branches

| Branch | Commit | vs `main` | Documents it adds |
|---|---|---|---|
| `main` | `378e865` | — | 19 (8 READMEs, 11 write-ups) |
| `differential_analysis` | `95c6ad0` | 97 behind, 2 ahead | 15 under `bic/` (2 READMEs, 13 write-ups), on this branch only |
| `diabetes_data` | `e9fb2db` | 74 behind, 1 ahead | none |
| `random_walk` | `86e1e2c` | 48 behind, 0 ahead | none (fully merged into `main`) |
| `story-networks` | `2d8fa73` | 30 behind, 0 ahead | none (merged into `main`, PR #1) |

**Read the `main` version.** Every document that exists on both `main` and a side branch is either identical
or an older snapshot on the side branch (see [Older copies](#older-copies-on-side-branches)).

## Authors at a glance

| Author | READMEs | Other write-ups | Lines | Where |
|---|---|---|---|---|
| [Vidal M. Arroyo](#vidal-m-arroyo) | 6, one of them shared (the root README) | 3 | 1,815 | `main` |
| [jspaul2003](#jean-sebastien-paul) | 2 | 13 | 1,566 | `differential_analysis` |
| [gandhimonil9823](#gandhimonil9823) | 1 | 8 | 700 | `main` |
| [Subarna Bhattacharya](#subarna-bhattacharya) | 1 | 0 | 555 | `main` |

Lines are counted per author by `git blame`, so a shared document is split. The root README has 107 lines:
43 by Subarna Bhattacharya, 42 by Vidal M. Arroyo, and 22 by Jimmy Zhen (the organisers' track brief).

## Vidal M. Arroyo

The network pipeline and everything built on it. All on `main`.

| Document | Lines | What it covers |
|---|---|---|
| [README.md](../README.md) | 107 | Landing page: the headline result, the three top stories, how an edge is made, quick start. The quick start (requirements, inputs, run, what to do if it fails) is mostly by Subarna Bhattacharya. The last 22 lines are the organisers' track brief. |
| [network/README.md](../network/README.md) | 1,285 | **The main document.** Eight parts: snapshot, research question, workflow, setup, inputs and outputs, methods for every pipeline step (1 to 21), validation, reuse. Also holds the contributors table, roadmap and "honest roadblocks". |
| [network/engine/README.md](../network/engine/README.md) | 50 | `exnet`, the R package that implements the edge rule: 3 S4 classes, the functions, an example, 25 tests. |
| [network/neo4j/README.md](../network/neo4j/README.md) | 170 | Exporting the networks to a Neo4j graph: graph model, import, example queries, expected counts, how to extend it. The visualiser itself is not started. |
| [network/inventory/README.md](../network/inventory/README.md) | 81 | Counts of MoTrPAC phosphosite and GlyGen data available for the 471 proteins and 450 metabolites. Counts only; nothing is integrated. |
| [network/resource/README.md](../network/resource/README.md) | 52 | Column definitions of the two shareable feature lists (`proteins_471.csv`, `metabolites_450.csv`). The lists are not stored in the repo. |
| [network/ENVIRONMENT.md](../network/ENVIRONMENT.md) | 92 | Software versions and MD5 fingerprints of every external input. Written by `00_environment.R` on each run. |
| [network/docs/FIGURES.md](../network/docs/FIGURES.md) | 35 | Captions for Fig. 1, the method figure and figures 17, 19 and 20. |
| [network/docs/PRESENTATION.md](../network/docs/PRESENTATION.md) | 21 | 11-slide outline mapped to the judging criteria, plus the live demo steps. |

## gandhimonil9823

The interaction database ("mnet") that supplies the network's edges. All on `main`.

| Document | Lines | What it covers |
|---|---|---|
| [network/interaction_db/README.md](../network/interaction_db/README.md) | 184 | How to build the database: STRING v12 protein–protein links, Rhea protein–metabolite links, lipid classes, phosphosite and glycosite annotations. Covers outputs, the score rule, configuration and licences. |
| [network/interaction_db/reports/](../network/interaction_db/reports/) (8 files) | 516 | One build report per phase: [scaffold](../network/interaction_db/reports/mnet_phase0_scaffold.md), [downloads](../network/interaction_db/reports/mnet_phase1_downloads.md), [metabolite mapping](../network/interaction_db/reports/mnet_phase2_metabolites.md), [lipids](../network/interaction_db/reports/mnet_phase3_lipids.md), [pre-fixes](../network/interaction_db/reports/mnet_phase4_prefixes.md), [Rhea layer](../network/interaction_db/reports/mnet_phase4_rhea.md), [PTM annotation](../network/interaction_db/reports/mnet_phase6_ptm.md), [final coverage](../network/interaction_db/reports/mnet_final_coverage.md). Start with the last two. |

## Subarna Bhattacharya

| Document | Lines | What it covers |
|---|---|---|
| [random_walk/README.md](../random_walk/README.md) | 499 | `random_walk(start, arm)`: one weighted walk from any node to 3 others, using endurance or resistance edge weights. Covers the method, a worked example (PPIB), how to use the walks, validation, choices considered, known limits and next steps. |

## Jean-Sebastien Paul

Candidate exerkine analyses from published MoTrPAC results. **On branch `differential_analysis` only**, so the
links go to GitHub. Base path: `bic/`.

| Document | Lines | What it covers |
|---|---|---|
| [summary.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/summary.md) | 199 | **Start here.** The whole line of work, 25 to 27 September: Table 3 screen, 16-candidate audit, CCN1, the wider exerkine screen, WARS1, and the novelty correction. |
| [README.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/README.md) | 24 | Contents of the `bic/` folder and the dataset used. |
| [motrpac-exploration/README.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/motrpac-exploration/README.md) | 100 | First hypothesis shortlist. Recommends type 2 diabetes, with fractalkine (CX3CL1) as lead candidate and BCAAs as contrast. Ends with 8 analysis rules. |
| [motrpac-exploration/table3/REPORT.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/motrpac-exploration/table3/REPORT.md) | 107 | All 28 exerkines of Chow et al. 2022, Table 3, looked up in MoTrPAC: 16 matched in plasma, 2 respond (fractalkine, lactate). |
| [motrpac-exploration/sixteen_candidates/REPORT.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/motrpac-exploration/sixteen_candidates/REPORT.md) | 117 | Significance audit of those 16 over 230 plasma tests. FGF21 and IL-15 are secondary, exploratory findings. |
| [motrpac-exploration/portal_c2.0/ACCESS_AND_ANALYSIS_PLAN.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/motrpac-exploration/portal_c2.0/ACCESS_AND_ANALYSIS_PLAN.md) | 64 | Which MoTrPAC data are public and which are restricted, and the project recommended with the public data. |
| [results/fractalkine/FINDINGS.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/results/fractalkine/FINDINGS.md) | 30 | CX3CL1: plasma rises during endurance exercise (2 of 10 tests). Endurance-specificity is not established. |
| [results/ccn1/FINDINGS.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/results/ccn1/FINDINGS.md) | 59 | CCN1: 9 significant exercise-vs-control results across muscle, adipose and plasma. |
| [results/exerkine_screen/FINDINGS.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/results/exerkine_screen/FINDINGS.md) | 79 | Screen of all 1,417 plasma Olink features: 142 genes rise in plasma, 15 also rise in tissue. Priorities: CD300LG, ANGPT2, WARS1. |
| [results/wars1/FINDINGS.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/results/wars1/FINDINGS.md) | 85 | WARS1 follow-up: plasma rises 10 min after resistance exercise, muscle RNA at 3.5 h. Which molecular form rises is unresolved. |
| [results/wars1/novelty_review/REVIEW.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/results/wars1/novelty_review/REVIEW.md) | 77 | Correction: WARS1 is already in MoTrPAC's current candidate table, so it should not be called novel. |
| [presentations/output/](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/tree/differential_analysis/bic/presentations/output) (3 files) | 164 | Text of a four-slide fractalkine and WARS1 deck in three near-identical versions: plain, with figures, both exercise modes. |
| [motrpac-exploration/portal_c2.0/package_NEWS.md](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1/blob/differential_analysis/bic/motrpac-exploration/portal_c2.0/package_NEWS.md) | 461 | A saved copy of the MoTrPAC R package's changelog, kept for version provenance. Not written by the team. |

## Older copies on side branches

| Document | Older snapshot on | Lines there vs `main` |
|---|---|---|
| `README.md` | `diabetes_data`, `differential_analysis`, `random_walk` / `story-networks` | 23 / 66 vs 107 |
| `network/README.md` | `differential_analysis` / `diabetes_data` / `random_walk` / `story-networks` | 375 / 790 / 1,156 / 1,274 vs 1,285 |
| `network/ENVIRONMENT.md` | `random_walk` / `story-networks` | 83 / 90 vs 92 |
| `random_walk/README.md` | `random_walk` | 531 vs 499 |

Every other shared document is identical to `main`. Each older copy is an earlier version from `main`'s own
history; none has been edited separately on its branch.

**Removed from `main`:** the music video folder `network/video/`, with its README, was removed on 2026-09-27
(commit `378e865`); the tooling is kept outside the repository. A 53-line early version of that README is
still on `story-networks`. The last version on `main` had 125 lines and is in the history
(`git show 378e865~1:network/video/README.md`).

## Gaps and inconsistencies

Found while reading the documents. Each one is checked against the files.

| # | What | Where |
|---|---|---|
| 1 | **No README** for the three R scripts in `data_scripts/` | branch `diabetes_data` |
| 2 | The `bic/` work is **not mentioned in any document on `main`**, and is not merged | branch `differential_analysis` |
| 3 | The contributors table omits Subarna Bhattacharya and jspaul2003 | `network/README.md`, section 8 |
| 4 | **Two MoTrPAC package versions** are cited: 0.2.4 by the network pipeline, 2.0.8 by the `bic/` analyses | `network/README.md`; `bic/README.md` |
| 5 | **Outdated network size:** 364 nodes / 764 edges, where `network/README.md` now reports 353 / 704 | `network/neo4j/README.md`; `random_walk/README.md` |

Two earlier gaps are closed. The scripts `lyrics_gen.R` and `t2d_consensus.R`, which were documented but never
in the repository, are no longer mentioned in `random_walk/README.md`. The conflicting default for `--arm` was
in `network/video/README.md`, which is removed.

## Reuse and next steps

### What can be reused today

| To do this | Use | Documented in |
|---|---|---|
| Install every Python package the project uses | `pip install -r requirements.txt` (Python 3.10 to 3.12) | [requirements.txt](../requirements.txt), which also lists the R packages and system tools |
| Rebuild every table, figure and interactive page | `bash network/run_all.sh` (about 10 minutes) | [network/README.md](../network/README.md), section 5 |
| Build a weighted network from your own responses and edges | the `exnet` R package | [network/engine/README.md](../network/engine/README.md) |
| Rebuild the interaction database, at any STRING threshold | `make` targets in `network/interaction_db/` | [network/interaction_db/README.md](../network/interaction_db/README.md) |
| Walk the network from any protein or metabolite | `random_walk(start, arm)` | [random_walk/README.md](../random_walk/README.md) |
| Load the networks into a graph database | `bash network/neo4j/run_local_neo4j.sh` | [network/neo4j/README.md](../network/neo4j/README.md) |
| Look up the 471 proteins and 450 metabolites elsewhere | `export_feature_lists.R` | [network/resource/README.md](../network/resource/README.md) |

### Next steps

Collected from the documents' own roadmaps and stated limits, plus the gaps above. "Source" is the document
each step comes from.

| Area | Next step | Source |
|---|---|---|
| Documentation | Close the 5 gaps above: add a README for `data_scripts/`, update the contributors table, correct the network size | this page |
| Documentation | Merge `differential_analysis` into `main`, or link the `bic/` work from a document on `main` | this page |
| Documentation | Agree on one MoTrPAC package version, or state in both places why they differ | this page |
| Network | Test the blood result against blood cell composition and a plasma-volume control | `network/README.md`, roadmap |
| Network | Replicate the muscle result in an independent T2D muscle proteome | `network/README.md`, roadmap |
| Network | Run the same tests on MoTrPAC training data when released | `network/README.md`, roadmap |
| Network | Reinstate a test of arm differences for single edges against measurement noise | `network/README.md`, roadmap |
| Network | Extend beyond the 471 proteins that the Olink panel allows | `network/README.md`, roadmap |
| Network | Decide the normalisation (step 12 is still an open item) | `network/README.md`, step 12 |
| Random walk | Run many walks per start node and arm, and compare them with a null model | `random_walk/README.md`, section 15 |
| Random walk | Start walks from disease-associated genes | `random_walk/README.md`, section 15 |
| Random walk | Move the walker into the pipeline as a numbered step with validation checks | `random_walk/README.md`, section 15 |
| Neo4j | Build the visualiser (not started), then choose which inventory data to integrate | `network/neo4j/README.md`; `network/inventory/README.md` |
| Interaction database | Review the suggested targets for the 20 metabolites with no edge (not auto-applied) | `mnet_final_coverage.md` |
| Interaction database | Replace the placeholder scores for catalysis, transport and lipid-class edges | `network/interaction_db/README.md` |
| Exerkines | Add a versioned disease gene or variant resource and an independent disease-expression dataset | `bic/summary.md`, outstanding work |
| Exerkines | Confirm the plasma signals of CX3CL1, CD300LG, ANGPT2 and WARS1 with an independent protein assay | `bic/results/*/FINDINGS.md` |
| Exerkines | Resolve which molecular form of WARS1 rises after exercise | `bic/results/wars1/FINDINGS.md` |
| Exerkines | Resolve the chemical identity of BAIBA | `bic/summary.md`, outstanding work |

## References

Every publication cited in the 34 documents, copied as the source document cites it.
**Not checked against the journals:** verify before manuscript use. Where a document gives only a link and a
description, that description is kept.

Cited in: **N** `network/README.md` · **R** `random_walk/README.md` · **I** `network/interaction_db/README.md` ·
**B** `bic/` documents on `differential_analysis`.

### Network methods

- El-Kebir M et al. (2015) xHeinz. *Bioinformatics* 31:3147–3155. [N]
- Hamilton WL, Ying R, Leskovec J (2017) Representation learning on graphs. *IEEE Data Eng. Bull.* [N, R]
- Tang J et al. (2015) LINE: large-scale information network embedding. *WWW*. [N, R]
- Kipf TN, Welling M (2016) Variational graph auto-encoders. arXiv:1611.07308. [N, R]
- Hinton G, Vinyals O, Dean J (2015) Distilling the knowledge in a neural network. arXiv:1503.02531. [N, R]
- Guo C et al. (2017) On calibration of modern neural networks. *ICML*. [N, R]
- Gretton A et al. (2012) A kernel two-sample test. *JMLR* 13:723–773. [N, R]
- Perozzi, Al-Rfou & Skiena (2014) DeepWalk. *KDD*. [R]
- Grover & Leskovec (2016) node2vec. *KDD*. [R]
- Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E (2008) Fast unfolding of communities in large networks. *J Stat Mech* 2008:P10008. [N]
- Tukey JW (1977) *Exploratory Data Analysis*. [N]
- Heller & Yekutieli (2014), repfdr. Name and year only. [N]

### Statistics and enrichment

- Wu D, Smyth GK (2012) Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133. [N]
- Ritchie ME et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43:e47. [N]
- Rivals I, Personnaz L, Taing L, Potier MC (2007) Enrichment or depletion of a GO category within a class of genes: which test? *Bioinformatics* 23:401–407. [N]
- Subramanian A et al. (2005) Gene set enrichment analysis. *PNAS* 102:15545–15550. [N]
- Timmons JA, Szkop KJ, Gallagher IJ (2015) Multiple sources of bias confound functional enrichment analysis of global -omics data. *Genome Biol* 16:186. [N]
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M (2022) Urgent need for consistent standards in functional enrichment analysis. *PLoS Comput Biol* 18:e1009935. [N]

### Databases and resources

- Szklarczyk D et al. (2023) The STRING database in 2023. *Nucleic Acids Res* 51:D638–D646. [N, I]
- Bansal P et al. (2022) Rhea, the reaction knowledgebase in 2022. *Nucleic Acids Res* 50:D693–D700. [N, I]
- Türei D et al. (2021) Integrated intra- and intercellular signaling knowledge for multicellular omics analysis (OmniPath). *Mol Syst Biol* 17:e9923. [N, I]
- The UniProt Consortium (2023) UniProt. *Nucleic Acids Res*. [I]
- Hastings et al. (2016) ChEBI. *Nucleic Acids Res*. [I]
- Aimo et al. (2015) SwissLipids. *Bioinformatics*. [I]
- Fahy E, Subramaniam S (2020) RefMet. *Nat. Methods* 17:1173. [N]
- York WS et al. (2020) GlyGen: computational and informatics resources for glycoscience. *Glycobiology* 30:72–73. [N]
- Liberzon A et al. (2011) Molecular signatures database (MSigDB) 3.0. *Bioinformatics* 27:1739–1740. [N]
- Milacic M et al. (2024) The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Res* 52:D672–D678. [N]
- Agrawal A et al. (2024) WikiPathways 2024. *Nucleic Acids Res* 52:D679–D689. [N]
- Gene Ontology Consortium (2023) The Gene Ontology knowledgebase in 2023. *Genetics* 224:iyad031. [N]
- Rath S et al. (2021) MitoCarta3.0. *Nucleic Acids Res* 49:D1541–D1547. [N]
- [Exerkine Atlas](https://exerkineatlas.org/), used for candidate nomination. [B]

### MoTrPAC

- MoTrPAC (2026) human acute-exercise papers: skeletal muscle, subcutaneous adipose tissue and blood (bioRxiv / PMC). [N]
- MoTrPAC landscape preprint, cited with two links: [doi:10.64898/2026.02.27.702183](https://doi.org/10.64898/2026.02.27.702183) and [PMC13184684](https://pmc.ncbi.nlm.nih.gov/articles/PMC13184684/). [B]

### Disease and ageing datasets

- Amar et al. (2024) *Cell Metab* 36:1411. doi:10.1016/j.cmet.2023.12.021. Source of the 9 disease sets. [N]
- Kjærgaard J et al. (2025) Personalized molecular signatures of insulin resistance and type 2 diabetes. *Cell* 188:4106–4122. [N]
- Needham EJ et al. (2024) Personalized phosphoproteomics of skeletal muscle insulin resistance and exercise links MINDY1 to insulin action. *Cell Metab* 36:2542–2559. [N]
- Larsen JK et al. (2023) High-throughput proteomics uncovers exercise training and type 2 diabetes–induced changes in human white adipose tissue. *Sci Adv* 9:eadi7548. [N]
- Ubaida-Mohien C et al. (2019) Discovery proteomics in aging human skeletal muscle. *eLife* 8:e49874. [N]
- Sun et al. (2023) *Nature*. UK Biobank plasma proteome, association with age. [N]
- Gadd et al. (2024) *Nat Aging*. UK Biobank plasma proteome, incident T2D. [N]
- Name and year only, all from the Amar et al. 2024 sets: Öhman 2021 and Chae 2018 (T2D muscle) · Coats 2018 (HCM heart) · Havlenova 2021 (heart failure, rat) · Park 2019 (MI, mouse) · Niu 2022 (NASH and cirrhosis liver) · Yuan 2020 (NAFLD liver) · Stocks 2022 (ob/ob mouse liver). [N]

### Lipids and type 2 diabetes

- Liu C et al. (2009) Lactate inhibits lipolysis in fat cells through activation of GPR81. *J Biol Chem*. [N]
- Ahmed K et al. (2010) An autocrine lactate loop mediates insulin-dependent inhibition of lipolysis through GPR81. *Cell Metab*. [N]
- Bergman BC et al. (2015) Serum sphingolipids and changes with exercise. *Am J Physiol Endocrinol Metab*. [N]
- Wigger L et al. (2017) Plasma dihydroceramides are diabetes susceptibility biomarker candidates. *Cell Rep*. [N]
- Adams SH et al. (2009) Plasma acylcarnitine profiles in type 2 diabetes. *J Nutr*. [N]
- Plasma-ceramide cohort studies: *J Lipid Res* 2021, and EPIC-Potsdam, *Nat Commun* 2022. No authors given. [N]

### Exerkines and candidate biology

- Chow et al. (2022) [Exerkines in health, resilience and disease](https://www.nature.com/articles/s41574-022-00641-2). Source of the 28 Table 3 candidates. [B]
- CX3CL1: [human-islet study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4209359/) · [Cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3717389/). [B]
- BCAAs: [bidirectional genetic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10827349/). [B]
- Kynurenine: [Agudelo et al.](https://pubmed.ncbi.nlm.nih.gov/25259918/) [B]
- Neprilysin: [PARADIGM-HF trial](https://www.nejm.org/doi/full/10.1056/NEJMoa1409077). [B]
- CCN1: [human exercise study](https://pubmed.ncbi.nlm.nih.gov/29924476/) · [human retinopathy study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10273100/) · [retinal mechanism study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11097549/). [B]
- CD300LG: [Lee-Ødegård et al., 2024](https://pubmed.ncbi.nlm.nih.gov/39190027/). [B]
- ANGPT2: [An et al., 2017](https://elifesciences.org/articles/24071). [B]
- FLT3LG: [FLT3LG study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9999360/). ANGPTL7: [ANGPTL7 study](https://doi.org/10.1371/journal.pone.0173024). [B]

### WARS1

- [Nguyen et al., 2023](https://pubmed.ncbi.nlm.nih.gov/36640342/). Direct and vesicular secretion. [B]
- [TLR/TREM-1 signaling](https://doi.org/10.3390/biom10091283). [B]
- [Gioelli et al., 2022](https://doi.org/10.1038/s41467-022-31904-1). Mini-WARS and NRP1. [B]
- [Tzima et al., 2005](https://doi.org/10.1074/jbc.C400431200). T2-WARS and VE-cadherin. [B]
- [Sun et al., 2024](https://doi.org/10.1007/s00018-023-05082-2). Insulin-receptor tryptophanylation. [B]
- [Uluvar et al., 2026](https://doi.org/10.1007/s00125-026-06800-8). Disease and tissue-QTL tables 16–17. [B]
- [Tokolyi et al., 2025](https://doi.org/10.1038/s41588-025-02096-3). Splicing QTL and hypertension. [B]
- [Nurbekov et al., 2015](https://genescells.ru/2313-1829/article/view/120500). Athlete overtraining pilot. [B]
- [Aas et al., 2023](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2023.1143966/full). Human muscle-cell vesicles. [B]
- [Abbasi et al., 2025](https://doi.org/10.1016/j.chest.2025.07.3648). CHEST conference abstract. [B]

## How this was generated

Compiled by Tia Kohir with Anthropic's Claude (Claude Code). Every listed document was read in full, except
`package_NEWS.md` (first 40 lines).

Refreshed the same day after `network/video/` was removed from `main` (`378e865`): document counts, line
counts, authorship, branch positions, older copies and gaps were recomputed from git. The documents that
changed on `main` in between were checked against their diffs, not read again in full.

To refresh the inventory, from the repo root:

```bash
git fetch origin
for ref in origin/main origin/diabetes_data origin/differential_analysis origin/random_walk origin/story-networks; do
  git ls-tree -r --name-only "$ref" | grep -iE '\.md$' | while read -r path; do
    lines=$(git show "$ref:$path" | wc -l)
    authors=$(git blame --line-porcelain "$ref" -- "$path" | grep '^author ' | sort | uniq -c | sort -rn | tr -s ' \n' ' ')
    echo "$ref | $path | $lines |$authors"
  done
done
```

A side-branch copy is an "older snapshot" when its blob hash appears in `git log origin/main -- <path>`.
Line counts and commits will change as the branches move; the date at the top says when this was true.
