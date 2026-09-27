# Team 2-PAC · Exercise as Medicine: endurance vs resistance, in the context of type 2 diabetes

**Endurance and resistance exercise move disease-linked proteins back toward healthy in different tissues —
endurance in muscle, resistance in blood — read off one physically gated, exercise-weighted multi-omic network.**

![Figure 1](network/docs/figures/fig1_story.png)

| | |
|---|---|
| **Question** | Do endurance and resistance exercise differ in how they move the molecules disturbed in type 2 diabetes (and ageing, its main risk factor)? In which tissue? |
| **Data** | MoTrPAC human acute exercise (RNA, protein, metabolites; adipose, blood, muscle) · STRING v12 and Rhea physical interactions · seven published T2D / insulin-resistance / ageing proteomes incl. UK Biobank |
| **Method** | every edge = a **physical link** (hard: STRING ≥ 700 or Rhea) weighted by the **dot product of the two molecules' exercise responses** (soft), per arm — see the figure below and [`network/engine`](network/engine) |
| **Top 3 stories** | **1. Muscle · T2D:** endurance reverses the T2D muscle proteome (Öhman 2021), resistance does not (difference p 0.007); the T2D proteins form a connected subgraph (p 0.039) · **2. Blood · ageing:** resistance reverses the 293 plasma proteins of ageing (UK Biobank; 0.47, p < 0.001) · **3. Blood · future T2D:** resistance reverses the 297 plasma proteins linked to future T2D (UK Biobank; 0.38, p < 0.001) |
| **Reproduce** | `Rscript network/requirements.R`, then `bash network/run_all.sh` (see [Quick start](#quick-start) for the inputs) — ~10 minutes, 25 engine tests + 44 validation checks, 139 outputs byte-identical to the reference run |
| **Read more** | [`network/README.md`](network/README.md): full documentation (snapshot, question, workflow, setup, inputs/outputs, methods, validation, reuse) |

### How an edge is made: hard (physical) x soft (exercise) weights

![Method](network/docs/figures/fig_method_hard_soft_edges.png)

A physical database decides **whether** two molecules are connected (panel a: A and E respond alike but are not
physically linked, so they get no edge). Their normalised exercise responses decide **how strongly**, in each arm
(panels b, c: w(A,B) = Σ z_A · z_B). The result is one network per arm on identical edges, so every difference
between endurance and resistance is a difference in exercise response, not in wiring (panel d).

### Quick start

**1. Requirements.** Install **R ≥ 4.4** first (the reference run used 4.4.3). R provides `Rscript`, which runs
almost every step; without it `run_all.sh` stops at step 00 with `Rscript: command not found`.

- **macOS:** download the installer from [cloud.r-project.org/bin/macos](https://cloud.r-project.org/bin/macos/)
  (the **arm64** `.pkg` for Apple Silicon M1–M4, the **x86_64** `.pkg` for Intel; `uname -m` tells you which) and
  run it. It puts `Rscript` in `/usr/local/bin`, which is already on the `PATH`.
- **Windows / Linux:** follow [cloud.r-project.org](https://cloud.r-project.org) (on Windows, add R's `bin` folder to the `PATH`).
- **pandoc** (step 17i) comes with [RStudio](https://posit.co/download/rstudio-desktop/); installing RStudio is the easiest way to get it.
- Building packages from source (only with `--exact`) needs a compiler; on macOS: `xcode-select --install`.

Open a **new** terminal and check with `Rscript --version`. Everything else is listed in two files:

| File | Covers |
|---|---|
| [`network/requirements.R`](network/requirements.R) | the 26 R packages (CRAN, Bioconductor 3.20, the MoTrPAC package from GitHub) with the reference versions; installs them and TinyTeX, then checks |
| [`network/requirements.txt`](network/requirements.txt) | Python ≥ 3.9 (standard library only, nothing to pip install) and the system tools: pandoc (step 17i; ships with RStudio / Positron / Quarto), pdflatex (step 13), `curl` (step 20f) |

Exact versions of the reference run: [`network/ENVIRONMENT.md`](network/ENVIRONMENT.md).

**2. Inputs that are not in the repo.** The MoTrPAC results come with the R package, and the public resources
(Rhea, RefMet / PubChem, GlyGen, the published disease tables) download on the first run and are cached. Three
inputs have to be supplied; point the environment variables at them (defaults are the authors' own folders):

| Variable | What | Needed by |
|---|---|---|
| `MNET_DIR` | the team's mnet resource (`edges.parquet`, `nodes.csv`, `proteins_ptm.csv`, `phosphosites.csv`, `glycosites.csv`, `motrpac_feature_site_map.csv`): ask the team, or build it with [`network/interaction_db`](network/interaction_db) (its `data/output/` folder) | step 02 onwards |
| `DISEASE_SCORES` | `disease_scores.csv.gz`, Amar et al. 2024 disease sets | steps 17i, 18d, 18c, 18bc, 19, 20 |
| `UBAIDA` | Ubaida-Mohien 2019 muscle ageing proteome (`ubaida_mohien_ 2019_elife_stat.csv`, from the Amar et al. repository) | step 20 |

Outputs go outside the repo: tables to `HACK_OUT` (default `~/Desktop/output/hackathon-2026-track1/network`),
figures to `HACK_FIG` (default `~/Desktop/output/hackathon`). All variables: [`network/README.md`](network/README.md), section 4.

**3. Run.**

```bash
git clone https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1.git
cd multiomics-hackathon-2026-track-1
Rscript network/requirements.R             # install every R package (--check: only report; --exact: reference versions)
export MNET_DIR=/path/to/mnet DISEASE_SCORES=/path/to/disease_scores.csv.gz UBAIDA=/path/to/ubaida_mohien_\ 2019_elife_stat.csv
bash network/run_all.sh                    # every step, engine tests, validation, reproducibility manifest
bash network/run_all.sh 14 17i             # or a range of steps (labels in network/run_all.sh)
Rscript network/engine/run_tests.R         # just the engine: toy example by hand, bad inputs, exact reproduction
open ~/Desktop/output/hackathon/17_interactive/17a_joint_network.html   # the interactive network
```

**If it fails.** The run stops at the first failing step and prints the end of its log (`$HACK_OUT/logs/<step>.log`).
- `Rscript: command not found` — R is not installed or not on the `PATH` (step 1).
- `there is no package called ...` — run `Rscript network/requirements.R --check` to see what is missing.
- a missing file under `MNET_DIR`, `DISEASE_SCORES` or `UBAIDA` — step 2.

> Molecular overlap or signature reversal does not demonstrate clinical benefit (see the brief below): our results
> are direction matches between one acute bout in healthy adults and disease / age signatures.

---

## Track brief (from the organisers)

### Stanford Multi-omics Hackathon 2026 Track 1

## Exercise as Medicine

*What relationships connect exercise-responsive biology with human disease?*

### Challenge

Integrate MoTrPAC results with disease genes, pathways, variants, biomarkers, or expression signatures to generate transparent, testable hypotheses.

### Data

Human acute-exercise results together with resources such as [Open Targets](https://platform.opentargets.org/), [GWAS Catalog](https://www.ebi.ac.uk/gwas/), [DisGeNET](https://www.disgenet.com/), [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [OMIM](https://omim.org/), and approved [CFDE](https://commonfund.nih.gov/dataecosystem) resources.

### Potential Outputs

- An exercise–disease evidence map
- Prioritization score
- Reproducible enrichment workflow
- Evidence explorer

> [!IMPORTANT]
> Molecular overlap or signature reversal does not demonstrate clinical benefit or justify medical recommendations. Team should select a single disease prior to the hackathon.
