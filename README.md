# Team 2-PAC · Exercise as Medicine: endurance vs resistance, in the context of type 2 diabetes

**Endurance and resistance exercise move disease-linked proteins back toward healthy in different tissues —
endurance in muscle, resistance in blood — read off one physically gated, exercise-weighted multi-omic network.**

![Figure 1](network/docs/figures/fig1_story.png)

| | |
|---|---|
| **Question** | Do endurance and resistance exercise differ in how they move the molecules disturbed in type 2 diabetes (and ageing, its main risk factor)? In which tissue? |
| **Data** | MoTrPAC human acute exercise (RNA, protein, metabolites; adipose, blood, muscle) · STRING v12 and Rhea physical interactions · seven published T2D / insulin-resistance / ageing proteomes incl. UK Biobank |
| **Method** | every edge = a **physical link** (hard: STRING ≥ 700 or Rhea) weighted by the **dot product of the two molecules' exercise responses** (soft), per arm — see the figure below and [`network/engine`](network/engine) |
| **Finding** | muscle: endurance reverses the T2D proteome, resistance does not (difference p 0.007), and the T2D proteins form a connected subgraph (p 0.039); blood: resistance reverses the plasma proteins of ageing and of future T2D (UK Biobank, p < 0.001) |
| **Reproduce** | `bash network/run_all.sh` — ~10 minutes, 25 engine tests + 44 validation checks, 136 outputs byte-identical to the reference run |
| **Read more** | [`network/README.md`](network/README.md): full documentation (snapshot, question, workflow, setup, inputs/outputs, methods, validation, reuse) |

### How an edge is made: hard (physical) x soft (exercise) weights

![Method](network/docs/figures/fig_method_hard_soft_edges.png)

A physical database decides **whether** two molecules are connected (panel a: A and E respond alike but are not
physically linked, so they get no edge). Their normalised exercise responses decide **how strongly**, in each arm
(panels b, c: w(A,B) = Σ z_A · z_B). The result is one network per arm on identical edges, so every difference
between endurance and resistance is a difference in exercise response, not in wiring (panel d).

### Quick start

```bash
git clone https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-1.git
cd multiomics-hackathon-2026-track-1
# R 4.4 + packages: see network/README.md section 4 (Setup)
bash network/run_all.sh                    # every step, engine tests, validation, reproducibility manifest
Rscript network/engine/run_tests.R         # just the engine: toy example by hand, bad inputs, exact reproduction
open ~/Desktop/output/hackathon/17_interactive/17a_joint_network.html   # the interactive network
```

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
