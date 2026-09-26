# Data inventory: MoTrPAC phosphoproteomics + GlyGen for our 471 proteins and 450 metabolites

## 1. Snapshot

| | |
|---|---|
| Question | Which phosphoproteomic (MoTrPAC, GlyGen) and glycosylation / other GlyGen (Homo sapiens) data exist for the 471 proteins and 450 metabolites in our graph, so the team can choose what to integrate into Neo4j? |
| Status | Run 2026-09-26 against GlyGen release 2.11.1 and MotrpacHumanPreSuspensionAnalysis v0.2.4. Counts only; nothing is integrated yet |
| Code | `glygen_protein_inventory.py` (GlyGen API, one call per protein), `glygen_motrpac_inventory.R` (MoTrPAC phospho, metabolite matching, coverage table), `export_phospho_features.R` (CSV of every MoTrPAC phospho feature ID) |
| Outputs (never in the repo) | `$HACK_OUT/inventory/`: `coverage_summary.csv`, `protein_inventory.csv` (one row per protein, every count), `metabolite_inventory.csv`, `glygen_protein_counts.csv`, `glygen_phosphosites.csv`, `glygen_glycosites.csv` (glycosylation site positions, used for the step 16c crosstalk), `motrpac_phospho_feature_ids.csv` (all 21,873 MoTrPAC phospho features with IDs and flags; `export_phospho_features.R`) |

## 2. How it was done

- **GlyGen** has 242 Homo sapiens datasets (117 protein, 54 proteoform, 69 glycan; 31 GB in total). Instead of
  downloading them, `glygen_protein_inventory.py` asks GlyGen's API for each of our proteins
  (`https://api.glygen.org/protein/detail/<canonical accession>/`), which returns every protein-level data
  section, and keeps only counts. All 471 genes map to a GlyGen canonical accession (via Entrez).
- **MoTrPAC phosphoproteomics**: muscle (0.5 / 4 / 24 h) and adipose (4 h only), both arms; no blood.
  Sites are mapped to our genes via `HUMAN_FEATURE_TO_GENE`; "responds" = adj. p < 0.05 in an
  exercise-vs-control (delta-delta) contrast; "known in GlyGen" = same canonical accession, position and
  residue as a GlyGen site (single-site features only, since multi-site IDs are separate features).
- **Metabolites** can meet GlyGen only as glycans (matched by ChEBI, PubChem compound ID or the
  stereo-free InChIKey block). GlyGen's human Rhea file lists enzymes only, so metabolite–enzyme links stay
  with our own step 5 (Rhea directly).

## 3. Results (2026-09-26)

### Proteins (471)

| Source | Data type | Proteins with data | Median per protein |
|---|---|---|---|
| MoTrPAC | phosphoproteomics, either tissue | 227 (48%) | 3 sites |
| MoTrPAC | muscle / adipose | 196 / 169 | 2 / 2 |
| MoTrPAC | a site responds (EE or RE vs control) | 88 (19%) | 1 |
| MoTrPAC | a measured site is already known in GlyGen | 168 (36%) | 2 |
| GlyGen | phosphorylation sites (UniProtKB, iPTMnet) | 385 (82%) | 4 |
| GlyGen | phosphosites with a known kinase | 121 (26%) | 2 |
| GlyGen | glycosylated at all (site known or protein-level) | 429 (91%) | 3 |
| GlyGen | glycosylation site with a known position | 307 (65%) | 4 |
| GlyGen | glycosylated, site unknown (protein-level; mostly O-GlcNAc Database) | 214 (45%) | — |
| GlyGen | reported sites with a glycan structure (GlyTouCan) | 281 (60%) | 4 sites, 2 glycans |
| GlyGen | N-linked / O-linked / O-GlcNAc sites | 179 / 267 / 202 | 3 / 3 / 2 |
| GlyGen | predicted / literature-mined sites | 135 / 42 | 1 / 1.5 |
| GlyGen | glycation sites | 1 | — |
| GlyGen | mutations / SNVs; mutagenesis | 422 (90%); 231 (49%) | 11; 4 |
| GlyGen | disease associations; biomarkers | 310 (66%); 23 (5%) | 4; 1 |
| GlyGen | PTM annotation; site annotation (active / binding) | 204 (43%); 349 (74%) | 1; 3 |
| GlyGen | expression normal tissue (Bgee); cancer (BioXpress) | 465; 471 | 13; 12 |
| GlyGen | pathways (Reactome, KEGG); reactions; enzyme (EC) | 471; 179; 163 | 3; 6; 1 |
| GlyGen | GO, function, structures (PDB), isoforms, orthologs, publications, cross-references | 470–471 each | — |

**Phosphorylation, site level:** MoTrPAC measures 909 unique phosphosite features on our proteins (739 in
muscle, 444 in adipose, 274 in both; 825 single-site, 787 confidently localised); 422 are sites GlyGen
already lists; 202 respond (86 in EE, 181 in RE; by tissue 192 muscle and 10 adipose measurements). GlyGen lists 2,699 phosphosites on our proteins, 398 with a
known kinase. Proteins with any phospho information: 398 (both sources 214, MoTrPAC only 13, GlyGen only
171, neither 73). 38 proteins have a responding MoTrPAC site and a GlyGen kinase annotation.

**Within the 286 proteins that have an edge in the gene network:** MoTrPAC phospho 134, GlyGen phospho
233, glycosylated 260.

### Metabolites (450)

| Data type | Metabolites |
|---|---|
| is a GlyGen glycan | 1 (glucose) |
| has a ChEBI ID (needed for most database matching) | 213 |
| Rhea enzyme among our 471 proteins (our step 5, for comparison) | 60 |

GlyGen is a glycan / glycoprotein resource: its small-molecule coverage is monosaccharides and glycans,
so it adds essentially nothing for our metabolites directly. Metabolites connect to GlyGen data only
through their enzymes (Rhea, step 5).

## 4. Read with care

- GlyGen's API detail does not include IntAct interactions for most proteins (use the IntAct file if
  interactions are needed); "reactions" come from Reactome and Rhea as GlyGen integrates them.
- Protein-level glycosylation records (no position) say a protein is glycosylated without a site; they
  cannot be joined to MoTrPAC sites or drawn as site nodes.
- MoTrPAC multi-site features (e.g. `P12345_S10sS12s`) are separate measurements from the single sites
  and were not matched to GlyGen.
- No glycoproteomics was measured in MoTrPAC; all glycosylation data are prior knowledge (GlyGen).
