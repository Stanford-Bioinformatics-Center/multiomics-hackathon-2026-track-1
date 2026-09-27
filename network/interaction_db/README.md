# Interaction network database: STRING + Rhea metabolites + PTM annotations

This pipeline builds the human interaction network used by the Track 1 network analysis:

- **Protein–protein:** STRING v12.0 (human, taxon 9606), mapped to UniProt IDs. The default network is **≥ 500**.
- **Protein–metabolite:** Rhea enzyme reactions (catalysis and transport), human enzymes only.
- **Lipid → lipid class:** links measured lipids to their class, which carries the Rhea enzymes.
- **Protein PTM annotations:** phosphosites (UniProt + OmniPath, with kinases), glycosites (UniProt), and a bridge from MoTrPAC phospho `feature_id`s to these sites.

All reference data is public and downloaded by the pipeline. **No data is stored in this repo** (`data/` is gitignored).

---

## How to run (full pipeline)

Run everything from this folder (`network/interaction_db/`). You need Python 3.10+, `make`, internet access (about 1 GB of downloads) and about 2 GB of free disk space.

```bash
cd network/interaction_db

# 1. Environment
make setup                 # create .venv, install requirements + the local package

# 2. Inputs (NOT in the repo; ask the team for the current copies) -> data/input/
#    - metabolite_chebi_ids.csv         measured metabolite list (required for the metabolite layer)
#    - motrpac_phospho_feature_ids.csv  MoTrPAC phospho feature_ids (optional; only for the MoTrPAC bridge)

# 3. Protein–protein network (STRING)
make data                  # download STRING + UniProt sources into data/raw/
make build-500             # STRING >= 500 network; the metabolite layer attaches to this one

# 4. Metabolite layer (Rhea)
make mnet-all              # download Rhea/ChEBI/SwissLipids -> metabolites -> lipids -> Rhea edges -> assemble
make mnet-qc               # quality checks (expect: all QC checks passed)

# 5. PTM annotations (phosphosites, glycosites, MoTrPAC bridge)
make mnet-ptm

# 6. Tests
make test
```

To run the metabolite layer one step at a time instead of `make mnet-all`:
`make mnet-download` → `mnet-metabolites` → `mnet-lipids` → `mnet-rhea` → `mnet-assemble` → `mnet-qc`.

**Input file columns**

| File | Required columns |
|---|---|
| `metabolite_chebi_ids.csv` | metabolite, refmet_name, refmet_id, super_class, main_class, pubchem_cid, inchi_key, chebi_id, chebi_all, chebi_method |
| `motrpac_phospho_feature_ids.csv` | feature_id, flanking_sequence, in_muscle, in_adipose |

If the MoTrPAC file is missing, `make mnet-ptm` still runs and just skips the bridge.

---

## Outputs (`data/output/`)

| File | One row per | Use it for |
|---|---|---|
| `edges_3col.parquet` | undirected pair (`node1, node2, score`) | **the watershed input**, in the same layout as the STRING network file |
| `edges.parquet` / `edges.csv` | pair × edge type | full detail: `edge_type` (ppi / catalysis / transport / lipid_is_a), source, Rhea evidence, confidence |
| `nodes.csv` | node | `node_type` (protein / metabolite / lipid_class), label, `is_measured`, ChEBI/RefMet/UniProt IDs, degree |
| `proteins_ptm.csv` | protein | phospho and glyco site counts and lists, upstream kinases, `is_kinase`, number of MoTrPAC sites |
| `phosphosites.csv` | phosphosite | site_id (`P46020_S758`), kinases, phosphatases, evidence, 15-residue sequence window |
| `kinase_substrate.csv` | kinase → site | kinase–substrate links (the input for kinase activity analysis, e.g. KSEA) |
| `glycosites.csv` | glycosite | type (N-linked / O-linked / O-GlcNAc / C-linked), evidence |
| `motrpac_feature_site_map.csv` | MoTrPAC feature × site | MoTrPAC `feature_id` → canonical site_id, tissue, network protein, known kinase |

The STRING networks themselves are written to `data/processed/string_v12.0_human_uniprot_network_ge{threshold}.parquet`.
Coverage reports are in `reports/`: start with `mnet_final_coverage.md` (metabolites) and `mnet_phase6_ptm.md` (PTMs).

**Important notes**
- **Scores:** only protein–protein scores are real STRING scores (see *The score rule* below). Scores for catalysis (900), transport (800) and lipid_is_a (500) are **placeholders**, set in `curation/edge_scores.csv`.
- **Node IDs:** proteins keep their STRING network ID (mostly UniProt accessions). Metabolites use `CHEBI:`, `REFMET:` or `MEAS:` prefixes, and lipid classes use `LIPIDCLASS:`. Use `node_type` rather than working out the type from the ID.
- **`degree` in `nodes.csv`** counts rows in `edges.csv`, so a pair with both catalysis and transport links counts twice.
- **Metabolite–metabolite links are deliberately not included** (only lipid → class).

---

## STRING protein–protein network

### The score rule

The scores are **not** the raw STRING v12 combined scores. About 31% of them are non-integer (multiples of 0.1). The rule:

> **`combined_score = physical_score` if the pair is in STRING's *physical* subnetwork, otherwise `0.9 × full_score`, rounded to 1 decimal; the threshold is applied last.**

Validated on the ≥ 700 network: the rule matches 99.92% of scores exactly, and the rebuilt network shares 98.3% of edges (edge Jaccard 0.983, node Jaccard 0.992). This is implemented as `score_mode: physical_else_scaled` in `build.py`; the raw behaviour is still available as `score_mode: raw`.

> **Caveat:** the rule was derived and validated only on edges scoring **≥ 700**. Edges between 500 and 700 follow the same rule, but there is no reference to check them against.

### Mapping (ENSP → UniProt)

STRING IDs (`9606.ENSP…`) are mapped with a fallback chain (`mapping_fallback: true`):

1. **UniProt accession:** from UniProt `HUMAN_9606_idmapping` (STRING rows), widened with the STRING `UniProt_AC` and `Ensembl_UniProt` aliases. A reviewed (Swiss-Prot) accession is preferred, and exactly one accession is chosen per ENSP.
2. **STRING `preferred_name`** (gene symbol), from `protein.info`.
3. **The raw STRING ID** (`9606.ENSP…`).

Isoform suffixes (`-\d+`) are stripped. Every protein in `protein.info` is included, so no STRING protein is dropped for lack of a UniProt accession.

### Thresholds

```bash
make build-500               # ..._ge500.parquet (default for the metabolite layer)
make build                   # ..._ge700.parquet (THRESHOLD defaults to 700)
make build THRESHOLD=400     # any threshold
```

The threshold is applied **last**, so a ≥ 500 build filtered to ≥ 700 is identical to the ≥ 700 build.

**Optional, if you have the reference file:** put it at `data/reference/reference_network_ge700.parquet`, then run `make compare` (diff against the reference), `make inspect` (profile it) or `make tune` (search over build settings). These reference-comparison reports are not included in this repo.

**UniProt release differences.** The pipeline always downloads the current UniProt release, so a few proteins can map to different IDs over time. For example, `CGB1` and `DGCR6` now map to `K7ELM3` and `Q14129`. The score rule itself doesn't depend on the release.

---

## Metabolite layer (`src/mnet/`)

Adds metabolites to the ≥ 500 STRING network (`mnet.ppi_threshold: 500`) without changing it. The protein–protein edges in the output are identical to the STRING file.

1. **Metabolites** (`mnet-metabolites`): measured ChEBI IDs are converted to the form Rhea uses (the charged form at pH 7.3). Ambiguous names get all their candidate IDs (e.g. Leucine/Isoleucine). `curation/manual_metabolite_mappings.csv` overrides the automatic mapping.
2. **Lipids** (`mnet-lipids`): names are parsed with Goslin and matched to SwissLipids. Each lipid is linked to a lipid class (PC, PE, TG, Cer, the ether classes, fatty acid chain-length classes, acylcarnitines…). Classes map to the general Rhea entries listed in `curation/lipid_class_to_chebi.csv`.
3. **Rhea edges** (`mnet-rhea`): human enzymes only, with enzyme IDs matched to the STRING node IDs. Reaction types are catalysis and transport; protein-type participants (e.g. `L-seryl-[protein]`) are excluded.
   - **Common metabolites** (H₂O, H⁺, ATP/ADP, NAD(P)H, CoA, ions…) are removed from enzyme links.
   - For **measured** common metabolites, a cofactor rule keeps the link only when the reaction actually makes or breaks down the metabolite (e.g. ATP synthase, adenylate kinase, NAD⁺ synthesis). Kinases, ATPases and GTPases don't count. The rules are in `curation/currency_*.csv`.
4. **Assemble + QC** (`mnet-assemble`, `mnet-qc`): writes `edges.*` and `nodes.csv`, then checks that the protein–protein edges match STRING exactly, that there are no duplicates or self-loops, that every edge node exists in `nodes.csv`, and that every input metabolite is represented.

## PTM annotations (`make mnet-ptm`, `src/mnet/ptm.py`)

- **Phosphosites:** the union of UniProt (reviewed, "Modified residue") and OmniPath enzyme–substrate data (filtered to sources that allow commercial use), with kinases and phosphatases. PhosphoSitePlus is **off** by default (non-commercial licence; `mnet.use_phosphositeplus`).
- **Glycosites:** UniProt "Glycosylation" features, with evidence marked experimental or predicted.
- **Site ID:** `<canonical UniProt>_<residue><position>` (e.g. `P46020_S758`).
- **MoTrPAC bridge:** splits multi-site features into individual sites and maps isoform sites to the main protein using the 15-residue sequence window. Sites that can't be matched uniquely are flagged `isoform_only`, and residue mismatches are flagged `residue_mismatch`; neither is dropped.

---

## Configuration

All settings are in `config.yaml` and can be overridden on the command line (`python -m string_network build --<option> …`). No paths or URLs are hard-coded anywhere else.

| option | default | meaning |
|---|---|---|
| `threshold` | 700 | STRING build threshold (`make build-500` overrides it) |
| `threshold_when` | `after` | apply the threshold before or after mapping/aggregation |
| `agg` | `max` | how duplicate undirected pairs are combined |
| `decimals` | 1 | round the final score |
| `score_mode` | `physical_else_scaled` | `raw` or the score rule above |
| `scale_factor` | 0.9 | multiplier for non-physical edges |
| `mapping_fallback` | `true` | UniProt → preferred_name → raw ID chain |
| `mapping_alias_sources` | `[UniProt_AC, Ensembl_UniProt]` | extra alias sources for UniProt coverage |
| `mnet.ppi_threshold` | 500 | which STRING build the metabolite layer attaches to |
| `mnet.measured_currency_policy` | `cofactor_rule` | how measured common metabolites are handled (`cofactor_rule` / `drop` / `keep_all`) |
| `mnet.use_phosphositeplus` | `false` | include PhosphoSitePlus (non-commercial licence) |

## Data sources and licences

Downloaded automatically. File sizes, SHA-256 hashes, download dates and release versions are logged in `data/raw/MANIFEST.txt`.

| source | used for | licence |
|---|---|---|
| STRING v12.0 (links, physical links, info, aliases) | protein–protein network | CC BY 4.0 |
| UniProt (idmapping, reviewed list, PTM features) | ID mapping, phospho/glyco sites | CC BY 4.0 |
| Rhea | enzyme–metabolite reactions | CC BY 4.0 |
| ChEBI | metabolite IDs, charge forms, hierarchy | CC BY 4.0 |
| SwissLipids | lipid classes | CC BY 4.0 |
| OmniPath (commercially licensed resources only) | kinase–substrate links | per resource |

Cite: STRING (Szklarczyk et al., NAR 2023), UniProt (The UniProt Consortium, NAR 2023), Rhea (Bansal et al., NAR 2022), ChEBI (Hastings et al., NAR 2016), SwissLipids (Aimo et al., Bioinformatics 2015), OmniPath (Türei et al., Mol Syst Biol 2021).

## Project layout

```
config.yaml            all settings
Makefile               setup / data / build / mnet-* / test targets
src/string_network/    STRING download, mapping, score rule, build, compare
src/mnet/              metabolites, lipids, Rhea, assembly, QC, PTM annotations
scripts/               reference profiling, tuning, score investigation (need the reference file)
curation/              editable mapping and scoring tables
tests/                 tests with small synthetic fixture files
reports/               coverage and mapping reports
data/                  input / raw / interim / processed / output (gitignored)
```
