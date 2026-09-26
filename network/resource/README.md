# Feature lists: 471 proteins and 450 metabolites

Reference lists of the features in the endurance-vs-resistance networks, for looking them up in other
datasets. Names and identifiers only; no measurements.

**The tables are not stored in this repository** (to keep the data as private as possible); only the
script that makes them is. To generate them, install the MoTrPAC R package (see `../README.md`, section
4), run pipeline steps 1, 1b, 1c and 2, then:

```bash
Rscript network/resource/export_feature_lists.R
```

The two CSVs below are written to `$HACK_RES` (default
`~/Desktop/output/hackathon-2026-track1/network/resource`), outside the repo. Share them directly with
teammates who need them rather than committing them.

## `proteins_471.csv`

The 471 genes measured as both RNA and protein in adipose, blood and muscle (MoTrPAC human
pre-suspension study; blood protein = OLINK panel, adipose and muscle protein = MS proteomics).

| Column | Meaning |
|---|---|
| `gene_symbol` | HGNC gene symbol |
| `entrez_gene` | NCBI Entrez gene ID (the ID used to match RNA and protein in the pipeline) |
| `uniprot` | UniProt accession(s) of the measured protein; several separated by `;` |
| `ensembl_gene` | Ensembl gene ID(s) |
| `in_string` | whether the protein appears in the curated STRING file (combined_score ≥ 700) |
| `network_degree` | number of STRING partners among the 471 (0 = isolated) |

Coverage: UniProt 471, Ensembl 471, in STRING 444.

## `metabolites_450.csv`

The 450 metabolites measured in adipose, blood (plasma) and muscle, matched across tissues by exact
RefMet name.

| Column | Meaning |
|---|---|
| `metabolite` | RefMet name, exactly as in the MoTrPAC data (punctuation in lipid names is meaningful: `/` = chain positions known, `_` = unknown) |
| `refmet_id` | RefMet identifier (Metabolomics Workbench) |
| `super_class`, `main_class` | RefMet chemical classes (14 super classes; 50 main classes) |
| `chebi_id` | one ChEBI ID (lowest number); `chebi_all` lists all, separated by `;` |
| `pubchem_cid` | PubChem compound ID |
| `inchi_key` | InChIKey (structure fingerprint) |
| `kegg_id` | KEGG compound ID, where the MoTrPAC package provides one |
| `platform_adipose`, `platform_blood`, `platform_muscle` | assay platform that measured it in each tissue (e.g. `metab-u-lrppos` = untargeted lipid reversed-phase, positive mode; `metab-t-*` = targeted panel) |

Coverage: RefMet ID 447, ChEBI 213, PubChem 208, InChIKey 208, KEGG 156. Blanks are mostly lipid
species defined without chain positions, which have no single structure; match those by RefMet name
or class. IDs were looked up on 2026-09-26 (see step 1c in `../README.md`).
