# mnet Phase 0 - scaffold

Adds the metabolite-network extension skeleton next to `string_network`. **No
data is processed in this phase.** The existing STRING pipeline is untouched.

## Created

- **Package `src/mnet/`**: `__init__`, `__main__`, `cli`, `download`, `rhea`,
  `chebi`, `metabolites`, `lipids`, `proteins`, `network`, `assemble`, `qc`.
  CLI: `python -m mnet {download,metabolites,lipids,rhea,proteins,network,assemble,qc}`
  (each is a Phase-0 stub that reports "not implemented yet" and returns 0).
- **Dirs**: `data/input/`, `data/interim/`, `data/output/` (with `.gitkeep`),
  `tests/mnet/fixtures/`. `.gitignore` updated to keep the placeholders while
  still ignoring downloaded/generated data.
- **Curation (committed)**: `curation/edge_scores.csv` (placeholder scores),
  `curation/currency_metabolites.csv` (headers only - filled in Phase 1),
  `curation/manual_metabolite_mappings.csv` (headers only).
- **Config**: new top-level `mnet:` section; new `urls:` keys `rhea_rhea2uniprot`,
  `rhea_directions`, `rhea_rdf`, `rhea_sparql`, `chebi_obo`, `swisslipids_lipids`.
  No existing keys renamed or removed.
- **Makefile**: `mnet-download`, `mnet-metabolites`, `mnet-lipids`, `mnet-rhea`,
  `mnet-assemble`, `mnet-qc`, `mnet-test`, `mnet-all`.
- **requirements.txt**: added `rdflib`, `SPARQLWrapper`, `pygoslin`, `networkx`.
- **README**: "Metabolite extension (`mnet`)" section with the phase table.
- **Tests**: `tests/mnet/test_scaffold.py` (imports, CLI wiring, config section,
  URL keys present).

## External URL verification (GET status)

| key | status | type | size | note |
|---|---|---|---|---|
| rhea_rhea2uniprot | 200 | tsv | 8.8 MB | OK |
| rhea_directions | 200 | tsv | 0.45 MB | OK |
| rhea_rdf | 200 | gzip | 8.7 MB | OK |
| chebi_obo | 200 | gzip | 47 MB | OK |
| rhea_sparql | 200 | html | - | endpoint reachable |
| swisslipids_lipids | 200 | **html, empty body** | 0 | **needs attention** |

**SwissLipids**: the documented API URL
`https://www.swisslipids.org/api/file.php?cast=normal&file=lipids.tsv` returns
HTTP 200 with an **empty body** (tried `cast`/no-`cast` and an alternate path).
It is only needed in Phase 3 (complex-lipid parsing), and `pygoslin` can parse
lipid shorthand without it, so this does not block Phase 0/1. **Flagged for the
user** - I did not guess a replacement URL.

## Acceptance checks

- `make build` -> byte-identical `..._ge700.parquet` (SHA-256 unchanged).
- `make test` -> **18 passed** (13 existing unchanged + 5 new mnet scaffold).
- `make mnet-test` -> 5 passed.

## Blocker for Phase 1

`data/input/metabolite_chebi_ids.csv` is **not present yet** (no `data/input/*.csv`).
Phase 1 (metabolite parsing) needs it. Please copy it in before Phase 1.
