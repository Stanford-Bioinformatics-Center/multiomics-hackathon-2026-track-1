# mnet Phase 1 - downloads & interim caches

UniProt release: UniProt Knowledgebase Release 2026_03 | STRING v12.0

## Source files (data/raw; sizes; see MANIFEST.txt for SHA-256 + date)

| file | size (bytes) |
|---|---|
| rhea2uniprot_sprot.tsv | 8,835,509 |
| rhea-directions.tsv | 446,712 |
| rhea2xrefs.tsv | 1,196,556 |
| chebi_pH7_3_mapping.tsv | 3,099,222 |
| chebi.obo.gz | 47,162,742 |
| structures.tsv.gz | 92,755,407 |
| swisslipids_lipids.tsv.gz | 76,772,233 |
| HUMAN_9606_idmapping.dat.gz | 38,464,542 |

Rhea participants: fetched from the SPARQL endpoint (https://sparql.rhea-db.org/sparql); chose SPARQL over parsing rhea.rdf.gz with rdflib for speed/memory. ChEBI structures/InChIKey from the 2025 `flat_files/structures.tsv.gz` (obo lacks InChIKey; formula/charge/smiles come from the obo `chemrof:` property_values).

## Interim parquet row counts

| table | rows |
|---|---|
| rhea_enzymes | 400,643 |
| rhea_directions | 18,611 |
| rhea_xrefs | 36,188 |
| rhea_participants | 90,372 |
| chebi_names | 218,822 |
| chebi_structures | 218,822 |
| uniprot_genes | 171,713 |
| chebi_structures with inchikey | 189,896 |

## Rhea

- master reactions with participants: **18,611**
- reactions with >=1 enzyme (Swiss-Prot): **14,163**
- reactions with >=1 HUMAN enzyme: **6,151**
- human enzymes (UniProt accessions): **4,140**
- of those already STRING node ids as-is: **3,575** (Phase 4 crosswalk peek; STRING nodes=13,844)
- transport reactions (same compound both sides): 1,712

- TrEMBL: disabled (`mnet.use_trembl=false`). Swiss-Prot covers most human enzymes; enabling it would add TrEMBL-only human enzymes (rhea2uniprot_trembl.tsv.gz, not downloaded).

## ChEBI participants in human reactions

- unique ChEBI ids in human reactions: **4,978**
- by compound_type: {'generic': 644, 'polymer': 75, 'small molecule': 4259}

### Top 30 ChEBI ids by number of human reactions

| chebi_id | name | #human_reactions |
|---|---|---|
| CHEBI:15378 | hydron | 3249 |
| CHEBI:15377 | water | 2117 |
| CHEBI:15379 | dioxygen | 707 |
| CHEBI:57287 | coenzyme A(4-) | 688 |
| CHEBI:30616 | ATP(4-) | 480 |
| CHEBI:43474 | hydrogenphosphate | 433 |
| CHEBI:456216 | ADP(3-) | 376 |
| CHEBI:58349 | NADP(3-) | 364 |
| CHEBI:57783 | NADPH(4-) | 360 |
| CHEBI:57540 | NAD(1-) | 321 |
| CHEBI:57945 | NADH(2-) | 283 |
| CHEBI:58210 | FMN(3-) | 247 |
| CHEBI:58223 | UDP(3-) | 244 |
| CHEBI:57618 | FMNH2(2-) | 238 |
| CHEBI:29101 | sodium(1+) | 223 |
| CHEBI:16526 | carbon dioxide | 202 |
| CHEBI:33019 | diphosphate(3-) | 200 |
| CHEBI:59789 | S-adenosyl-L-methionine zwitterion | 193 |
| CHEBI:57856 | S-adenosyl-L-homocysteine zwitterion | 182 |
| CHEBI:456215 | adenosine 5'-monophosphate(2-) | 141 |
| CHEBI:16240 | hydrogen peroxide | 137 |
| CHEBI:16810 | 2-oxoglutarate(2-) | 132 |
| CHEBI:57288 | acetyl-CoA(4-) | 131 |
| CHEBI:28938 | ammonium | 114 |
| CHEBI:60377 | cytidine 5'-monophosphate(2-) | 111 |
| CHEBI:58343 | adenosine 3',5'-bismonophosphate(4-) | 103 |
| CHEBI:30031 | succinate(2-) | 102 |
| CHEBI:30823 | oleate | 95 |
| CHEBI:57387 | oleoyl-CoA(4-) | 87 |
| CHEBI:57692 | FAD(3-) | 85 |

## Currency metabolites

- rows written to curation/currency_metabolites.csv: 46 (charge forms present in Rhea participants)
- concepts NOT found in Rhea participants: ['oxidized flavodoxin', 'reduced flavodoxin']

## Metabolite CSV check

- rows: 450 | columns ok: True
- with ChEBI: 195

## SwissLipids

- status: ok (779,257 rows, cols=['Lipid ID', 'Level', 'Name', 'Abbreviation*']...)

