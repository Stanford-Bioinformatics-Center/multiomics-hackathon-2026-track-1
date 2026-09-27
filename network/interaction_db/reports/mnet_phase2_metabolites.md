# mnet Phase 2 - metabolite mapping

## Correction A - PPI network is the >= 500 build

- ppi path: `/Users/monilgandhi/Desktop/string_database/data/processed/string_v12.0_human_uniprot_network_ge500.parquet`
- rows: 320,811 | nodes: 18,299 | min score: 500.0
- human Rhea enzymes: 4,140; of those already ge500 node ids: **4,016**

## Correction B - compound_type from Rhea classes

Rhea's compound class is now the source of truth (was ChEBI-formula heuristic).

| compound_type | Phase 1 (formula, human uniq ChEBI) | Phase 2 (Rhea class, all uniq ChEBI) |
|---|---|---|
| small molecule | 3644 | 12935 |
| generic | 624 | 1131 |
| polymer | 2 | 195 |
| other | - | 0 |
- participants with reactive-part/underlying ChEBI stored: 7,717

## Correction C - measured currency / high-degree metabolites

- measured metabolites flagged is_currency: 6
  ['ADP', 'AMP', 'ATP', 'FAD', 'GTP', 'NAD+']
- measured metabolites in > 150 human reactions: ['ADP', 'ATP', 'NAD+']
- (no filtering applied yet - that is Phase 4)

## Coverage by super_class

| super_class | total | mapped | in human Rhea | unmapped |
|---|---|---|---|---|
| (none) | 3 | 2 | 2 | 1 |
| Alkaloids | 14 | 14 | 5 | 0 |
| Benzenoids | 7 | 7 | 1 | 0 |
| Carbohydrates | 3 | 3 | 3 | 0 |
| Fatty Acyls | 83 | 68 | 52 | 15 |
| Glycerolipids | 1 | 1 | 0 | 0 |
| Glycerophospholipids | 7 | 7 | 1 | 0 |
| Nucleic acids | 28 | 28 | 25 | 0 |
| Organic acids | 66 | 66 | 58 | 0 |
| Organic nitrogen compounds | 6 | 5 | 4 | 1 |
| Organoheterocyclic compounds | 3 | 3 | 2 | 0 |
| Prenol Lipids | 2 | 2 | 1 | 0 |
| Sphingolipids | 5 | 5 | 4 | 0 |
| Sterol Lipids | 7 | 7 | 7 | 0 |
| **TOTAL** | 235 | 218 | 165 | 17 |

## Counts per mapping_method / normalization_method

- mapping_method: {'charge_normalized': 126, 'exact': 69, 'unmapped': 17, 'class_level': 16, 'name': 3, 'manual': 3, 'inchikey': 1}
- normalization_method: {'ph7_3_mapping': 188, 'none': 36, 'best_rhea_form': 11}

## Unmapped and low-confidence

- unmapped (17): ['Aminobutyric acid', 'C1-DeoxyCer 18:0;O/24:1', 'CAR 14:1', 'CAR 16:1', 'CAR 18:1', 'CAR 18:2', 'CAR 4:0;OH', 'CAR 5:0;OH', 'CAR 5:1', 'CAR 8:1', 'Hydroxydecanoic acid', 'Hydroxydodecanoic acid', 'Hydroxyoctanoic acid', 'Hydroxytetradecanoic acid', 'Nonadecenoic acid', 'Tetradecadienoic acid', 'Tridecylamine']
- low confidence (16): ['Docosapentaenoic acid', 'Docosatetraenoic acid', 'Docosatrienoic acid', 'Docosenoic acid', 'Dodecadienoic acid', 'Eicosadienoic acid', 'Eicosapentaenoic acid', 'Eicosatetraenoic acid', 'Eicosatrienoic acid', 'Eicosenoic acid', 'Hexadecenoic acid', 'Octadecadienoic acid', 'Octadecatetraenoic acid', 'Octadecatrienoic acid', 'Octadecenoic acid', 'Tetracosenoic acid']

## 10 charge-normalization examples (original -> Rhea)

| metabolite | original | -> Rhea | Rhea name |
|---|---|---|---|
| 2-Hydroxy-3-methylbutyric acid | CHEBI:60631 | CHEBI:744019 | (S)-2-hydroxy-3-methylbutyrate |
| 2-Oxoglutaric acid | CHEBI:30915 | CHEBI:16810 | 2-oxoglutarate(2-) |
| 3-Hydroxybenzoic acid | CHEBI:30764 | CHEBI:16193 | 3-hydroxybenzoate |
| 3-Hydroxybutyric acid | CHEBI:20067 | CHEBI:37054 | 3-hydroxybutyrate |
| 3-Indolepropionic acid | CHEBI:43580 | CHEBI:82916 | 3-(1H-indol-3-yl)propanoate |
| 3-Methyl-2-oxovaleric acid | CHEBI:35932 | CHEBI:28654 | 3-methyl-2-oxovalerate |
| 5alpha-Androstane-3alpha-ol-17-one sulfate | CHEBI:83037 | CHEBI:133003 | androsterone sulfate(1-) |
| 9S-HODE | CHEBI:34496 | CHEBI:77852 | 9(S)-HODE(1-) |
| ADP | CHEBI:16761 | CHEBI:456216 | ADP(3-) |
| AMP | CHEBI:16027 | CHEBI:456215 | adenosine 5'-monophosphate(2-) |

- obsolete ChEBI replacements applied: 0
- log lines (ties / merges / obsolete): 3 (see interim log)

