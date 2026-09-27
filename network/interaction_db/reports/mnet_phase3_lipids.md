# mnet Phase 3 - lipids

- rows in scope: 310 | parse success: 274 (88%) | parse fail: 36
- SwissLipids hits: 201 | exact_chebi assigned: 65
- lipid_is_a edges: 309 | class nodes: 20

## Counts by level
{'species': 137, 'sn_position': 52, 'molecular_species': 40, nan: 36, 'full_structure': 29, 'complete_structure': 12, 'structure_defined': 4}

## Counts by class (top 20)

| class | n |
|---|---|
| TG | 59 |
| FA-long-chain | 52 |
| PC | 44 |
| SM | 31 |
| PE | 23 |
| acylcarnitine | 19 |
| PC-ether | 14 |
| Cer | 13 |
| LPC | 11 |
| PE-ether | 10 |
| DG | 8 |
| FA-medium-chain | 5 |
| SPB | 4 |
| LPE | 4 |
| FA-short-chain | 4 |
| MG | 3 |
| LPC-ether | 2 |
| FA-very-long-chain | 2 |
| LPE-ether | 1 |
| SPBP | 1 |

## Class -> ChEBI (with human Rhea reactions)

| class_node | n_chebi | n_human_rhea_reactions |
|---|---|---|
| LIPIDCLASS:SPB | 12 | 144 |
| LIPIDCLASS:FA-long-chain | 2 | 68 |
| LIPIDCLASS:LPC | 8 | 65 |
| LIPIDCLASS:FA-very-long-chain | 2 | 64 |
| LIPIDCLASS:Cer | 14 | 64 |
| LIPIDCLASS:FA-medium-chain | 2 | 63 |
| LIPIDCLASS:LPC-ether | 6 | 63 |
| LIPIDCLASS:PC-ether | 5 | 56 |
| LIPIDCLASS:PC | 3 | 49 |
| LIPIDCLASS:DG | 7 | 46 |
| LIPIDCLASS:LPE | 14 | 40 |
| LIPIDCLASS:LPE-ether | 4 | 34 |
| LIPIDCLASS:PE | 4 | 32 |
| LIPIDCLASS:PE-ether | 4 | 32 |
| LIPIDCLASS:SPBP | 10 | 21 |
| LIPIDCLASS:TG | 2 | 14 |
| LIPIDCLASS:SM | 4 | 11 |
| LIPIDCLASS:MG | 2 | 3 |
| LIPIDCLASS:acylcarnitine | 1 | 2 |
| LIPIDCLASS:FA-short-chain | 1 | 0 |

## Classes with 0 human Rhea reactions (1)
['LIPIDCLASS:FA-short-chain']

## Parse failures (36)
['2-Hydroxy-3-methylbutyric acid', '2-Hydroxybutyric acid', '2-Hydroxyglutaric acid', '3-Hydroxybutyric acid', '3-Hydroxyphenyl valeric acid', 'Arachidic acid', 'Azelaic acid', 'Behenic acid', 'Capric acid', 'Glutaric acid', 'Hydroxydecanoic acid', 'Hydroxydodecanoic acid', 'Hydroxyhexadecanoic acid', 'Hydroxyoctanoic acid', 'Hydroxytetradecanoic acid', 'Lauric acid', 'Lignoceric acid', 'Linoleoyl-EA', 'Maleic acid', 'Margaric acid', 'Myristic acid', 'Myristoleic acid', 'NAGly 11:0', 'NAGly 12:0', 'Oleamide', 'Oleoyl-EA', 'PGE2 ethanolamide', 'Palmitoleic acid', 'Palmitoleoyl-EA', 'Palmitoyl-EA', 'Pentadecylic acid', 'Sebacic acid', 'Stearic acid', 'Suberic acid', 'Tricosylic acid', 'gamma-Aminobutyric acid']

