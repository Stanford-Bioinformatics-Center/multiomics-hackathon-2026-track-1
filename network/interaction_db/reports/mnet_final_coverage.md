# mnet final coverage

- nodes: 22,836 {'protein': 18378, 'metabolite': 4438, 'lipid_class': 20}
- edges: 338,238 {'ppi': 320811, 'catalysis': 15066, 'transport': 2054, 'lipid_is_a': 307}
- components: 72 | largest component: 22,722 nodes

## Measured metabolites by super_class

| super_class | total | >=1 edge | protein-reachable | in largest comp |
|---|---|---|---|---|
| (lipid/none) | 230 | 229 | 229 | 229 |
| Alkaloids | 14 | 5 | 5 | 5 |
| Benzenoids | 7 | 1 | 1 | 1 |
| Carbohydrates | 3 | 3 | 3 | 3 |
| Fatty Acyls | 68 | 67 | 66 | 67 |
| Glycerolipids | 1 | 1 | 1 | 1 |
| Glycerophospholipids | 7 | 7 | 7 | 7 |
| Nucleic acids | 28 | 25 | 25 | 25 |
| Organic acids | 66 | 58 | 58 | 58 |
| Organic nitrogen compounds | 6 | 4 | 4 | 4 |
| Organoheterocyclic compounds | 3 | 2 | 2 | 2 |
| Prenol Lipids | 2 | 1 | 1 | 1 |
| Sphingolipids | 5 | 5 | 5 | 5 |
| Sterol Lipids | 7 | 7 | 7 | 7 |
| **TOTAL** | 447 | 415 | 414 | 415 |

## Proteins
- protein nodes: 18,378; Rhea enzymes new (in_string=False): 79

## Lipids by lipid_class (measured)
{'TG': 58, 'PC': 40, 'SM': 31, 'PE': 21, 'PC-ether': 14, 'Cer': 13, 'LPC': 10, 'PE-ether': 10, 'FA-long-chain': 10, 'DG': 8, 'LPE': 4, 'MG': 3, 'LPC-ether': 2, 'LPE-ether': 1, 'FA-very-long-chain': 1, 'FA-medium-chain': 1}

## Top 10 hubs per node type
- protein: P04637(753), P17612(751), P22694(741), P22612(735), P60709(646), P01375(641), P62873(639), P31749(596), P00533(593), P05231(589)
- metabolite: CHEBI:57379(118), CHEBI:29985(106), CHEBI:7896(104), CHEBI:16810(98), REFMET:RM0152034(87), CHEBI:57288(85), CHEBI:30031(84), CHEBI:30823(83), CHEBI:30616(79), CHEBI:456215(74)
- lipid_class: LIPIDCLASS:FA-long-chain(151), LIPIDCLASS:PC(137), LIPIDCLASS:PC-ether(113), LIPIDCLASS:FA-medium-chain(97), LIPIDCLASS:FA-very-long-chain(97), LIPIDCLASS:DG(93), LIPIDCLASS:Cer(79), LIPIDCLASS:PE(77), LIPIDCLASS:LPC(74), LIPIDCLASS:TG(70)

## Suggested targets for degree-0 metabolites (20) - NOT auto-applied
- CHEBI:133185 (Acisoga) -> CHEBI:233601 (S-palmitoyl-N-acetylcysteamine, degree 2)
- CHEBI:16020 (1-Methyladenosine) -> CHEBI:21891 (N(6)-methyladenosine, degree 1)
- CHEBI:16193 (3-Hydroxybenzoic acid) -> CHEBI:17879 (4-hydroxybenzoate, degree 2)
- CHEBI:17563 (Phthalic acid) -> CHEBI:30031 (succinate(2-), degree 84)
- CHEBI:176478 (N-Acetylisoputreanine) -> CHEBI:35924 (peroxol, degree 10)
- CHEBI:32979 (Phenyllactic acid) -> CHEBI:11009 ((R)-3-phenyllactate, degree 1)
- CHEBI:34698 (Diethyl phthalate) -> CHEBI:76478 (dibutyrin, degree 5)
- CHEBI:35280 (Proline betaine) -> CHEBI:17750 (glycine betaine, degree 7)
- CHEBI:38905 (Tributylamine) -> CHEBI:16269 (N,N-dimethylaniline, degree 3)
- CHEBI:46468 (alpha-Tocoquinone) -> CHEBI:142491 (R-H, degree 25)
- CHEBI:58071 (Monoethylhexyl phthalic acid) -> CHEBI:16150 (benzoate, degree 8)
- CHEBI:58310 (alpha-N-Phenylacetylglutamine) -> CHEBI:132071 (N-arachidonoyl-L-alaninate, degree 3)
- CHEBI:5832 (Lenticin) -> CHEBI:17750 (glycine betaine, degree 7)
- CHEBI:61887 (N-Acetylglycine) -> CHEBI:29746 (glycocholate, degree 14)
- CHEBI:64755 (Edetic acid) -> CHEBI:76942 (beta-citrylglutamate(4-), degree 2)
- CHEBI:72991 (Imidazolepropionic acid) -> CHEBI:30089 (acetate, degree 61)
- CHEBI:74076 (Ile-Pro) -> CHEBI:90039 (cyclo(L-His-L-Pro), degree 3)
- CHEBI:744019 (2-Hydroxy-3-methylbutyric acid) -> CHEBI:134178 (all-trans-4-hydroxyretinoate, degree 10)
- CHEBI:82916 (3-Indolepropionic acid) -> CHEBI:30089 (acetate, degree 61)
- CHEBI:89825 (5-Hydroxytryptophol) -> CHEBI:55515 (N-acetyltryptamine, degree 1)

## Measured metabolites with NO protein connection (33)
- CHEBI:133185 (Acisoga): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:157782 (S-Methylcysteine S-oxide): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:16020 (1-Methyladenosine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:16193 (3-Hydroxybenzoic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:165212 (bis(2-Ethylhexyl)phthalic acid): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:166577 (2-Acetylpyrrolidine): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:17563 (Phthalic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:176478 (N-Acetylisoputreanine): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:18123 (Trigonelline): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:222828 (Monoethylglycinexylidide): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:25858 (Paraxanthine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:27732 (Caffeine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:28821 (Piperine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:32979 (Phenyllactic acid): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:34595 (Benzyl butyl phthalate): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:34698 (Diethyl phthalate): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:35280 (Proline betaine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:38905 (Tributylamine): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:46468 (alpha-Tocoquinone): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:58071 (Monoethylhexyl phthalic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:58310 (alpha-N-Phenylacetylglutamine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:5832 (Lenticin): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:61887 (N-Acetylglycine): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:64755 (Edetic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:72991 (Imidazolepropionic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:74076 (Ile-Pro): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:744019 (2-Hydroxy-3-methylbutyric acid): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:77761 (2-Piperidinone): not in human Rhea (diet/xenobiotic or no human enzyme)
- CHEBI:82916 (3-Indolepropionic acid): in Rhea but only as currency/cofactor (edges dropped)
- CHEBI:89825 (5-Hydroxytryptophol): not in human Rhea (diet/xenobiotic or no human enzyme)
- MEAS:C1_DeoxyCer_18_0_O_24_1 (C1-DeoxyCer 18:0;O/24:1): unmapped
- REFMET:RM0073571 (Aminobutyric acid): unmapped
- REFMET:RM0149929 (Tridecylamine): unmapped

## Evidence for including unmeasured metabolites
- measured metabolites in largest component WITH unmeasured: 92.8%
- WITHOUT unmeasured metabolites: 92.8%

