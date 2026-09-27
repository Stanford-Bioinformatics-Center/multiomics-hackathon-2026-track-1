# mnet Phase 6 - PTM annotation

- protein nodes: 18,378
- with >=1 phosphosite: 8,334 (45.3%); >=1 glycosite: 4,189 (22.8%)
- protein nodes with no UniProt accession: 12

## Phase 6.1 before -> after
- isoform sites mapped: 0 / 5,796 (0%) -> 5,047 mapped, 682 isoform_only (88% of isoforms)
- phosphosites total: 40,077 (UniProt only) -> 47,409 (UniProt-only=33,085 / OmniPath-only=7,344 / both=6,980)
- MoTrPAC % in phosphosites.csv: ~36.7% -> 43.1%; % with known kinase: ~9.0% -> 11.4%
- proteins with >=1 phosphosite: -> 8,334 (45.3%)
- dedupe: phosphosites duplicate site_ids 12 -> 0; kinase_substrate duplicate enzyme-site-direction 21 -> 0
- window15 centre mismatches: 2 (UniProt annotation-vs-sequence quirks, e.g. O75478_S6, Q15154_S159); tissue column: empty -> populated
- 5 isoform_only examples (window not uniquely found in canonical seq -> divergent isoform region):
    - O00159-3_S6s (base O00159, S6)
    - O00159-3_Y3y (base O00159, Y3)
    - O00429-2_S561s (base O00429, S561)
    - O00499-10_S270s (base O00499, S270)
    - O00499-10_S315s (base O00499, S315)

## Phosphosites (47,409)  [UNION of UniProt + OmniPath]
- by source: UniProt-only=33,085 | OmniPath-only=7,344 | both=6,980
- by residue: {'S': 36520, 'T': 7551, 'Y': 3338}
- by evidence: {'experimental': 29968, 'predicted/similarity': 10095, 'curated (OmniPath)': 7344, 'other': 2}
- with a known kinase: 15,524
- top kinases by substrate sites: {'CDK1': 1262, 'GSK3B': 1221, 'CDK2': 1216, 'PRKACA': 1079, 'PRKCA': 1043, 'CSNK2A1': 985, 'MAPK1': 966, 'SRC': 947, 'MAPK3': 748, 'MAPK14': 742, 'AKT1': 578, 'MAPK8': 559, 'ATM': 525, 'PRKCB': 509, 'EGF': 501, 'CDK5': 459, 'RPS6KA3': 425, 'CSNK2A2': 401, 'PLK1': 390, 'PRKCD': 364}

## Glycosites (16,904)
- by type: {'N-linked': 14937, 'O-linked': 1688, 'O-GlcNAc': 179, 'C-linked': 79, 'other': 21}
- by source: {'UniProt': 16904}

## MoTrPAC bridge
- feature rows parsed: 21,873; failed: 0
- mapping_status: {'canonical_ok': 18101, 'isoform_mapped': 5047, 'isoform_only': 682, 'residue_mismatch': 292}
- % MoTrPAC sites whose protein is a network node: 96.8%
- % found in phosphosites.csv: 43.1%
- % with a known kinase: 11.4%
- by tissue (n | %in_db | %known_kinase):
    - adipose: 4,309 | 51.1% | 10.5%
    - both: 5,189 | 71.2% | 22.1%
    - muscle: 14,624 | 30.7% | 7.8%

## Runtime / files
- runtime: 7s | UniProt release: UniProt Knowledgebase Release 2026_03
- proteins_ptm.csv: 1,626,004 bytes
- phosphosites.csv: 4,629,606 bytes
- glycosites.csv: 1,796,061 bytes
- kinase_substrate.csv: 4,882,918 bytes
- motrpac_feature_site_map.csv: 2,198,839 bytes
