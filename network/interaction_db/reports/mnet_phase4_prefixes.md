# mnet Phase 4 pre-fixes

## 1. Non-fatty-acids removed from lipid scope
- Removed: **Aminobutyric acid** (amino acid, not a fatty acid); kept in curation/manual_metabolite_mappings.csv. The 4 hydroxy fatty acids (Hydroxyoctanoic/decanoic/dodecanoic/tetradecanoic) are genuine FAs and kept.

## 2. FA chain-length classes: before -> after (distinct human reactions)

| class | chebi_ids | before | after |
|---|---|---|---|
| LIPIDCLASS:FA-short-chain | nan | 60 | 0 |
| LIPIDCLASS:FA-medium-chain | CHEBI:28868;CHEBI:59558 | 60 | 63 |
| LIPIDCLASS:FA-long-chain | CHEBI:28868;CHEBI:57560 | 8 | 68 |
| LIPIDCLASS:FA-very-long-chain | CHEBI:28868;CHEBI:58950 | 60 | 64 |
- specific chain-length anions now used (medium CHEBI:59558, long CHEBI:57560, very-long CHEBI:58950) PLUS the shared general fatty-acid anion CHEBI:28868; short-chain has no specific anion in human Rhea.

## 3. Ether classes

- LIPIDCLASS:PC-ether: 56 human reactions, chebis=CHEBI:30909;CHEBI:36702;CHEBI:36707;CHEBI:57643;CHEBI:77286
- LIPIDCLASS:PE-ether: 32 human reactions, chebis=CHEBI:60520;CHEBI:64612;CHEBI:76168;CHEBI:77290
- LIPIDCLASS:LPC-ether: 63 human reactions, chebis=CHEBI:30909;CHEBI:36702;CHEBI:36707;CHEBI:57875;CHEBI:58168;CHEBI:77287
- LIPIDCLASS:LPE-ether: 34 human reactions, chebis=CHEBI:64381;CHEBI:65213;CHEBI:76168;CHEBI:77288
(alkyl/alkenyl generics used where present; see how_found/confidence in curation/lipid_class_to_chebi.csv)

## 4. Shared generic ChEBI (claimed by >1 class node)

- shared generic ChEBIs: **21**
  - CHEBI:28868 (fatty acid anion) claimed by 3 class nodes
  - CHEBI:30909 (1-alkyl-sn-glycero-3-phosphocholine) claimed by 3 class nodes
  - CHEBI:77286 (1-(Z)-alk-1-enyl-2-acyl-sn-glycero-3-phosphocholine) claimed by 3 class nodes
  - CHEBI:36702 (2-acyl-1-alkyl-sn-glycero-3-phosphocholine) claimed by 3 class nodes
  - CHEBI:60520 (1-alkyl-2-acyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 3 class nodes
  - CHEBI:76168 (1-alkyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 3 class nodes
  - CHEBI:77290 (1-(Z)-alk-1-enyl-2-acyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 3 class nodes
  - CHEBI:17950 (beta-D-galactosyl-(1->4)-beta-D-glucosyl-(1<->1)-N-acylsphingosine) claimed by 2 class nodes
  - CHEBI:22801 (beta-D-glucosyl-N-acylsphingosine) claimed by 2 class nodes
  - CHEBI:31488 (N-acylsphinganine) claimed by 2 class nodes
  - CHEBI:52639 (N-acylsphingosine) claimed by 2 class nodes
  - CHEBI:57674 (N-acylsphingosine 1-phosphate(2-)) claimed by 2 class nodes
  - CHEBI:57875 (2-acyl-sn-glycero-3-phosphocholine) claimed by 2 class nodes
  - CHEBI:58168 (1-O-acyl-sn-glycero-3-phosphocholine) claimed by 2 class nodes
  - CHEBI:77287 (1-(Z)-alk-1-enyl-sn-glycero-3-phosphocholine) claimed by 2 class nodes
  - CHEBI:36707 (2-acetyl-1-alkyl-sn-glycero-3-phosphocholine) claimed by 2 class nodes
  - CHEBI:64381 (1-acyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 2 class nodes
  - CHEBI:65213 (2-acyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 2 class nodes
  - CHEBI:77288 (1-(Z)-alk-1-enyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 2 class nodes
  - CHEBI:57643 (1,2-diacyl-sn-glycero-3-phosphocholine) claimed by 2 class nodes
  - CHEBI:64612 (1,2-diacyl-sn-glycero-3-phosphoethanolamine zwitterion) claimed by 2 class nodes
- total (class, shared-generic) claims across classes: **49** (each generates Rhea catalysis edges to EACH claiming class node in Phase 4, marked shared_generic=True)

