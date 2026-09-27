# mnet Phase 4 - Rhea protein-metabolite layer

- edges (configured policy `cofactor_rule`): 17,130
  - catalysis: 15,076 | transport: 2,054
- edge counts by measured_currency_policy: {'cofactor_rule': 17130, 'drop': 16915, 'keep_all': 19493}

## Currency protein-edge counts: before -> after (cofactor rule)

| metabolite | before | after |
|---|---|---|
| ATP | 871 | 79 |
| ADP | 782 | 39 |
| AMP | 120 | 74 |
| GTP | 244 | 34 |
| NAD+ | 107 | 51 |
| FAD | 0 | 4 |

### Example kept reactions
- ATP: ['RHEA:11332', 'RHEA:11600', 'RHEA:12973', 'RHEA:14357', 'RHEA:14433']
- ADP: ['RHEA:11460', 'RHEA:11600', 'RHEA:12973', 'RHEA:13749', 'RHEA:13893']
- AMP: ['RHEA:10040', 'RHEA:10412', 'RHEA:11460', 'RHEA:11800', 'RHEA:12973']
- GTP: ['RHEA:13549', 'RHEA:13665', 'RHEA:15229', 'RHEA:17473', 'RHEA:22484']
- NAD+: ['RHEA:11800', 'RHEA:16301', 'RHEA:18629', 'RHEA:19149', 'RHEA:21360']
- FAD: ['RHEA:13729', 'RHEA:17237', 'RHEA:67492', 'RHEA:73147']
- protein nodes: 4,139 (new/in_string=False: 119)
- crosswalk methods: {'string_direct': 4016, 'via_ensp': 81, 'new_node': 38, 'via_gene': 5}
- metabolite/class nodes with an edge: 4,171 (measured: 2,990 edges)
- shared_generic edges: 818
