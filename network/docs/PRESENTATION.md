# Presentation outline (Team 2-PAC) — mapped to the judging criteria

Judging: **Presentation** (scientific novelty, impact, content) · **Presentation / GitHub** (collaboration, methods &
approach) · **GitHub** (documentation & integrity, technical complexity).

| # | Slide | Show | Say (one line) | Criterion |
|---|---|---|---|---|
| 1 | Title — Team 2-PAC | Fig. 1 thumbnail | "Three stories: endurance fixes the T2D muscle proteome; resistance fixes the plasma proteins of ageing and future T2D." | content |
| 2 | The question | Track brief + T2D | "Exercise as medicine: which kind, for which tissue? We picked T2D, as the brief asks." | content, impact |
| 3 | The data | MoTrPAC 3 tissues x 3 omes; 7 disease / ageing proteomes incl. UK Biobank | "Tissue-matched disease data only — no liver-vs-muscle shortcuts." | methods |
| 4 | **The method** | Fig. M | "Physical databases decide WHETHER two molecules connect; exercise decides HOW STRONGLY — the same edges in both arms, so every difference is biology, not wiring." | novelty, methods |
| 5 | The network | Figure 17 live (arm-specific edges on) | "353 molecules, 704 edges; red = endurance-only, blue = resistance-only." | technical complexity |
| 6 | **Story 1 · Muscle T2D** | Fig. 1a-b | "Endurance reverses the T2D muscle proteome (p 0.007 vs resistance); the T2D proteins form a connected subgraph (p 0.039)." | impact |
| 7 | **Story 2 · Blood ageing** | Fig. 1c | "In plasma it flips: resistance reverses the 293 proteins of ageing (UK Biobank, p < 0.001)." | novelty, impact |
| 8 | **Story 3 · Blood future T2D** | Fig. 1d, Fig. 1e | "And the 297 plasma proteins that predict T2D years ahead (UK Biobank, p < 0.001) — resistance again. Supporting: the pooled Kjærgaard muscle cohorts (post hoc) show the same muscle gap." | impact, integrity |
| 9 | Honest limits | README roadblocks | "Direction matches, one acute bout, healthy adults — not treatment. What didn't work, and why we dropped it." | integrity |
| 10 | Reproducible, not a black box | run_all + tests + manifest | "One command, 10 minutes: 25 engine tests, 44 checks, 139 outputs byte-identical." | documentation, complexity |
| 11 | Team + next steps | Contributors table, roadmap | "Who did what; blood-cell-composition control and training data next." | collaboration |

Live demo (slide 5): open `17a_joint_network.html` → tick "arm-specific edges only" → click SRC (the hub) → switch the
node colour to "T2D change".
