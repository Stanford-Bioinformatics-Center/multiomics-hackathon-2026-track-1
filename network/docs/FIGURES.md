# Figures and captions

Titles inside figures are descriptive; the claims and method detail are here. Figure files are written to
`$HACK_FIG` (default `~/Desktop/output/hackathon`); Fig. 1 and Fig. M are also committed in `network/docs/figures/`.

**Fig. 1 | Endurance and resistance move disease-linked proteins back in different tissues** (`21_story_figure.R`).
(a) The 68 muscle proteins altered in T2D (Öhman 2021, p < 0.05; Amar et al. 2024 sets) plotted by their T2D z against
their mean normalised muscle response to endurance or resistance exercise (MoTrPAC, 0.5-24 h). Reversal = −Spearman
correlation; permutation p (10,000 shuffles of the T2D labels). Endurance moves proteins that are lower in T2D up and
higher in T2D down (reversal 0.31, p 0.012); resistance does not (−0.07, p 0.56); endurance − resistance p 0.007.
(b) The T2D proteins plus the network nodes linked to two or more of them (42 nodes; more connected than degree-matched
random seed sets, p 0.039). Edges coloured red / blue when strong (top 25% of |w| over both arms) after endurance /
resistance only; node outline purple / orange = lower / higher in T2D, black = connector. (c, d) UK Biobank plasma
proteins associated with age (Sun 2023; 293 proteins, Bonferroni) or with future T2D (Gadd 2024; 297 proteins) against the
mean normalised blood-protein response (same Olink platform): resistance reverses both (0.47 and 0.38, p < 0.001),
endurance weakly or not (0.15, p 0.013; 0.01, p 0.90). (e) Every tissue-matched disease / ageing test (step 20), grouped
by readout tissue; filled = p < 0.05; * pooled cohorts, post hoc.

**Fig. M | How an edge is made** (`docs/make_method_figure.R`). (a) Hard layer: a physical database (STRING ≥ 700 for
protein pairs; Rhea for enzyme-metabolite pairs) decides whether an edge exists; A and E respond alike (dot product 0.57)
but are not linked, so they get no edge. (b) Soft layer: normalised exercise responses per arm; the metabolite is doubled
into the RNA and protein slots. (c) The weight of A-B is the dot product of the two vectors, per arm (0.50 endurance,
0.38 resistance). (d) One weighted network per arm on identical edges, and the arm-specific edges (strong ≥ 0.3 in one arm
only). All numbers are computed by the exnet engine from its tested toy example.

**Figure 17 (interactive pages)** — joint, gene and metabolite networks; filters for omes / tissues / times / arm;
modules with pathway names and CAMERA-PR tests; MoTrPAC phosphosite and glycosylation tags; T2D layers; arm-specific
edges (red endurance only / blue resistance only; joint network at 25%: 39 vs 96 edges).

**Figures 19a-c** — the three T2D stories on the 68 T2D-altered proteins: (a) their own connected pieces (not beyond
chance), (b) T2D proteins + connectors (p 0.039), (c) protein level with responding muscle phosphosites.
Minimal PTM tags: one tag per category with its count.

**Figures 20a-e** — (a) every disease set and test; (b, c) the chosen story of new T2D (UK Biobank incident T2D, blood)
and ageing (UK Biobank age, blood); (d, e) the supporting muscle sets (Kjærgaard pooled, post hoc; Ubaida-Mohien).
