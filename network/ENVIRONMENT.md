# Environment and input fingerprints

Recorded by `network/00_environment.R` on 2026-09-26. Re-run it (or `network/run_all.sh`) and compare: the same
versions and input checksums should give identical outputs (checked with `network/99_manifest.R`).

## Software

| Tool | Version |
|---|---|
| R | R version 4.4.3 (2025-02-28) |
| Python | Python 3.12.4 |
| pandoc | 3.10 |
| TinyTeX (pdflatex) | pdfTeX 3.141592653-2.6-1.40.29 (TeX Live 2026) |
| platform | aarch64-apple-darwin20 |

## R packages used by the pipeline

| Package | Version |
|---|---|
| data.table | 1.18.2.1 |
| igraph | 2.2.2 |
| ggplot2 | 3.5.2 |
| ggrepel | 0.9.8 |
| ggforce | 0.5.0 |
| ggnewscale | 0.5.2 |
| scales | 1.4.0 |
| nanoparquet | 0.4.3 |
| jsonlite | 2.0.0 |
| visNetwork | 2.1.4 |
| htmlwidgets | 1.6.4 |
| htmltools | 0.5.9 |
| rmarkdown | 2.30 |
| tinytex | 0.58 |
| limma | 3.62.2 |
| readxl | 1.4.5 |
| checkmate | 2.3.4 |
| testthat | 3.3.2 |
| pkgload | 1.5.0 |
| patchwork | 1.3.2 |
| exnet | 1.0.0 |
| TMSig | 1.0.0 |
| fgsea | 1.32.4 |
| MotrpacHumanPreSuspensionAnalysis | 0.2.4 |
| MotrpacHumanPreSuspensionData | 0.0.1.4 |
| MotrpacHumanPreSuspension | 0.0.1.2 |
| Matrix | 1.7.5 |

## External inputs (MD5 fingerprints)

| Source | File | Bytes | Modified | MD5 |
|---|---|---|---|---|
| team mnet resource (MNET_DIR) | edges.parquet | 1,414,570 | 2026-09-26 17:26 | de2639bb258dfe50348952c8c33f684b |
| team mnet resource (MNET_DIR) | glycosites.csv | 1,796,061 | 2026-09-26 17:26 | 293bef5d2ff4fd7b07a70b2945debaaf |
| team mnet resource (MNET_DIR) | kinase_substrate.csv | 4,882,918 | 2026-09-26 17:26 | 526d9a5d3021ec12f9ab16312f2fcac4 |
| team mnet resource (MNET_DIR) | mnet_final_coverage.md | 6,617 | 2026-09-26 17:26 | 4c104ed2ac849b594e6e2af336a887b6 |
| team mnet resource (MNET_DIR) | mnet_phase6_ptm.md | 2,595 | 2026-09-26 17:26 | 341b2519519fcd5564928853efd6fee5 |
| team mnet resource (MNET_DIR) | motrpac_feature_site_map.csv | 2,198,839 | 2026-09-26 17:26 | 714ab0c8ebc932ed5f24ee69643f8204 |
| team mnet resource (MNET_DIR) | nodes.csv | 1,667,845 | 2026-09-26 17:26 | 9dfb61b9914c8642505d3bd407adfbab |
| team mnet resource (MNET_DIR) | phosphosites.csv | 4,629,606 | 2026-09-26 17:26 | c0180bff376f5c2249e80a4009ccd534 |
| team mnet resource (MNET_DIR) | proteins_ptm.csv | 1,626,004 | 2026-09-26 17:26 | cf772af47e699cd47ea5072b5f54839f |
| team mnet resource (MNET_DIR) | README.md | 3,233 | 2026-09-26 17:26 | ab6a128f64b73a2d61f7ed98d819c72a |
| legacy curated STRING file (STRING_PARQUET; used only with EDGE_SOURCE=legacy) | Metabolomics_database_watershed_template_data_p_value_string_network_ge700.parquet | 585,329 | 2026-09-26 10:30 | 3a1fc7a9776d76553b8ac01cb152fdb3 |
| Amar et al. 2024 disease proteomics sets, processed in Venus week 6 (DISEASE_SCORES; step 18d) | disease_scores.csv.gz | 995,726 | 2026-09-24 10:22 | d4a098ad84762f24f7f4335579b69420 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | MOESM3.xlsx | 7,733,778 | 2026-09-26 22:27 | e993d67cf5833b4c96fcc3780fb137d2 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc1.xlsx | 10,851,954 | 2026-09-26 21:28 | f3a74f2cfbe065085921edffd5d39dc9 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc2.xlsx | 12,975,402 | 2026-09-26 21:28 | f5c8604836d1b2577b7d7dbbbaa40781 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc3.xlsx | 11,073 | 2026-09-26 21:28 | 3f4396042f86a76b0e6cc9ae964a52c7 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | COMMIT_SHA | 41 | 2026-09-26 21:26 | d821d8607d06e21f8630e12f88828c0b |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | Exprs_adipose_clean.txt | 3,004,748 | 2026-09-26 21:26 | 434e4eef00ce02aa0821760588f608f9 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | File_names.txt | 10,322 | 2026-09-26 21:26 | aad8dc9f38a10770762b69a60dcc09e1 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | Limma_and_GeneSet_analysis.R | 68,635 | 2026-09-26 21:26 | 1728905219bcf34a67e1ac26f61c1ab3 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | Preprocessing_and_Figure1.R | 17,592 | 2026-09-26 21:26 | 923e73d5cc787dcfe6598dc316622131 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | README.md | 1,112 | 2026-09-26 21:26 | 817af329f7cf5a5a40f9d374fd7e61d3 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc2.xlsx | 11,532 | 2026-09-26 21:28 | 63b0aaf2183f8c634f2382b5e4541a1f |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc3.xlsx | 2,191,560 | 2026-09-26 21:28 | 961f0d66d97b64a9ccb11b4cf0b29228 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc4.xlsx | 13,509 | 2026-09-26 21:28 | 542ef3d259746943afa6771dcb77865f |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc5.xlsx | 29,658,815 | 2026-09-26 21:28 | 89a190b3665c2dd50e2d487e34ff09ae |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | mmc6.xlsx | 2,099,402 | 2026-09-26 21:28 | 77a87cbabeeda62adf19f127c3f99b33 |
| published disease tables for step 20 (DISEASE_EXT; Kjærgaard 2025, Needham 2024, Larsen 2023, Sun 2023, Gadd 2024; step 20f) | MOESM3.xlsx | 34,411,394 | 2026-09-26 22:28 | aeb2acdc5586034e91197658a10a7d70 |
| Ubaida-Mohien 2019 muscle ageing proteome, Amar et al. repository copy (UBAIDA; step 20) | ubaida_mohien_ 2019_elife_stat.csv | 3,075,476 | 2026-08-05 20:59 | 2b48b6cb930d00dcaed765d3f8a9b74d |
| Rhea cache (HACK_EXT), release 142 / 2026-09-02 | chebi_pH7_3_mapping.tsv | 3,099,222 | 2026-09-26 13:02 | adf14687a450d8366f9b364101202635 |
| Rhea cache (HACK_EXT), release 142 / 2026-09-02 | rhea-kegg.reaction.gz | 1,556,033 | 2026-09-26 13:02 | 325907b498cfb6c89aa0fa6ff8605bab |
| Rhea cache (HACK_EXT), release 142 / 2026-09-02 | rhea-release.properties | 53 | 2026-09-26 13:02 | 210aad03260f7186ef97170962e1666e |
| Rhea cache (HACK_EXT), release 142 / 2026-09-02 | rhea2uniprot_sprot.tsv | 8,835,509 | 2026-09-26 13:02 | 828d8e03c8f1549eb759591e3010ba14 |
| Rhea cache (HACK_EXT), release 142 / 2026-09-02 | uniprot_human_swissprot_accessions.txt | 144,669 | 2026-09-26 13:12 | f664ecdd637b545af27fa0eb300721c7 |
| GlyGen API cache (inventory/glygen_cache), service release now: 2.11.1 | 471 per-protein files | 4,359,891 | 2026-09-26 16:49 | 22782fe766108e325d9bb3b0a4119e9f |
| MoTrPAC packages (installed) | MotrpacHumanPreSuspensionAnalysis | NA | NA | version 0.2.4 |
| MoTrPAC packages (installed) | MotrpacHumanPreSuspensionData | NA | NA | version 0.0.1.4 |

Full `sessionInfo()` is written to `$HACK_OUT/00_environment/session_info.txt` (not committed).
