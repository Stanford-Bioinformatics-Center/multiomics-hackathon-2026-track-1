#!/usr/bin/env Rscript
# Figure 7 — Secretome
#
# Panels:  FIG7A  exerkine candidates: tissue direction against plasma
#                 direction                            (also writes ST7a)
#          FIG7B  CX3CL1 (fractalkine) across muscle and adipose transcript
#                 and Olink
#          FIG7C  CCN1 across muscle and adipose transcript, blood Olink and
#                 muscle protein
#          FIG7D  CCN1 against its receptor transcripts, within tissue, arm
#                 and timepoint
# Tables:  ST7a   exerkine candidate list               (written by FIG7A)
#
# There is no tables/ST7.R: ST7a is FIG7A's own numbers, written from the frame
# the panel was about to select its rows from, and config/table_map.json records
# it against this script.
#
# The COMPARTMENTS extracellular scores FIG7A reads are vendored under
# figure_7/sources/; see docs/external_dependencies.md.
#
# Needs consortium data access. FIG7D reads the adipose and muscle transcript
# qc matrices through load_qc().
#
#   Rscript figures/landscape/FIG7.R          every panel
#   Rscript figures/landscape/FIG7.R FIG7A    one panel

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(ggh4x)
  library(ggplot2)
  library(patchwork)
  library(WGCNA)
  # load_differential_analysis(), which plot_single_feature() calls, resolves its
  # lazy-loaded DA objects by name through the search path, so the Analysis
  # package has to be attached and not only namespace-qualified.
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG7_ED9_helpers.R"))

# ---- shared ----------------------------------------------------------------

# There is no shared loader here. The one load more than one panel makes — the
# four-ome secretome DA (FIG7A, FIG7D) — is memoised in helpers/FIG7_ED9.R
# instead, because ED9A and ED9B read it with the same arguments. FIG7D's
# load_qc() is this figure's alone and is read once.

# ---- FIG7A — exerkine candidates against the plasma direction --------------

# Genes whose plasma Olink protein rises after exercise and whose tissue
# feature (muscle or adipose transcript, protein or phosphosite; blood
# transcript) is also differentially abundant, restricted to genes with a
# COMPARTMENTS extracellular score of at least 4. One column per tissue-ome,
# each cell the tissue direction against the plasma direction. Row annotation:
# the exercise arm and the extracellular score.
#
# This panel also writes ST7a, the full candidate list at every score and
# plasma direction, from the same frame the heatmap selects its rows from.
# ED9B draws the rows this panel leaves out.

fig7a <- function() {
  panel_init("FIG7A")

  candidates <- exerkine_candidate_frame()
  write_st7a(candidates, plasma_da_proteins(secretome_da()))

  main <- candidates[candidates$score >= COMPARTMENTS_MAIN_SCORE &
                       candidates$blood_direction == "Up", , drop = FALSE]
  m <- candidate_heatmap_matrix(main)

  # Plasma is Up on every row, so only the two Up cells, Mixed and Absent occur.
  colors <- CANDIDATE_CELL_COLORS[c("Tissue Up & Blood Up", "Tissue Down & Blood Up",
                                    "Mixed", "Absent")]

  heatmap <- candidate_heatmap(
    m, colors,
    candidate_row_annotation(
      main, rownames(m), score_range = c(3.5, 5),
      baseline = 3.5, border = FALSE,
      axis_param = list(at = c(3.5, 4, 4.5, 5), labels = c("3.5", "4", "4.5", "5"))
    )
  )

  export_panel(heatmap, "FIG7A")
}

# ---- ST7a — exerkine candidate list ----------------------------------------

# Every candidate tissue x assay x contrast, before the score and plasma
# direction cut FIG7A draws. It is FIG7A's own input rather than a table of its
# own, which is why it is written here.
write_st7a <- function(candidates, plasma) {
  export_table(
    exerkine_candidate_table(candidates, plasma)[, table_spec("ST7a")$columns, drop = FALSE],
    "ST7a"
  )
}

# ---- FIG7B — CX3CL1 trajectories -------------------------------------------

# One secreted protein followed through three measurements: the transcript in
# muscle and in adipose, and the circulating protein by Olink. Exercise-group
# means with 95% confidence intervals across the acute timepoints, three
# sub-plots in a row under one collected legend.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.

fig7b <- function() {
  panel_init("FIG7B")
  export_panel(single_feature_plot("FIG7B"), "FIG7B")
}

# ---- FIG7C — CCN1 trajectories ---------------------------------------------

# The same protein through four measurements: the transcript in muscle and in
# adipose, the circulating protein by Olink, and the muscle protein by mass
# spectrometry. Two sub-plots per row under one collected legend.
#
# CCN1 is also ED3C, drawn the other way: sex-split, muscle and adipose only.
#
# One more entry of the single-feature catalog in single_feature_plots.R.

fig7c <- function() {
  panel_init("FIG7C")
  export_panel(single_feature_plot("FIG7C"), "FIG7C")
}

# ---- FIG7D — CCN1 against its receptor transcripts -------------------------

# Biweight midcorrelation of the CCN1 transcript with fifteen integrin, VEGF
# receptor and co-receptor transcripts, computed over the participants of one
# tissue x exercise arm x post-exercise timepoint. Stars are nominal p; a cell
# is outlined when the receptor is itself differentially abundant against
# control in that arm and timepoint.

DA_LABEL <- "DA gene\n(adj p < 0.05)"

fig7d <- function() {
  panel_init("FIG7D")

  cor_df <- ccn1_receptor_correlations(
    MotrpacHumanPreSuspensionData::load_qc(), secretome_da())

  strip_colors <- strip_themed(
    background_y = elem_list_rect(fill = unname(TISSUE_COLORS[levels(factor(cor_df$Tissue))])),
    background_x = elem_list_rect(fill = unname(EXERCISE_COLORS[levels(cor_df$group)]))
  )

  p <- ggplot(cor_df, aes(x = Timepoint, y = gene_symbol, fill = bicor)) +
    geom_tile(color = "#fcfcfb", linewidth = 0.5) +
    geom_tile(data = cor_df[cor_df$is_DA, , drop = FALSE],
              aes(color = DA_LABEL), fill = NA, linewidth = 0.4) +
    geom_text(aes(label = sig_stars), color = "black", size = 4, vjust = 0.75) +
    scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#e34948",
                         midpoint = 0, limits = c(-1, 1), name = "bicor") +
    scale_color_manual(name = NULL, values = stats::setNames("black", DA_LABEL)) +
    guides(color = guide_legend(override.aes = list(fill = "grey95"))) +
    scale_x_discrete(limits = CCN1_TIMEPOINTS,
                     labels = SECRETOME_TIMEPOINT_LABELS[CCN1_TIMEPOINTS]) +
    facet_grid2(Tissue ~ group, strip = strip_colors,
                scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )

  export_panel(p, "FIG7D")
}

run_panels(list(
  FIG7A = fig7a,
  FIG7B = fig7b,
  FIG7C = fig7c,
  FIG7D = fig7d
))
