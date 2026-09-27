# FIG7_ED9.R — the secretome computation Figure 7 and Extended Data 9 share.
#
# One file because both figures read it:
#   secretome_da()                FIG7A, FIG7D, ED9A, ED9B
#   exerkine_candidate_frame()    FIG7A, ED9A, ED9B
#   candidate_heatmap_matrix()    FIG7A, ED9B
#   candidate_heatmap()           the same two
#   exerkine_candidate_table()    FIG7A, for ST7a
#   ccn1_receptor_correlations()  FIG7D
#   ccn1_cross_tissue()           ED9D
#
# Ported from precovid-analyses PR #106 (b410508),
# figures/landscape/figure_8/lanscape_secretome_github_091426.R.
#
# Sourced after lib/panel_export.R and lib/highlights.R, which define
# panel_source() and highlight_feature(). Attaches no packages: everything is
# namespace-qualified, so the script that sources this file decides what it
# attaches.

SECRETOME_FDR <- 0.05

# Tissue-side omes. Blood contributes its transcriptome only; blood Olink is the
# plasma side. Order here is the row order of the candidate frame, and so of
# the heatmaps.
SECRETOME_TISSUE_ASSAYS <- list(
  muscle  = c("transcript-rna-seq", "prot-pr", "prot-ph"),
  adipose = c("transcript-rna-seq", "prot-pr", "prot-ph"),
  blood   = "transcript-rna-seq"
)

# COMPARTMENTS integrated channel, best score over the four extracellular
# locations. FIG7A shows genes at or above COMPARTMENTS_MAIN_SCORE; ED9B the rest.
COMPARTMENTS_LOCATIONS <- c("Extracellular region", "Extracellular space",
                            "Extracellular exosome", "Extracellular vesicle")
COMPARTMENTS_MAIN_SCORE <- 4

CCN1_TRANSCRIPT_ID <- "ENSG00000142871.18"
CCN1_RECEPTOR_GENES <- c("KDR", "ITGAV", "ITGA5", "ITGB3", "ITGB1", "ITGA2B", "ITGA6",
                         "ITGAD", "ITGAM", "ITGB2", "ITGB5", "LRP1", "SDC4", "TLR2",
                         "TLR4")
CCN1_TIMEPOINTS <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
CCN1_GROUPS <- c(ADUEndur = "EE", ADUResist = "RE", ADUControl = "CON")
CCN1_MIN_SAMPLES <- 5

SECRETOME_TIMEPOINT_LABELS <- c(
  pre_exercise      = "Pre",
  during_20_min     = "D20M",
  during_40_min     = "D40M",
  post_10_min       = "P10M",
  post_15_30_45_min = "P15-45M",
  post_3.5_4_hr     = "P3.5-4H",
  post_24_hr        = "P24H"
)

ST7A_TIMEPOINT_LABELS <- c(
  during_20_min     = "During_20min",
  during_40_min     = "During_40min",
  post_10_min       = "Post_10min",
  post_15_30_45_min = "Post_15_30_45min",
  post_3.5_4_hr     = "Post_3.5_4hr",
  post_24_hr        = "Post_24hr"
)

# Temporal concordance, per gene: a tissue event in the early window with a
# plasma event at or after 15-45 min, or a tissue event at 3.5-4 h with a plasma
# event at 3.5-4 h or 24 h. The early window is every timepoint under an hour,
# during and post, as b410508 groups them.
TEMPORAL_EARLY_TISSUE <- c("during_20_min", "during_40_min", "post_10_min",
                           "post_15_30_45_min")
TEMPORAL_EARLY_PLASMA <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
TEMPORAL_MID_TISSUE   <- "post_3.5_4_hr"
TEMPORAL_MID_PLASMA   <- c("post_3.5_4_hr", "post_24_hr")

EXERCISE_COLORS <- c(
  EE   = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS["ADUEndur"]),
  RE   = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS["ADUResist"]),
  CON  = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS["ADUControl"]),
  Both = "black"
)

TISSUE_DISPLAY <- c(adipose = "Adipose", muscle = "Muscle", blood = "Blood")

TISSUE_COLORS <- c(
  Muscle  = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS["muscle"]),
  Adipose = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS["adipose"]),
  Blood   = unname(MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS["blood"])
)

# ---- differential analysis --------------------------------------------------

#' The four omes the secretome reads, exercise-vs-control contrasts only, one
#' frame in the legacy's block order (SECRETOME_TISSUE_ASSAYS, then blood Olink),
#' each block in its table's own row order. Gene symbols come from
#' HUMAN_FEATURE_TO_GENE by assay and feature_id, phosphosites included.
#'
#' Lazy and memoised here rather than in either figure script: FIG7A, FIG7D,
#' ED9A and ED9B all want the same frame with no arguments, and one run of
#' either script asks for it twice.
secretome_da <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
      stop("secretome_da() needs MotrpacHumanPreSuspensionAnalysis on the search ",
           "path: declare it in the panel's data_packages in config/panel_map.json.",
           call. = FALSE)
    }
    tables <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = c("transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol"),
      selected_tissues = "all",
      single_matrix = FALSE,
      verbose = FALSE
    )
    blocks <- c(
      unlist(lapply(names(SECRETOME_TISSUE_ASSAYS), function(tissue) {
        lapply(SECRETOME_TISSUE_ASSAYS[[tissue]], function(assay) tables[[tissue]][[assay]])
      }), recursive = FALSE),
      list(tables[["blood"]][["prot-ol"]])
    )
    da <- do.call(rbind, lapply(blocks, as.data.frame))
    da <- da[da$contrast_type == "exercise_with_controls" &
               da$contrast_category %in% c("EE-CON", "RE-CON"), , drop = FALSE]
    da$tissue <- as.character(da$tissue)
    da$assay <- as.character(da$assay)
    da$feature_id <- as.character(da$feature_id)
    da$Timepoint <- as.character(da$Timepoint)

    mapping <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
    mapping <- mapping[mapping$assay %in% unique(da$assay), c("assay", "feature_id", "gene_symbol")]
    idx <- match(paste(da$assay, da$feature_id), paste(mapping$assay, mapping$feature_id))
    da$gene_symbol <- as.character(mapping$gene_symbol)[idx]

    da$Exercise <- sub("-CON$", "", as.character(da$contrast_category))
    da$Direction <- ifelse(da$logFC > 0, "Up", "Down")
    rownames(da) <- NULL
    cache <<- da
    cache
  }
})

#' Plasma Olink proteins differentially abundant after exercise.
plasma_da_proteins <- function(da) {
  p <- da[da$tissue == "blood" & da$assay == "prot-ol" &
            da$adj_p_value < SECRETOME_FDR & !is.na(da$gene_symbol), , drop = FALSE]
  rownames(p) <- NULL
  p
}

#' One row per plasma protein: Up or Down when every significant contrast
#' agrees, Mixed otherwise, and the logFC of largest magnitude.
plasma_direction <- function(plasma) {
  by_gene <- split(plasma$logFC, plasma$gene_symbol)
  data.frame(
    gene_symbol = names(by_gene),
    blood_direction = vapply(by_gene, function(x) {
      if (all(x > 0)) "Up" else if (all(x < 0)) "Down" else "Mixed"
    }, character(1)),
    blood_logFC = vapply(by_gene, function(x) x[which.max(abs(x))], numeric(1)),
    stringsAsFactors = FALSE
  )
}

#' Best extracellular COMPARTMENTS score per gene symbol, from the vendored
#' subset of the integrated channel.
#'
#' No guard. The file is vendored under figure_7/sources/ rather than produced
#' by a fit, so panel_source() resolves it from the repo and
#' ACUTE_COMPARTMENTS_TSV only points the reader at a newer copy.
compartments_extracellular <- function(path = NULL) {
  if (is.null(path)) {
    path <- panel_source("figure_7", "compartments_extracellular.tsv",
                         env_var = "ACUTE_COMPARTMENTS_TSV")
  }
  x <- read.csv(path, sep = "\t", check.names = FALSE, quote = "",
                stringsAsFactors = FALSE)
  x <- x[x$location %in% COMPARTMENTS_LOCATIONS, , drop = FALSE]
  tapply(x$score, x$gene_symbol, max)
}

#' Tissue features differentially abundant after exercise whose gene symbol is
#' a plasma DA protein: the exerkine candidate frame every secretome panel and
#' ST7a reads. One row per tissue x assay x contrast.
exerkine_candidates <- function(da, compartments) {
  plasma <- plasma_da_proteins(da)

  blocks <- lapply(names(SECRETOME_TISSUE_ASSAYS), function(tissue) {
    lapply(SECRETOME_TISSUE_ASSAYS[[tissue]], function(assay) {
      da[da$tissue == tissue & da$assay == assay, , drop = FALSE]
    })
  })
  tis <- do.call(rbind, unlist(blocks, recursive = FALSE))
  tis <- tis[tis$adj_p_value < SECRETOME_FDR &
               tis$gene_symbol %in% plasma$gene_symbol, , drop = FALSE]

  arms <- tapply(tis$Exercise, tis$gene_symbol, function(x) {
    if (length(unique(x)) > 1) "Both" else unique(x)
  })
  tis$Exercise2 <- as.character(arms[tis$gene_symbol])

  tis$score <- as.numeric(compartments[tis$gene_symbol])
  if (anyNA(tis$score)) {
    stop("no extracellular COMPARTMENTS score for: ",
         paste(sort(unique(tis$gene_symbol[is.na(tis$score)])), collapse = ", "),
         call. = FALSE)
  }

  direction <- plasma_direction(plasma)
  idx <- match(tis$gene_symbol, direction$gene_symbol)
  tis$blood_direction <- direction$blood_direction[idx]
  tis$blood_logFC <- direction$blood_logFC[idx]

  rownames(tis) <- NULL
  tis
}

#' The one candidate frame FIG7A, ED9A and ED9B draw from, and ST7a is written
#' off.
#'
#' Memoised, and memoised HERE rather than in a figure script: its callers sit
#' in both FIG7.R and ED9.R, it takes no arguments, so every caller wants the
#' same frame, and one secretome_da() pass serves a whole figure.
exerkine_candidate_frame <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    cache <<- exerkine_candidates(secretome_da(), compartments_extracellular())
    cache
  }
})

# ---- candidate heatmaps (FIG7A, ED9B) --------------------------------------

CANDIDATE_COLUMN_LABELS <- c(
  "muscle-transcript-rna-seq"  = "Muscle - Transcript",
  "muscle-prot-pr"             = "Muscle - Protein",
  "muscle-prot-ph"             = "Muscle - Phospho",
  "adipose-transcript-rna-seq" = "Adipose - Transcript",
  "adipose-prot-pr"            = "Adipose - Protein",
  "adipose-prot-ph"            = "Adipose - Phospho",
  "blood-transcript-rna-seq"   = "Blood - Transcript"
)

CANDIDATE_CELL_COLORS <- c(
  "Tissue Up & Blood Up"     = "#C84C05",
  "Tissue Down & Blood Up"   = "#FFB200",
  "Tissue Up & Blood Down"   = "#DDEB9D",
  "Tissue Down & Blood Down" = "#638C6D",
  "Mixed"                    = "purple",
  "Absent"                   = "white"
)

CANDIDATE_CELL_LEGEND <- c(
  "Tissue Up & Blood Up"     = "Blood Up + Tissue Up",
  "Tissue Down & Blood Up"   = "Blood Up + Tissue Down",
  "Tissue Up & Blood Down"   = "Blood Down + Tissue Up",
  "Tissue Down & Blood Down" = "Blood Down + Tissue Down",
  "Mixed"                    = "Mixed",
  "Absent"                   = "Absent"
)

candidate_cell_label <- function(tissue_logFC, blood_direction) {
  ifelse(tissue_logFC > 0 & blood_direction == "Up",   "Tissue Up & Blood Up",
  ifelse(tissue_logFC < 0 & blood_direction == "Up",   "Tissue Down & Blood Up",
  ifelse(tissue_logFC > 0 & blood_direction == "Down", "Tissue Up & Blood Down",
  ifelse(tissue_logFC < 0 & blood_direction == "Down", "Tissue Down & Blood Down",
         "Mixed"))))
}

#' Gene x (tissue - ome) matrix of direction labels. A gene with both
#' directions in one column is "Mixed"; a column it was never significant in
#' is "Absent". Rows follow first appearance in the frame, columns follow
#' CANDIDATE_COLUMN_LABELS.
candidate_heatmap_matrix <- function(candidates) {
  column <- paste(candidates$tissue, candidates$assay, sep = "-")
  unknown <- setdiff(unique(column), names(CANDIDATE_COLUMN_LABELS))
  if (length(unknown) > 0) {
    stop("tissue-assay with no heatmap column: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  label <- candidate_cell_label(candidates$logFC, candidates$blood_direction)

  genes <- unique(candidates$gene_symbol)
  columns <- names(CANDIDATE_COLUMN_LABELS)[names(CANDIDATE_COLUMN_LABELS) %in% column]
  m <- matrix("Absent", nrow = length(genes), ncol = length(columns),
              dimnames = list(genes, unname(CANDIDATE_COLUMN_LABELS[columns])))
  cells <- split(label, list(candidates$gene_symbol, column), drop = TRUE, sep = "\t")
  for (key in names(cells)) {
    parts <- strsplit(key, "\t", fixed = TRUE)[[1]]
    gene <- parts[1]
    col <- unname(CANDIDATE_COLUMN_LABELS[parts[2]])
    values <- unique(cells[[key]])
    m[gene, col] <- if (length(values) > 1) "Mixed" else values
  }
  m
}

#' Exercise arm and extracellular score beside each heatmap row, indexed by
#' the matrix's row names.
candidate_row_annotation <- function(candidates, genes, score_range, ...) {
  meta <- unique(candidates[, c("gene_symbol", "Exercise2", "score")])
  if (anyDuplicated(meta$gene_symbol)) {
    stop("a candidate gene carries more than one arm or score: ",
         paste(meta$gene_symbol[duplicated(meta$gene_symbol)], collapse = ", "),
         call. = FALSE)
  }
  meta <- meta[match(genes, meta$gene_symbol), , drop = FALSE]
  fill <- circlize::colorRamp2(score_range, c("white", "green"))
  ComplexHeatmap::rowAnnotation(
    Exercise = meta$Exercise2,
    Score = ComplexHeatmap::anno_barplot(
      meta$score,
      gp = grid::gpar(fill = fill(meta$score)),
      ylim = score_range,
      ...
    ),
    col = list(Exercise = EXERCISE_COLORS[c("EE", "RE", "Both")])
  )
}

candidate_heatmap <- function(m, colors, row_annotation) {
  ComplexHeatmap::Heatmap(
    m,
    name = "LogFC-Blood Direction",
    col = colors,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Blood & Tissue Direction",
      at = names(colors),
      labels = unname(CANDIDATE_CELL_LEGEND[names(colors)]),
      legend_gp = grid::gpar(fill = colors)
    ),
    right_annotation = row_annotation,
    border = TRUE
  )
}

# ---- ST7a ------------------------------------------------------------------

#' Genes whose tissue event precedes or coincides with a plasma event.
temporally_concordant_genes <- function(candidates, plasma) {
  tissue_tp <- split(candidates$Timepoint, candidates$gene_symbol)
  plasma_tp <- split(plasma$Timepoint, plasma$gene_symbol)
  genes <- names(tissue_tp)
  hit <- vapply(genes, function(g) {
    tp <- tissue_tp[[g]]
    pp <- plasma_tp[[g]]
    (any(tp %in% TEMPORAL_EARLY_TISSUE) && any(pp %in% TEMPORAL_EARLY_PLASMA)) ||
      (any(tp %in% TEMPORAL_MID_TISSUE) && any(pp %in% TEMPORAL_MID_PLASMA))
  }, logical(1))
  genes[hit]
}

#' ST7a, one row per candidate tissue x assay x contrast.
exerkine_candidate_table <- function(candidates, plasma) {
  label <- ST7A_TIMEPOINT_LABELS[candidates$Timepoint]
  if (anyNA(label)) {
    stop("timepoint with no ST7a label: ",
         paste(unique(candidates$Timepoint[is.na(label)]), collapse = ", "),
         call. = FALSE)
  }
  temporal <- temporally_concordant_genes(candidates, plasma)
  plasma_direction <- ifelse(candidates$blood_logFC > 0, "Up", "Down")
  data.frame(
    gene_symbol = candidates$gene_symbol,
    assay = candidates$assay,
    tissue = candidates$tissue,
    contrast = paste0(candidates$Exercise, "_", unname(label)),
    Modality = candidates$Exercise2,
    Extracellular_score = candidates$score,
    Tissue_direction = candidates$Direction,
    Plasma_direction = plasma_direction,
    Directionality_concordant = ifelse(candidates$Direction == plasma_direction, "Yes", "No"),
    Temporally_concordant = ifelse(candidates$gene_symbol %in% temporal, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# ---- CCN1 correlations (FIG7D, ED9D) --------------------------------------

#' A tissue's normalised transcript matrix and the metadata rows it has columns
#' for, vial labels as character.
tissue_transcripts <- function(qc, tissue) {
  ome <- qc[[tissue]][["transcript-rna-seq"]]
  norm <- ome$qc_norm
  meta <- ome$sample_metadata
  meta$vialLabel <- as.character(meta$vialLabel)
  meta <- meta[meta$vialLabel %in% colnames(norm), , drop = FALSE]
  if (!CCN1_TRANSCRIPT_ID %in% rownames(norm)) {
    stop("CCN1 (", CCN1_TRANSCRIPT_ID, ") is not in the ", tissue,
         " transcript qc matrix", call. = FALSE)
  }
  list(norm = norm, meta = meta)
}

#' Biweight midcorrelation of CCN1 with each receptor transcript, within one
#' tissue, one exercise arm and one timepoint, over that cell's participants.
#' `p` is nominal; `p_adj` is Benjamini-Hochberg over the tissue.
ccn1_receptor_correlations <- function(qc, da) {
  mapping <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  mapping <- mapping[mapping$assay == "transcript-rna-seq", c("feature_id", "gene_symbol")]
  mapping$feature_id <- as.character(mapping$feature_id)
  mapping$gene_symbol <- as.character(mapping$gene_symbol)
  target_ids <- mapping$feature_id[mapping$gene_symbol %in% CCN1_RECEPTOR_GENES]

  rows <- list()
  for (tissue in c("adipose", "muscle")) {
    tt <- tissue_transcripts(qc, tissue)
    for (grp in names(CCN1_GROUPS)[CCN1_GROUPS %in% c("EE", "RE")]) {
      for (tp in CCN1_TIMEPOINTS) {
        cell <- tt$meta[tt$meta$Timepoint == tp & tt$meta$randomGroupCode == grp, , drop = FALSE]
        if (nrow(cell) < CCN1_MIN_SAMPLES) next
        expr <- tt$norm[, cell$vialLabel, drop = FALSE]
        targets <- intersect(target_ids, rownames(expr))
        if (length(targets) == 0) next
        ccn1 <- as.numeric(expr[CCN1_TRANSCRIPT_ID, ])
        target_mat <- t(as.matrix(expr[targets, , drop = FALSE]))
        res <- WGCNA::bicorAndPvalue(target_mat, ccn1, use = "pairwise.complete.obs")
        rows[[length(rows) + 1]] <- data.frame(
          Tissue = unname(TISSUE_DISPLAY[tissue]),
          tissue = tissue,
          feature_id = targets,
          group = unname(CCN1_GROUPS[grp]),
          Timepoint = tp,
          n = nrow(cell),
          bicor = as.numeric(res$bicor),
          p = as.numeric(res$p),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$gene_symbol <- mapping$gene_symbol[match(out$feature_id, mapping$feature_id)]
  out$p_adj <- stats::ave(out$p, out$tissue, FUN = function(p) stats::p.adjust(p, "BH"))

  transcript_da <- da[da$assay == "transcript-rna-seq" & da$tissue %in% c("adipose", "muscle"),
                      c("tissue", "feature_id", "Exercise", "Timepoint", "adj_p_value"),
                      drop = FALSE]
  key_out <- paste(out$tissue, out$feature_id, out$group, out$Timepoint)
  key_da <- paste(transcript_da$tissue, transcript_da$feature_id,
                  transcript_da$Exercise, transcript_da$Timepoint)
  if (anyDuplicated(key_da)) {
    stop("more than one exercise-vs-control row per tissue x feature x arm x timepoint",
         call. = FALSE)
  }
  out$adj_p_value <- transcript_da$adj_p_value[match(key_out, key_da)]
  out$is_DA <- !is.na(out$adj_p_value) & out$adj_p_value < SECRETOME_FDR
  out$sig_stars <- ifelse(out$p < 0.001, "***",
                   ifelse(out$p < 0.01, "**",
                   ifelse(out$p < 0.05, "*", "")))
  out$Timepoint <- factor(out$Timepoint, levels = CCN1_TIMEPOINTS)
  out$group <- factor(out$group, levels = c("EE", "RE"))
  rownames(out) <- NULL
  out
}

#' CCN1 transcript in adipose and in muscle for every participant with both
#' samples at a timepoint, with the Pearson correlation per exercise group.
ccn1_cross_tissue <- function(qc) {
  mus <- tissue_transcripts(qc, "muscle")
  adi <- tissue_transcripts(qc, "adipose")
  points <- list()
  for (tp in unique(mus$meta$Timepoint)) {
    m <- mus$meta[mus$meta$Timepoint == tp, , drop = FALSE]
    a <- adi$meta[adi$meta$Timepoint == tp, , drop = FALSE]
    shared <- intersect(m$pid, a$pid)
    if (length(shared) == 0) next
    m <- m[match(shared, m$pid), , drop = FALSE]
    a <- a[match(shared, a$pid), , drop = FALSE]
    points[[tp]] <- data.frame(
      Timepoint = tp,
      pid = shared,
      group = unname(CCN1_GROUPS[as.character(m$randomGroupCode)]),
      Adipose = as.numeric(adi$norm[CCN1_TRANSCRIPT_ID, a$vialLabel]),
      Muscle = as.numeric(mus$norm[CCN1_TRANSCRIPT_ID, m$vialLabel]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, points)
  rownames(out) <- NULL
  out
}

#' Pearson r and p per timepoint x group, and the caption line per timepoint.
ccn1_cross_tissue_stats <- function(points) {
  cells <- split(points, list(points$Timepoint, points$group), drop = TRUE)
  stats <- do.call(rbind, lapply(cells, function(d) {
    test <- stats::cor.test(d$Adipose, d$Muscle)
    data.frame(Timepoint = d$Timepoint[1], group = d$group[1], n = nrow(d),
               r = unname(test$estimate), p = test$p.value,
               stringsAsFactors = FALSE)
  }))
  stats <- stats[order(stats$Timepoint, match(stats$group, names(EXERCISE_COLORS))), ]
  stats$label <- sprintf("%s: r = %.3f, p = %.3g", stats$group, stats$r, stats$p)
  rownames(stats) <- NULL
  stats
}

