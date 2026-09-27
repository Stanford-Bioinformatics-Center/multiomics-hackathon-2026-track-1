# =====================================================================================================
# R/figure_style.R — one figure style for every final figure (modelled on the MoTrPAC landscape papers)
# =====================================================================================================
# WHAT IT GIVES: theme_motrpac() (small sans-serif type, thin black axes, no minor grid, light strips),
# the project palettes (arms, disease direction, edge types), and save_figure() (PNG at 300 dpi + PDF).
# Panels are tagged a, b, c ... (bold, lower case) with patchwork::plot_annotation(tag_levels = "a"), and the
# figure title follows the "Fig. N | short descriptive title" pattern. Titles stay descriptive; claims and
# method detail go in the captions (network/docs/FIGURES.md), not inside the figure.
# HOW TO USE:  source(file.path(HERE, "R", "figure_style.R"))
# =====================================================================================================
suppressMessages({ library(ggplot2) })

# Colours used everywhere, so a colour means the same thing in every figure.
ARM_COL <- c(endurance = "#D7301F", resistance = "#2B8CBE", both = "#7B3294", neither = "#BDBDBD")   # arms (PTM tags use the same hues)
DIR_COL <- c(lower = "#5E3C99", higher = "#E66100")                                                     # lower / higher in the disease (or with age)
TYPE_COL <- c("protein - protein" = "#737373", "metabolite - metabolite" = "#1B7837", "metabolite - protein" = "#8C510A")
DIV_PAL <- c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")                                        # diverging fill for responses / weights

#' The figure theme.
#' @param base `numeric(1)` base font size in points (default 8, as in journal figures at print size).
#' @return a ggplot2 theme.
theme_motrpac <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = base + 1, hjust = 0),
          plot.subtitle = element_text(size = base, colour = "grey25"),
          plot.tag = element_text(face = "bold", size = base + 4),
          axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          axis.text = element_text(size = base - 1, colour = "black"), axis.title = element_text(size = base),
          strip.background = element_rect(fill = "grey92", colour = NA), strip.text = element_text(size = base, face = "bold"),
          legend.title = element_text(size = base - 1, face = "bold"), legend.text = element_text(size = base - 1),
          legend.key.size = grid::unit(3, "mm"), panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"))
}

#' A blank theme for network panels (no axes), same typography.
theme_motrpac_void <- function(base = 8) {
  theme_void(base_size = base, base_family = "Helvetica") +
    theme(plot.title = element_text(face = "bold", size = base + 1, hjust = 0), plot.tag = element_text(face = "bold", size = base + 4),
          legend.title = element_text(size = base - 1, face = "bold"), legend.text = element_text(size = base - 1), legend.key.size = grid::unit(3, "mm"))
}

#' Save a figure as PNG (300 dpi) and PDF with the same name.
#' @param plot a ggplot / patchwork object.
#' @param path `character(1)` output path without extension.
#' @param width,height `numeric(1)` size in mm (journal figure widths: 89 single, 183 double column).
save_figure <- function(plot, path, width = 183, height = 120) {
  ggsave(paste0(path, ".png"), plot, width = width, height = height, units = "mm", dpi = 300, bg = "white")
  # (the base pdf device: deterministic apart from its date stamps, which fix_pdf_dates() pins; cairo_pdf varies run to run)
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height, units = "mm", device = grDevices::pdf, bg = "white", useDingbats = FALSE)
  fix_pdf_dates(paste0(path, ".pdf"))
  message("-> ", path, ".png / .pdf")
}

#' Make a PDF byte-identical across reruns: the PDF device writes the current time into /CreationDate and
#' /ModDate ("D:YYYYMMDDHHMMSS..."); overwrite those 14 digits in place with the fixed build date
#' (SOURCE_DATE_EPOCH, as run_all.sh sets it). Same length, so the PDF's internal byte offsets stay valid.
#' @param file `character(1)` path of the PDF.
fix_pdf_dates <- function(file) {
  stamp <- format(as.POSIXct(as.numeric(Sys.getenv("SOURCE_DATE_EPOCH", "1790380800")), origin = "1970-01-01", tz = "UTC"), "%Y%m%d%H%M%S")
  b <- readBin(file, "raw", file.info(file)$size); key <- charToRaw("Date (D:"); new <- charToRaw(stamp)
  hits <- which(vapply(seq_len(length(b) - length(key) - 14), function(i) all(b[i:(i + length(key) - 1)] == key), logical(1)))
  for (h in hits) b[(h + length(key)):(h + length(key) + 13)] <- new
  writeBin(b, file)
}
