#!/usr/bin/env Rscript
# video/choose_your_own_adventure.R — draws the start-node choices (called by choose_your_own_adventure.py).
# Arguments: <start_menu_edges.csv> <output png>. One column per story group; each choice is a coloured dot (colour =
# story, validated palette: blue / orange / aqua / violet, every dot also labelled), its first-step neighbours small
# grey dots around it with small grey names, grey links. Descriptive title only (no claims in the figure).
suppressMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
E <- fread(args[1])
GROUPS <- c("T2D muscle", "T2D blood", "T2D + ageing blood", "Ageing blood")
COL <- c("T2D muscle" = "#2a78d6", "T2D blood" = "#eb6834", "T2D + ageing blood" = "#1baf7a", "Ageing blood" = "#4a3aa7")
E[, group := factor(group, levels = GROUPS)]

# layout: the groups present, left to right; each group in columns of up to 3 choices, top to bottom in menu order
ch <- unique(E[, .(choice, group)])
ch[, group := droplevels(group)]
ch[, k := seq_len(.N), by = group]
ch[, ncol := as.integer(ceiling(.N / 3)), by = group]
starts_x <- cumsum(c(0, head(ch[, .(w = max(ncol)), by = group]$w, -1) * 1.1 + 0.4))
names(starts_x) <- as.character(unique(ch$group))
ch[, `:=`(x = starts_x[as.character(group)] + ((k - 1) %/% 3) * 1.1, y = -((k - 1) %% 3) * 1.25)]
R <- 0.36
N <- merge(E, ch[, .(choice, x0 = x, y0 = y)], by = "choice")
N[, `:=`(i = seq_len(.N), n = .N), by = choice]
# neighbours fan out over the UPPER half only (the name sits below the dot): one straight up, more spread 15-165 deg
N[, a := ifelse(n == 1, pi / 2, pi * (15 + 150 * (i - 1) / pmax(n - 1, 1)) / 180)][, `:=`(x1 = x0 + R * cos(a), y1 = y0 + R * sin(a))]
# one colour per choice (Vidal: red, blue, green, yellow, orange, purple; no story heading), placed so that every pair
# of grid neighbours is easy to tell apart (checked with the dataviz palette validator over all 720 placements: worst
# neighbour pair CVD dE 15.3, normal-vision dE 20.8); grid order = row by row. Every dot also carries its name.
HUES <- c(red = "#e34948", yellow = "#eda100", blue = "#2a78d6", green = "#008300", orange = "#eb6834", purple = "#4a3aa7")
setorder(ch, -y, x); ch[, dot := unname(HUES)[(seq_len(.N) - 1) %% length(HUES) + 1]]

p <- ggplot() +
  geom_segment(data = N, aes(x0, y0, xend = x1, yend = y1), colour = "#c9c8c3", linewidth = 0.5) +
  geom_point(data = N, aes(x1, y1), colour = "#9a9993", size = 1.6) +
  geom_text(data = N, aes(x1, y1, label = neighbour), colour = "#6e6d68", size = 2.1,
            vjust = -0.9, hjust = ifelse(cos(N$a) > 0.3, 0.2, ifelse(cos(N$a) < -0.3, 0.8, 0.5))) +
  geom_point(data = ch, aes(x, y), colour = ch$dot, size = 6.5) +
  geom_point(data = ch, aes(x, y), colour = "#fcfcfb", size = 6.5, shape = 1, stroke = 0.9) +   # 2px surface ring
  geom_label(data = ch, aes(x, y - 0.13, label = choice), fill = "#fcfcfb", colour = "#0b0b0b", label.size = 0,
             fontface = "bold", size = 3.6, vjust = 1, label.padding = unit(0.08, "lines")) +

  coord_equal(clip = "off") +
  labs(title = "Choose Your Own Adventure") +                     # title only: no hint of what the choice is for
  theme_void(base_size = 12) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        plot.title = element_text(face = "bold", size = 22, colour = "#0b0b0b", hjust = 0.5, margin = margin(b = 16)),
        plot.subtitle = element_text(size = 10.5, colour = "#52514e", hjust = 0.5, margin = margin(b = 14)),
        legend.position = "bottom", legend.text = element_text(colour = "#0b0b0b"), legend.title = element_text(colour = "#52514e"),
        plot.margin = margin(18, 24, 14, 24))
ggsave(args[2], p, width = max(6.5, 2.2 * (max(ch$x) - min(ch$x) + 2.2)), height = 6.4, dpi = 200, bg = "#fcfcfb")
