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

# layout: groups left to right; T2D muscle (the largest) in two sub-columns; choices top to bottom in menu order
ch <- unique(E[, .(choice, group)])
ch[, k := seq_len(.N), by = group]
ch[, ncol := ifelse(group == "T2D muscle", 2L, 1L)]
gx <- c("T2D muscle" = 0, "T2D blood" = 2.2, "T2D + ageing blood" = 3.4, "Ageing blood" = 4.6)
ch[, `:=`(x = gx[as.character(group)] + ((k - 1) %% ncol) * 1.1, y = -((k - 1) %/% ncol) * 1.25)]
R <- 0.36
N <- merge(E, ch[, .(choice, x0 = x, y0 = y)], by = "choice")
N[, `:=`(i = seq_len(.N), n = .N), by = choice]
# neighbours fan out over the UPPER half only (the name sits below the dot): one straight up, more spread 15-165 deg
N[, a := ifelse(n == 1, pi / 2, pi * (15 + 150 * (i - 1) / pmax(n - 1, 1)) / 180)][, `:=`(x1 = x0 + R * cos(a), y1 = y0 + R * sin(a))]
heads <- ch[, .(x = min(x) + (max(x) - min(x)) / 2, y = 0.62), by = group]

p <- ggplot() +
  geom_segment(data = N, aes(x0, y0, xend = x1, yend = y1), colour = "#c9c8c3", linewidth = 0.5) +
  geom_point(data = N, aes(x1, y1), colour = "#9a9993", size = 1.6) +
  geom_text(data = N, aes(x1, y1, label = neighbour), colour = "#6e6d68", size = 2.1,
            vjust = -0.9, hjust = ifelse(cos(N$a) > 0.3, 0.2, ifelse(cos(N$a) < -0.3, 0.8, 0.5))) +
  geom_point(data = ch, aes(x, y, colour = group), size = 6.5) +
  geom_point(data = ch, aes(x, y), colour = "#fcfcfb", size = 6.5, shape = 1, stroke = 0.9) +   # 2px surface ring
  geom_label(data = ch, aes(x, y - 0.13, label = choice), fill = "#fcfcfb", colour = "#0b0b0b", label.size = 0,
             fontface = "bold", size = 3.6, vjust = 1, label.padding = unit(0.08, "lines")) +
  geom_text(data = heads, aes(x, y, label = group, colour = group), fontface = "bold", size = 4.2, show.legend = FALSE) +
  scale_colour_manual(values = COL, breaks = GROUPS, name = "Story") +
  coord_equal(clip = "off") +
  labs(title = "Choose Your Own Adventure",
       subtitle = "Start nodes for the Team 2-PAC song: the strongest markers of our three stories. Grey = where the walk can go first.") +
  theme_void(base_size = 12) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        plot.title = element_text(face = "bold", size = 22, colour = "#0b0b0b", hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10.5, colour = "#52514e", hjust = 0.5, margin = margin(b = 14)),
        legend.position = "bottom", legend.text = element_text(colour = "#0b0b0b"), legend.title = element_text(colour = "#52514e"),
        plot.margin = margin(18, 24, 14, 24))
ggsave(args[2], p, width = 11, height = 6.4, dpi = 200, bg = "#fcfcfb")
