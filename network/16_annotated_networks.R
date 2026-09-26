#!/usr/bin/env Rscript
# =====================================================================================================
# 16_annotated_networks.R — STEP 16: THE JOINT DIFFERENCE NETWORK (15b) ANNOTATED WITH PHOSPHO / GLYCOSYLATION
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   A first look at what the next data layers would add to the joint network: which proteins in figure 15b
#   (edges = endurance minus resistance) also have phosphosites that respond to exercise (MoTrPAC), and
#   which are known to be glycosylated (GlyGen). A sample visualisation for choosing what to integrate.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Figure 15b is redrawn unchanged (same class-grouped layout, class bubbles, edge colours and widths =
#   w_EE - w_RE, node sizes, metabolite triangles) and only the protein circles are recoloured:
#     16a, MoTrPAC phosphoproteomics (measured; muscle 0.5 / 4 / 24 h and adipose 4 h; both arms):
#         not measured · measured, no site responds · a site responds after endurance only · after
#         resistance only · after both ("responds" = adj. p < 0.05 in that arm's exercise-vs-control
#         contrast at any time point, either tissue). Labels: the proteins with the most responding sites,
#         with "responding / measured" site counts.
#     16b, GlyGen glycosylation (prior knowledge; release 2.11.1; nothing measured in MoTrPAC):
#         no glycosylation record · glycosylated, site unknown (protein-level evidence) · N-linked sites
#         only · O-linked sites only (incl. O-GlcNAc) · both N- and O-linked sites. Labels: the proteins
#         with the most glycosylation sites, with site and glycan-structure counts.
#     16c, both together + site-level crosstalk: fill = the 16a phospho category, ring = the 16b glycosylation
#         category, and a diamond on proteins where a measured phosphosite is the SAME residue as a known
#         O-glycosylation (mostly O-GlcNAc) site (same canonical protein, position, residue; multi-site
#         features contribute each residue); gold diamond = such a residue responds to exercise.
#
# HOW TO RUN
#   After step 15 and the inventory (network/inventory/glygen_protein_inventory.py, then
#   glygen_motrpac_inventory.R):   Rscript network/16_annotated_networks.R   (about 20 seconds)
#
# DATA AND PROVENANCE
#   Network: steps 14 / 15. Phospho: MotrpacHumanPreSuspensionAnalysis v0.2.4 (MUSCLE/ADIPOSE_PROT_PH_DA,
#   HUMAN_FEATURE_TO_GENE). Glycosylation: GlyGen release 2.11.1 via the inventory (protein_inventory.csv).
#
# TECH STACK:  R 4.4; data.table, ggplot2, ggrepel, ggforce, ggnewscale (second fill / colour scales in 16c).
#
# INPUTS
#   $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 15_class_layout.csv, 02_nodes_string.csv
#   $HACK_OUT/inventory/protein_inventory.csv   (gly_sites_N, gly_sites_O, gly_protein_level_no_site, glycans_at_sites)
#   $HACK_OUT/inventory/glygen_glycosites.csv   (glycosylation site positions, for 16c)
# OUTPUTS
#   $HACK_FIG/16a_joint_edge_difference_phospho.png, 16b_joint_edge_difference_glycosylation.png,
#            16c_joint_edge_difference_phospho_glyco_crosstalk.png (never in the repo)
#   $HACK_OUT/16_crosstalk_sites.csv      per residue measured as a MoTrPAC phosphosite AND listed by GlyGen as an
#                                         O-glycosylation site: protein, residue, features, responds (EE / RE / both / no),
#                                         glycosylation subtype, evidence and source
#   $HACK_OUT/16_protein_annotation.csv   per network protein: phospho sites measured / responding (EE, RE),
#                                         glycosylation category and counts
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   304 proteins and 60 metabolites drawn (as 15b). The script stops if the layout, nodes or edges differ
#   from steps 14 / 15, or if a protein has no annotation row.
#
# KNOWN LIMITS
#   Node colour is a protein-level summary: site-level detail (which site, which tissue and time) is in
#   16_protein_annotation.csv and the MoTrPAC tables. Glycosylation is database knowledge, not measured in
#   this study. Differences between arms are descriptive (not tested).
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(ggforce) })

# Folders (override with HACK_OUT / HACK_FIG); figures and data never go into the repo.
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Same seed, label count and styling as step 15.
SEED <- 20260926; NLAB <- 12; ALPHA <- 0.05
INK <- "#1A1A1A"
TYPE_LTY <- c("protein - protein" = "solid", "metabolite - metabolite" = "42", "metabolite - protein" = "11")
SHAPES <- c(protein = 21, metabolite = 24)
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"

# ---- network and layout (as figure 15b) --------------------------------------------------------------------
E <- fread(file.path(OUT, "14_joint_edges.csv"))[, .(a = node_a, b = node_b, edge_type = factor(edge_type, levels = names(TYPE_LTY)), w_diff)]
N <- fread(file.path(OUT, "14_joint_nodes.csv"))[, .(node, node_type, class, strength_EE, strength_RE)]
N <- fread(file.path(OUT, "15_class_layout.csv"))[, .(node, x, y)][N, on = "node"]
stopifnot(!anyNA(N$x), all(c(E$a, E$b) %in% N$node))
N[, abs_delta := abs(strength_EE - strength_RE)]

# ---- protein annotation: MoTrPAC phosphoproteomics ----------------------------------------------------------
genes <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol)]
data("HUMAN_FEATURE_TO_GENE", package = "MotrpacHumanPreSuspensionAnalysis")
f2g <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay == "prot-ph", .(feature_id = as.character(feature_id), entrez_gene = as.character(entrez_gene))])
# per feature: responds after EE / RE (any time point) in each tissue, then across tissues
ph <- rbindlist(lapply(c("MUSCLE", "ADIPOSE"), function(t) {
  data(list = paste0(t, "_PROT_PH_DA"), package = "MotrpacHumanPreSuspensionAnalysis")
  d <- as.data.table(get(paste0(t, "_PROT_PH_DA")))[contrast_type == "exercise_with_controls"]
  d[, .(sig_EE = any(adj_p_value < ALPHA & grepl("^Endur", contrast_short)), sig_RE = any(adj_p_value < ALPHA & grepl("^Resist", contrast_short))),
    by = .(feature_id = as.character(feature_id))]
}))[, .(sig_EE = any(sig_EE), sig_RE = any(sig_RE)), by = feature_id]
ph <- f2g[ph, on = "feature_id", nomatch = 0]
# per protein: sites measured, responding in EE, in RE, in either
php <- ph[, .(ph_measured = .N, ph_resp_EE = sum(sig_EE), ph_resp_RE = sum(sig_RE), ph_resp_any = sum(sig_EE | sig_RE)), by = entrez_gene]
A <- php[genes, on = "entrez_gene"]
for (v in c("ph_measured", "ph_resp_EE", "ph_resp_RE", "ph_resp_any")) set(A, which(is.na(A[[v]])), v, 0L)

# ---- protein annotation: GlyGen glycosylation --------------------------------------------------------------
inv <- file.path(OUT, "inventory", "protein_inventory.csv")
if (!file.exists(inv)) stop("run network/inventory/glygen_protein_inventory.py and glygen_motrpac_inventory.R first (", inv, " missing)")
gl <- fread(inv, colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gly_N = gly_sites_N, gly_O = gly_sites_O, gly_sites = gly_sites_unique,
                                                                   gly_protein_level = gly_protein_level_no_site, glycans = glycans_at_sites)]
A <- gl[A, on = "entrez_gene"]
stopifnot(!anyNA(A$gly_sites))

# ---- categories -----------------------------------------------------------------------------------------------
PH_LEV <- c("not measured", "measured, no site responds", "a site responds after endurance only",
            "a site responds after resistance only", "sites respond after both")
PH_COL <- setNames(c("grey85", "white", "#D6604D", "#4393C3", "#762A83"), PH_LEV)
A[, ph_cat := fcase(ph_measured == 0, PH_LEV[1], ph_resp_any == 0, PH_LEV[2], ph_resp_EE > 0 & ph_resp_RE == 0, PH_LEV[3],
                    ph_resp_RE > 0 & ph_resp_EE == 0, PH_LEV[4], default = PH_LEV[5])]
GL_LEV <- c("no glycosylation record", "glycosylated, site unknown (protein-level)", "N-linked sites only",
            "O-linked sites only (incl. O-GlcNAc)", "N- and O-linked sites")
GL_COL <- setNames(c("grey85", "#E6F5D0", "#1B9E77", "#D95F02", "#7570B3"), GL_LEV)
A[, gl_cat := fcase(gly_N > 0 & gly_O > 0, GL_LEV[5], gly_N > 0, GL_LEV[3], gly_O > 0, GL_LEV[4], gly_protein_level > 0, GL_LEV[2], default = GL_LEV[1])]
# attach to the drawn proteins and save
N <- A[, !"entrez_gene"][N, on = c(gene_symbol = "node")]; setnames(N, "gene_symbol", "node")
stopifnot(!anyNA(N[node_type == "protein", ph_cat]), !anyNA(N[node_type == "protein", gl_cat]))
fwrite(N[node_type == "protein", .(protein = node, ph_measured, ph_resp_EE, ph_resp_RE, ph_resp_any, ph_cat, gly_sites, gly_N, gly_O,
                                   gly_protein_level, glycans, gl_cat)], file.path(OUT, "16_protein_annotation.csv"))

# ---- drawing (figure 15b + recoloured proteins) -----------------------------------------------------------------
theme_net <- function(base = 8) theme_classic(base_size = base) %+replace% theme(
  text = element_text(colour = INK, size = base), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
  axis.title = element_blank(), legend.position = "bottom", legend.box = "vertical",
  legend.title = element_text(colour = INK, size = base * 0.8), legend.text = element_text(colour = INK, size = base * 0.74),
  legend.background = element_blank(), legend.key = element_blank(),
  plot.title = element_text(colour = INK, face = "bold", size = base * 1.18, hjust = 0, margin = margin(b = base * 0.55)),
  plot.caption = element_text(colour = "grey35", size = base * 0.72, hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
# edge segments at the grouped layout, smallest differences drawn first; colour limit at the 95th percentile
D <- copy(E)[, `:=`(x = N$x[match(a, N$node)], y = N$y[match(a, N$node)], xend = N$x[match(b, N$node)], yend = N$y[match(b, N$node)])]
lim_d <- as.numeric(quantile(abs(D$w_diff), 0.95)); D <- D[order(abs(w_diff))]
# class bubbles and names (as 15b)
cn <- N[node_type == "metabolite", .(x = mean(x), y = max(y) + 0.035), by = class]
base_plot <- function() ggplot() +
  geom_mark_hull(data = N[node_type == "metabolite"], aes(x, y, group = class), colour = "grey45", fill = "grey60", alpha = 0.08,
                 linewidth = 0.25, expand = unit(2.2, "mm"), radius = unit(2, "mm"), concavity = 3) +
  geom_segment(data = D, aes(x, y, xend = xend, yend = yend, colour = w_diff, linewidth = abs(w_diff), linetype = edge_type), lineend = "round") +
  # metabolites: grey triangles, as 15b
  geom_point(data = N[node_type == "metabolite"], aes(x, y, size = abs_delta, shape = node_type), fill = "grey80", colour = "grey25", stroke = 0.22) +
  geom_text_repel(data = cn, aes(x, y, label = class), size = 2.4, colour = "grey30", fontface = "bold.italic", box.padding = 0.15,
                  min.segment.length = 0.3, segment.size = 0.15, segment.colour = "grey50", max.overlaps = Inf, seed = SEED) +
  scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-lim_d, lim_d), oob = scales::squish,
                         breaks = c(-lim_d, 0, lim_d), labels = c("higher in resistance\n(w_RE > w_EE)", "same", "higher in endurance\n(w_EE > w_RE)"),
                         name = "edge difference  w_EE − w_RE",
                         guide = guide_colourbar(order = 10, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5)) +
  scale_linetype_manual(values = TYPE_LTY, name = "edge type") + scale_linewidth(range = c(0.1, 1.6), guide = "none") +
  scale_size(range = c(0.6, 4), name = "node strength difference (absolute)") + scale_shape_manual(values = SHAPES, name = "node") +
  coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.05), clip = "off") + theme_net()

# one annotated figure: proteins filled by category, labels for the top proteins by a count
annotated <- function(cat_col, levels, cols, legend_name, lab_rank, lab_text, title, caption, file) {
  P <- N[node_type == "protein"]
  # legend entries with the number of drawn proteins in each category
  n_cat <- P[, .N, by = c(cat_col)]; lbl <- sprintf("%s (%d)", levels, n_cat$N[match(levels, n_cat[[cat_col]])]); lbl <- sub("\\(NA\\)", "(0)", lbl)
  P[, cat := factor(get(cat_col), levels = levels)]
  P[, rank_value := lab_rank(P)]
  L <- P[rank_value > 0][order(-rank_value)][seq_len(min(NLAB, .N))]
  L[, lab := lab_text(L)]
  q <- base_plot() +
    geom_point(data = P, aes(x, y, size = abs_delta, shape = node_type, fill = cat), colour = "grey20", stroke = 0.3) +
    scale_fill_manual(values = cols, labels = setNames(lbl, levels), drop = FALSE, name = legend_name) +
    geom_text_repel(data = L, aes(x, y, label = lab), size = 2.1, colour = "grey10", min.segment.length = 0.1, segment.size = 0.15,
                    box.padding = 0.3, max.overlaps = Inf, seed = SEED) +
    guides(fill = guide_legend(order = 1, ncol = 3, override.aes = list(shape = 21, size = 3)),
           colour = guide_colourbar(order = 2, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5),
           linetype = guide_legend(order = 3, override.aes = list(colour = "grey30", linewidth = 0.6)),
           shape = guide_legend(order = 4, override.aes = list(fill = "grey80", size = 2.6)), size = guide_legend(order = 5)) +
    labs(title = title, caption = caption)
  ggsave(file.path(FIG, file), q, width = 11, height = 7.8, dpi = 300, bg = "white")
  message("-> ", file.path(FIG, file))
}

# 16a: MoTrPAC phosphoproteomics
annotated("ph_cat", PH_LEV, PH_COL, "protein: MoTrPAC phosphosites",
          lab_rank = function(d) d$ph_resp_any + d$ph_measured / 1000,
          lab_text = function(d) sprintf("%s  %d/%d", d$node, d$ph_resp_any, d$ph_measured),
          title = "Joint protein-metabolite network, endurance minus resistance edge weights, proteins annotated with MoTrPAC phosphoproteomics",
          caption = paste0("Phosphosites measured in muscle (0.5, 4, 24 h) and adipose (4 h); responds = adj. p < 0.05 in that arm's exercise-vs-control contrast at any time point. ",
                           "Labels: proteins with the most responding sites, responding / measured."),
          file = "16a_joint_edge_difference_phospho.png")
# 16b: GlyGen glycosylation
annotated("gl_cat", GL_LEV, GL_COL, "protein: GlyGen glycosylation",
          lab_rank = function(d) d$gly_sites + d$glycans / 1000,
          lab_text = function(d) sprintf("%s  %d sites, %d glycans", d$node, d$gly_sites, d$glycans),
          title = "Joint protein-metabolite network, endurance minus resistance edge weights, proteins annotated with GlyGen glycosylation",
          caption = paste0("GlyGen release 2.11.1 (database knowledge, not measured in this study). Sites with a known position; 'site unknown' = protein-level evidence, mostly the O-GlcNAc Database. ",
                           "Labels: proteins with the most glycosylation sites, with glycan structures observed at those sites."),
          file = "16b_joint_edge_difference_glycosylation.png")
# ---- 16c: both layers + site-level crosstalk -------------------------------------------------------------------
# Crosstalk residue = a serine / threonine / tyrosine that MoTrPAC measured as a phosphosite AND GlyGen lists as an
# O-glycosylation site (mostly O-GlcNAc) on the same canonical protein, position and residue. The two modifications
# compete for the same hydroxyl group, so these are candidate phospho / O-GlcNAc switch sites.
gsf <- file.path(OUT, "inventory", "glygen_glycosites.csv")
if (!file.exists(gsf)) stop("re-run network/inventory/glygen_protein_inventory.py (", gsf, " missing)")
og <- fread(gsf)[type == "O-linked" & residue %in% c("Ser", "Thr", "Tyr")]
og <- og[, .(gly_subtypes = paste(sort(unique(subtype)), collapse = ";"), gly_evidence = paste(sort(unique(category)), collapse = ";"),
             gly_sources = paste(sort(unique(unlist(strsplit(source, ";")))), collapse = ";")),
         by = .(glygen_ac, base = sub("-.*$", "", glygen_ac), position, res1 = substr(residue, 1, 1))]
# MoTrPAC features split into their sites (a multi-site feature contributes each of its residues)
fs <- ph[, .(feature_id, entrez_gene, sig_EE, sig_RE, acc = sub("_.*$", "", feature_id), sites = sub("^[^_]*_", "", feature_id))]
fs <- fs[, .(res1 = regmatches(sites, gregexpr("[STY]", sites))[[1]], position = as.integer(regmatches(sites, gregexpr("[0-9]+", sites))[[1]]),
             n_sites = lengths(regmatches(sites, gregexpr("[STY][0-9]+", sites)))), by = .(feature_id, entrez_gene, sig_EE, sig_RE, acc)]
# match: plain accessions (UniProt canonical sequence) on the GlyGen canonical base; isoform accessions only if GlyGen's canonical is that isoform
xt <- rbind(merge(fs[!grepl("-", acc)], og, by.x = c("acc", "position", "res1"), by.y = c("base", "position", "res1")),
            merge(fs[grepl("-", acc)], og[, !"base"], by.x = c("acc", "position", "res1"), by.y = c("glygen_ac", "position", "res1")), fill = TRUE)
xt <- genes[xt, on = "entrez_gene"]
# one row per residue: features measuring it, and whether any of them responds (EE, RE); single-site evidence flagged
XT <- xt[, .(features = paste(sort(unique(feature_id)), collapse = ";"), n_features = uniqueN(feature_id),
             single_site_feature = any(n_sites == 1), responds_EE = any(sig_EE), responds_RE = any(sig_RE),
             gly_subtypes = gly_subtypes[1], gly_evidence = gly_evidence[1], gly_sources = gly_sources[1]),
         by = .(protein = gene_symbol, residue = paste0(res1, position), position)][order(protein, position)]
XT[, responds := fcase(responds_EE & responds_RE, "both", responds_EE, "endurance", responds_RE, "resistance", default = "no")]
fwrite(XT[, !"position"], file.path(OUT, "16_crosstalk_sites.csv"))
# per protein: crosstalk residues and how many respond
xp <- XT[, .(xt_n = .N, xt_resp = sum(responds != "no"), xt_label = paste0(protein[1], "  ", paste(residue, collapse = ","))), by = protein]
xp[, xt_label := fifelse(xt_n > 4, sprintf("%s  %d residues", protein, xt_n), xt_label)]
xp[, xt_label := sprintf("%s (%d/%d respond)", xt_label, xt_resp, xt_n)]
P3 <- xp[N[node_type == "protein"], on = c(protein = "node")]; setnames(P3, "protein", "node")
P3[is.na(xt_n), `:=`(xt_n = 0L, xt_resp = 0L)]
P3[, `:=`(ph_f = factor(ph_cat, levels = PH_LEV), gl_f = factor(gl_cat, levels = GL_LEV))]
# ring colours for glycosylation (darker than the 16b fills so a ring stays visible)
GL_RING <- setNames(c("grey55", "#A6D96A", "#1B9E77", "#D95F02", "#7570B3"), GL_LEV)
nlab <- function(df, col, lev) { n <- df[, .N, by = c(col)]; sub("\\(NA\\)", "(0)", sprintf("%s (%d)", lev, n$N[match(lev, n[[col]])])) }
# the crosstalk marker: a diamond just above-right of the protein, gold if a crosstalk residue responds
MK <- P3[xt_n > 0][, `:=`(mx = x + 0.009, my = y + 0.014, mk = factor(fifelse(xt_resp > 0, "a shared residue responds (either arm)", "shared residue(s), none respond"),
                                                                   levels = c("a shared residue responds (either arm)", "shared residue(s), none respond")))]
q <- base_plot() +
  ggnewscale::new_scale_colour() +
  geom_point(data = P3, aes(x, y, size = abs_delta, shape = node_type, fill = ph_f, colour = gl_f), stroke = 1.05) +
  scale_fill_manual(values = PH_COL, labels = setNames(nlab(P3, "ph_cat", PH_LEV), PH_LEV), drop = FALSE, name = "fill: MoTrPAC phosphosites",
                    guide = guide_legend(order = 1, ncol = 3, override.aes = list(shape = 21, size = 3, colour = "grey20", stroke = 0.3))) +
  scale_colour_manual(values = GL_RING, labels = setNames(nlab(P3, "gl_cat", GL_LEV), GL_LEV), drop = FALSE, name = "ring: GlyGen glycosylation",
                      guide = guide_legend(order = 2, ncol = 3, override.aes = list(shape = 21, size = 3, fill = "white", stroke = 1.3))) +
  ggnewscale::new_scale_fill() +
  geom_point(data = MK, aes(mx, my, fill = mk), shape = 23, size = 2.1, colour = "grey10", stroke = 0.35) +
  scale_fill_manual(values = c("a shared residue responds (either arm)" = "#FFD700", "shared residue(s), none respond" = "white"), drop = FALSE,
                    name = "diamond: phosphosite = O-glycosylation site", guide = guide_legend(order = 3, override.aes = list(size = 3))) +
  geom_text_repel(data = MK, aes(x, y, label = xt_label), size = 2.05, colour = "grey10", min.segment.length = 0.1, segment.size = 0.15,
                  box.padding = 0.35, max.overlaps = Inf, seed = SEED) +
  guides(linetype = guide_legend(order = 11, override.aes = list(colour = "grey30", linewidth = 0.6)),
         shape = guide_legend(order = 12, override.aes = list(fill = "grey80", colour = "grey25", stroke = 0.3, size = 2.6)), size = guide_legend(order = 13)) +
  labs(title = "Joint protein-metabolite network, endurance minus resistance edge weights: MoTrPAC phosphosites, GlyGen glycosylation and shared residues",
       caption = paste0("Fill: phosphosites (muscle 0.5/4/24 h, adipose 4 h; responds = adj. p < 0.05 vs control). Ring: GlyGen glycosylation (release 2.11.1, database knowledge). ",
                        "Diamond: a measured phosphosite is also a known O-glycosylation (mostly O-GlcNAc) site on the same residue; labels list those residues."))
ggsave(file.path(FIG, "16c_joint_edge_difference_phospho_glyco_crosstalk.png"), q, width = 11.5, height = 8.6, dpi = 300, bg = "white")
message("-> ", file.path(FIG, "16c_joint_edge_difference_phospho_glyco_crosstalk.png"))
cat(sprintf("crosstalk residues: %d on %d proteins (%d drawn in the network); responding: %d\n", nrow(XT), uniqueN(XT$protein),
            sum(P3$xt_n > 0), sum(XT$responds != "no")))

# Show the category counts.
print(N[node_type == "protein", .N, by = ph_cat][order(-N)]); print(N[node_type == "protein", .N, by = gl_cat][order(-N)])
