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
#     16c, both together + site-level crosstalk: fill = the 16a phospho category, ring = glycosylated or not
#         (any type, site known or protein-level; the types are in 16b), and a diamond on proteins where a measured phosphosite is the SAME residue as a known
#         O-glycosylation (mostly O-GlcNAc) site (same canonical protein, position, residue; multi-site
#         features contribute each residue); gold diamond = such a residue responds to exercise.
#
# HOW TO RUN
#   After step 15 and the inventory (network/inventory/glygen_protein_inventory.py, then
#   glygen_motrpac_inventory.R):   Rscript network/16_annotated_networks.R   (about 20 seconds)
#
# DATA AND PROVENANCE
#   Network: steps 14 / 15. Phospho: MotrpacHumanPreSuspensionAnalysis v0.2.4 (MUSCLE/ADIPOSE_PROT_PH_DA,
#   HUMAN_FEATURE_TO_GENE). Glycosylation: UniProt sites via the team's mnet resource (resources/mo_annotation) + GlyGen
#   release 2.11.1 protein-level evidence and glycans (inventory); crosstalk through mnet's MoTrPAC site bridge.
#
# TECH STACK:  R 4.4; data.table, ggplot2, ggrepel, ggforce, ggnewscale (second fill / colour scales in 16c).
#
# INPUTS
#   $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 15_class_layout.csv, 02_nodes_string.csv
#   $HACK_OUT/inventory/protein_inventory.csv   (gly_sites_N, gly_sites_O, gly_protein_level_no_site, glycans_at_sites)
#   $HACK_OUT/inventory/glygen_glycosites.csv   (GlyGen O-glycosylation sites, for 16c)
#   resources/mo_annotation (MNET_DIR): proteins_ptm.csv, glycosites.csv, motrpac_feature_site_map.csv (mnet PTM)
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
# Create the figure folder if needed.
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
# Same seed, label count and styling as step 15.
SEED <- 20260926; NLAB <- 12; ALPHA <- 0.05
# (ALPHA = the adjusted p-value below which a phosphosite counts as responding)
# Text colour, edge line types by type, node shapes (protein circle, metabolite triangle), difference colours.
INK <- "#1A1A1A"
TYPE_LTY <- c("protein - protein" = "solid", "metabolite - metabolite" = "42", "metabolite - protein" = "11")
SHAPES <- c(protein = 21, metabolite = 24)
COL_RE <- "#2166AC"; COL_SAME <- "grey85"; COL_EE <- "#B2182B"

# ---- network and layout (as figure 15b) --------------------------------------------------------------------
# Edges from step 14: the two end nodes, the edge type (in a fixed order) and the endurance-minus-resistance difference.
E <- fread(file.path(OUT, "14_joint_edges.csv"))[, .(a = node_a, b = node_b, edge_type = factor(edge_type, levels = names(TYPE_LTY)), w_diff)]
# Nodes from step 14: type, metabolite class and strength in each arm.
N <- fread(file.path(OUT, "14_joint_nodes.csv"))[, .(node, node_type, class, strength_EE, strength_RE)]
# Add each node's position in the step 15 class-grouped layout (so the picture matches figure 15b).
N <- fread(file.path(OUT, "15_class_layout.csv"))[, .(node, x, y)][N, on = "node"]
# Safety check: every node has a position and every edge end is a node.
stopifnot(!anyNA(N$x), all(c(E$a, E$b) %in% N$node))
# Node size in the figures = size of the strength difference between the arms (as in 15b).
N[, abs_delta := abs(strength_EE - strength_RE)]

# ---- protein annotation: MoTrPAC phosphoproteomics ----------------------------------------------------------
# Our 471 genes (Entrez ID and symbol) from step 2.
genes <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gene_symbol)]
# MoTrPAC's table linking each measured feature to its gene; keep the phosphosite features, as text IDs.
data("HUMAN_FEATURE_TO_GENE", package = "MotrpacHumanPreSuspensionAnalysis")
f2g <- unique(as.data.table(HUMAN_FEATURE_TO_GENE)[assay == "prot-ph", .(feature_id = as.character(feature_id), entrez_gene = as.character(entrez_gene))])
# per feature: responds after EE / RE (any time point) in each tissue, then across tissues
ph <- rbindlist(lapply(c("MUSCLE", "ADIPOSE"), function(t) {
  # load this tissue's phosphoproteomics results and keep the exercise-vs-control (delta-delta) contrasts
  data(list = paste0(t, "_PROT_PH_DA"), package = "MotrpacHumanPreSuspensionAnalysis")
  d <- as.data.table(get(paste0(t, "_PROT_PH_DA")))[contrast_type == "exercise_with_controls"]
  # per feature: TRUE if any endurance (or resistance) time point has adj. p < ALPHA
  d[, .(sig_EE = any(adj_p_value < ALPHA & grepl("^Endur", contrast_short)), sig_RE = any(adj_p_value < ALPHA & grepl("^Resist", contrast_short))),
    by = .(feature_id = as.character(feature_id))]
}))[, .(sig_EE = any(sig_EE), sig_RE = any(sig_RE)), by = feature_id]
# attach the gene to each phospho feature (features without one of our genes' IDs are dropped)
ph <- f2g[ph, on = "feature_id", nomatch = 0]
# per protein: sites measured, responding in EE, in RE, in either
# (each MoTrPAC phospho feature counts as one "site" here; a feature can span more than one residue)
php <- ph[, .(ph_measured = .N, ph_resp_EE = sum(sig_EE), ph_resp_RE = sum(sig_RE), ph_resp_any = sum(sig_EE | sig_RE)), by = entrez_gene]
# one row per network gene (genes without phospho data get empty values)
A <- php[genes, on = "entrez_gene"]
# genes without phospho data: counts set to 0
for (v in c("ph_measured", "ph_resp_EE", "ph_resp_RE", "ph_resp_any")) set(A, which(is.na(A[[v]])), v, 0L)

# ---- protein annotation: glycosylation (mnet UniProt sites + GlyGen protein-level evidence and glycans) -------
# Site-level glycosylation from the team's mnet resource (resources/mo_annotation, UniProt glycosites); GlyGen
# (network/inventory) adds protein-level evidence (mostly O-GlcNAc Database, no site) and glycan structures.
# Where the mnet annotation files are (override with MNET_DIR; the default is outside the repo).
MNET_DIR <- Sys.getenv("MNET_DIR", unset = path.expand("~/Desktop/output/hackathon/resources/mo_annotation"))
# The GlyGen protein inventory (written by network/inventory); stop with a message if it is missing.
inv <- file.path(OUT, "inventory", "protein_inventory.csv")
if (!file.exists(inv)) stop("run network/inventory/glygen_protein_inventory.py and glygen_motrpac_inventory.R first (", inv, " missing)")
# Each gene's UniProt accession(s), one row per accession (a gene can have several, ";"-separated).
g_acc <- fread(file.path(OUT, "02_nodes_string.csv"), colClasses = list(character = "entrez_gene"))[, .(acc = unlist(strsplit(uniprot, ";"))), by = entrez_gene]
# mnet per-protein counts (isoform suffix removed from the accession): N-linked sites, O-linked sites (in this
# file these EXCLUDE O-GlcNAc), O-GlcNAc sites and all glycosites.
mptm <- fread(file.path(MNET_DIR, "proteins_ptm.csv"))[, .(acc = sub("-[0-9]+$", "", node_id), gly_N = n_N_linked, gly_O = n_O_linked, gly_OG = n_O_GlcNAc, gly_sites = n_glycosites)]
# match to our genes; if a gene has several accessions, keep the one with the most glycosites
mg <- merge(g_acc, mptm, by = "acc")[order(-gly_sites)][!duplicated(entrez_gene), !"acc"]
# GlyGen: protein-level glycosylation evidence without a site, and glycan structures at sites
gl <- fread(inv, colClasses = list(character = "entrez_gene"))[, .(entrez_gene, gly_protein_level = gly_protein_level_no_site, glycans = glycans_at_sites)]
# combine the mnet site counts with the GlyGen columns (one row per gene in the inventory)
gl <- mg[gl, on = "entrez_gene"]
# genes without mnet glycosylation records: site counts set to 0
for (v in c("gly_N", "gly_O", "gly_OG", "gly_sites")) set(gl, which(is.na(gl[[v]])), v, 0L)
# add the glycosylation columns to the per-gene phospho table
A <- gl[A, on = "entrez_gene"]
# Safety check: every gene has a glycosylation site count.
stopifnot(!anyNA(A$gly_sites))

# ---- categories -----------------------------------------------------------------------------------------------
# Phospho categories and colours (the categories go into the saved table; PH_COL is not used by the current
# figures, which draw tags instead of colouring the protein circles).
PH_LEV <- c("not measured", "measured, no site responds", "a site responds after endurance only",
            "a site responds after resistance only", "sites respond after both")
PH_COL <- setNames(c("grey85", "white", "#D6604D", "#4393C3", "#762A83"), PH_LEV)
# each protein's phospho category: not measured / nothing responds / endurance only / resistance only / both
A[, ph_cat := fcase(ph_measured == 0, PH_LEV[1], ph_resp_any == 0, PH_LEV[2], ph_resp_EE > 0 & ph_resp_RE == 0, PH_LEV[3],
                    ph_resp_RE > 0 & ph_resp_EE == 0, PH_LEV[4], default = PH_LEV[5])]
# Glycosylation categories and colours (again, GL_COL is not used by the current figures).
GL_LEV <- c("no glycosylation record", "glycosylated, site unknown (protein-level)", "N-linked sites only",
            "O-linked sites only (incl. O-GlcNAc)", "N- and O-linked sites")
GL_COL <- setNames(c("grey85", "#E6F5D0", "#1B9E77", "#D95F02", "#7570B3"), GL_LEV)
# each protein's glycosylation category, from its N-linked, O-linked and protein-level counts
A[, gl_cat := fcase(gly_N > 0 & gly_O > 0, GL_LEV[5], gly_N > 0, GL_LEV[3], gly_O > 0, GL_LEV[4], gly_protein_level > 0, GL_LEV[2], default = GL_LEV[1])]
# attach to the drawn proteins and save
# (joined by gene symbol; metabolite nodes get empty annotation columns)
N <- A[, !"entrez_gene"][N, on = c(gene_symbol = "node")]; setnames(N, "gene_symbol", "node")
# Safety check: every drawn protein has both categories.
stopifnot(!anyNA(N[node_type == "protein", ph_cat]), !anyNA(N[node_type == "protein", gl_cat]))
# Save the per-protein annotation table.
fwrite(N[node_type == "protein", .(protein = node, ph_measured, ph_resp_EE, ph_resp_RE, ph_resp_any, ph_cat, gly_sites, gly_N, gly_O,
                                   gly_protein_level, glycans, gl_cat)], file.path(OUT, "16_protein_annotation.csv"))

# ---- drawing (figure 15b + PTM tags on the proteins) -----------------------------------------------------------------
# Figure theme: small text, no axes, legend at the bottom, grey caption under the figure.
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
# The shared base figure (15b without the protein circles), reused by 16a, 16b and 16c.
base_plot <- function() ggplot() +
  # class bubbles
  geom_mark_hull(data = N[node_type == "metabolite"], aes(x, y, group = class), colour = "grey45", fill = "grey60", alpha = 0.08,
                 linewidth = 0.25, expand = unit(2.2, "mm"), radius = unit(2, "mm"), concavity = 3) +
  # edges: colour and width = difference, line type = edge type
  geom_segment(data = D, aes(x, y, xend = xend, yend = yend, colour = w_diff, linewidth = abs(w_diff), linetype = edge_type), lineend = "round") +
  # metabolites: grey triangles, as 15b
  geom_point(data = N[node_type == "metabolite"], aes(x, y, size = abs_delta, shape = node_type), fill = "grey80", colour = "grey25", stroke = 0.22) +
  # class names
  geom_text_repel(data = cn, aes(x, y, label = class), size = 2.4, colour = "grey30", fontface = "bold.italic", box.padding = 0.15,
                  min.segment.length = 0.3, segment.size = 0.15, segment.colour = "grey50", max.overlaps = Inf, seed = SEED, max.time = 60, max.iter = 1e4) +
  # blue - grey - red difference scale centred at zero
  scale_colour_gradient2(low = COL_RE, mid = COL_SAME, high = COL_EE, midpoint = 0, limits = c(-lim_d, lim_d), oob = scales::squish,
                         breaks = c(-lim_d, 0, lim_d), labels = c("higher in resistance\n(w_RE > w_EE)", "same", "higher in endurance\n(w_EE > w_RE)"),
                         name = "edge difference  w_EE − w_RE",
                         guide = guide_colourbar(order = 10, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5)) +
  # line types, widths, node sizes, shapes
  scale_linetype_manual(values = TYPE_LTY, name = "edge type") + scale_linewidth(range = c(0.1, 1.6), guide = "none") +
  scale_size(range = c(0.6, 4), name = "node strength difference (absolute)") + scale_shape_manual(values = SHAPES, name = "node") +
  # plotting area and theme
  coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.05), clip = "off") + theme_net()

# ---- 16c: both layers + site-level crosstalk -------------------------------------------------------------------
# Crosstalk residue = a serine / threonine / tyrosine that MoTrPAC measured as a phosphosite AND GlyGen lists as an
# O-glycosylation site (mostly O-GlcNAc) on the same canonical protein, position and residue. The two modifications
# compete for the same hydroxyl group, so these are candidate phospho / O-GlcNAc switch sites.
# O-glycosylation sites (Ser / Thr / Tyr) from mnet (UniProt) and GlyGen (other databases), as canonical site IDs
# The GlyGen site file (from network/inventory); stop with a message if it is missing.
gsf <- file.path(OUT, "inventory", "glygen_glycosites.csv")
if (!file.exists(gsf)) stop("re-run network/inventory/glygen_protein_inventory.py (", gsf, " missing)")
# mnet (UniProt) O-glycosylation sites on Ser / Thr / Tyr, as canonical site IDs (accession_residue+position)
og_m <- fread(file.path(MNET_DIR, "glycosites.csv"))[residue %in% c("S", "T", "Y") & glyco_type %in% c("O-linked", "O-GlcNAc"), .(site_id, gly_type = glyco_type, gly_source = "UniProt (mnet)")]
# GlyGen O-linked sites, rewritten into the same site-ID format (accession without isoform, one-letter residue, position)
og_g <- fread(gsf)[type == "O-linked" & residue %in% c("Ser", "Thr", "Tyr"),
                   .(site_id = paste0(sub("-.*$", "", glygen_ac), "_", substr(residue, 1, 1), position), gly_type = fifelse(subtype == "O-GlcNAcylation", "O-GlcNAc", "O-linked"), gly_source = paste0("GlyGen: ", source))]
# one row per site: all glycosylation subtypes and sources that report it
og <- unique(rbind(og_m, og_g))[, .(gly_subtypes = paste(sort(unique(gly_type)), collapse = ";"), gly_sources = paste(sort(unique(gly_source)), collapse = "; ")), by = site_id]
# MoTrPAC features -> canonical sites through mnet's isoform-aware bridge (residue mismatches excluded)
br <- unique(fread(file.path(MNET_DIR, "motrpac_feature_site_map.csv"))[mapping_status != "residue_mismatch" & site_id != "", .(feature_id, site_id, n_sites)])
# phospho features -> canonical sites -> keep only sites that are also O-glycosylation sites
xt <- merge(ph[, .(feature_id, entrez_gene, sig_EE, sig_RE)], br, by = "feature_id")[og, on = "site_id", nomatch = 0]
xt <- genes[xt, on = "entrez_gene", nomatch = 0]   # our 471 proteins only
# one row per residue: features measuring it, and whether any of them responds (EE, RE); single-site evidence flagged
XT <- xt[, .(features = paste(sort(unique(feature_id)), collapse = ";"), n_features = uniqueN(feature_id),
             single_site_feature = any(n_sites == 1), responds_EE = any(sig_EE), responds_RE = any(sig_RE),
             gly_subtypes = gly_subtypes[1], gly_sources = gly_sources[1]),
         by = .(protein = gene_symbol, residue = sub("^[^_]*_", "", site_id), position = as.integer(sub("^[^_]*_[STY]", "", site_id)))][order(protein, position)]
# summarise the arm pattern in one word
XT[, responds := fcase(responds_EE & responds_RE, "both", responds_EE, "endurance", responds_RE, "resistance", default = "no")]
# save (the helper position column is left out)
fwrite(XT[, !"position"], file.path(OUT, "16_crosstalk_sites.csv"))
# per protein: crosstalk residues and how many respond
xp <- XT[, .(xt_n = .N, xt_resp = sum(responds != "no"), xt_label = paste0(protein[1], "  ", paste(residue, collapse = ","))), by = protein]
# long lists (more than 4 residues) are shortened to a count
xp[, xt_label := fifelse(xt_n > 4, sprintf("%s  %d residues", protein, xt_n), xt_label)]
# add how many of the residues respond
xp[, xt_label := sprintf("%s (%d/%d respond)", xt_label, xt_resp, xt_n)]

# ---- PTM tags: a stalk from each protein ending in a symbol (as in the interactive pages) --------------------
# Tag vocabulary. MoTrPAC phosphosites (measured): one pin per responding site (up to 6, then "+n"), red = responds
# after endurance only, blue = after resistance only, purple = after both (adj. p < 0.05, any tissue / time).
# Glycosylation in SNFG symbols (database knowledge): blue square = N-linked (GlcNAc), yellow square = O-linked
# mucin-type (GalNAc), blue square with a white dot = O-GlcNAc, white square = glycosylated, site unknown.
# Crosstalk: a diamond, gold if a shared phospho / O-glycosylation residue responds, white otherwise.
# The tag types (in legend order), their fill colours and their point shapes (21 circle, 22 square, 23 diamond).
TAG_LEV <- c("phosphosite responds after endurance only", "phosphosite responds after resistance only", "phosphosite responds after both",
             "N-linked glycosylation (GlcNAc)", "O-linked glycosylation, mucin-type (GalNAc)", "O-GlcNAc", "glycosylated, site unknown",
             "phospho = O-glyco residue, responds", "phospho = O-glyco residue, none respond")
TAG_FILL <- setNames(c("#E41A1C", "#377EB8", "#984EA3", "#0072BC", "#FFD400", "#0072BC", "#FFFFFF", "#FFD700", "#FFFFFF"), TAG_LEV)
TAG_SHAPE <- setNames(c(21, 21, 21, 22, 22, 22, 22, 23, 23), TAG_LEV)
# responding MoTrPAC sites per drawn protein, with the arm pattern (sites responding in both arms first)
# (one row per responding phospho feature)
site_arm <- genes[ph[sig_EE | sig_RE, .(feature_id, entrez_gene, arm = fifelse(sig_EE & sig_RE, TAG_LEV[3], fifelse(sig_EE, TAG_LEV[1], TAG_LEV[2])))], on = "entrez_gene", nomatch = 0]
# sort key: sites responding in both arms first, then endurance only, then resistance only
site_arm[, ord := match(arm, TAG_LEV[c(3, 1, 2)])]
# build the tag list for a set of layers ("phospho", "glyco", "crosstalk")
build_tags <- function(layers) {
  # the drawn proteins, and an empty list that collects the tags of each requested layer
  P <- N[node_type == "protein"]
  out <- list()
  if ("phospho" %in% layers) {
    # responding sites of the drawn proteins, numbered per protein in sort order
    s <- site_arm[gene_symbol %in% P$node][order(gene_symbol, ord)][, k := seq_len(.N), by = gene_symbol]
    # at most 6 pins per protein
    out$p <- s[k <= 6, .(node = gene_symbol, tag = arm, n = NA_integer_, more = NA_integer_)]
    # proteins with more than 6 responding sites get one extra "+n" row (no symbol)
    mo <- s[, .(extra = .N - 6L), by = gene_symbol][extra > 0]
    if (nrow(mo)) out$m <- mo[, .(node = gene_symbol, tag = NA_character_, n = NA_integer_, more = extra)]
  }
  # glycosylation: one square per type present, with its site count; "site unknown" only when no site is known
  if ("glyco" %in% layers) out$g <- rbind(P[gly_N > 0, .(node, tag = TAG_LEV[4], n = gly_N)], P[gly_O > 0, .(node, tag = TAG_LEV[5], n = gly_O)],
                                          P[gly_OG > 0, .(node, tag = TAG_LEV[6], n = gly_OG)], P[gly_sites == 0 & gly_protein_level > 0, .(node, tag = TAG_LEV[7], n = NA_integer_)])[, more := NA_integer_]
  # crosstalk: one diamond per protein, gold if any shared residue responds, with the residue count
  if ("crosstalk" %in% layers) out$x <- xp[protein %in% P$node, .(node = protein, tag = fifelse(xt_resp > 0, TAG_LEV[8], TAG_LEV[9]), n = xt_n, more = NA_integer_)]
  # all tags in one table
  Tg <- rbindlist(out, fill = TRUE)
  # fan the tags clockwise from the upper right; y is stretched so the stalks look equally long on the page
  # number each protein's tags (1, 2, ...) so each gets its own angle
  Tg[, k := seq_len(.N), by = node]
  # add the protein's position
  Tg <- N[, .(node, x0 = x, y0 = y)][Tg, on = "node"]
  # first tag at 70 degrees, each next one 32 degrees further clockwise; L = stalk length; ASP = y stretch
  ang <- (70 - (Tg$k - 1) * 32) * pi / 180; L <- 0.021; ASP <- 1.95
  # hx, hy = where the tag symbol sits; cx, cy = a little further out, where its count is written
  Tg[, `:=`(hx = x0 + L * cos(ang), hy = y0 + L * ASP * sin(ang), cx = x0 + 1.55 * L * cos(ang), cy = y0 + 1.55 * L * ASP * sin(ang))]
  # return the tag table
  Tg
}
# one tagged figure: 15b with grey protein circles and the chosen tag layers
tagged <- function(layers, labels, title, caption, file) {
  # the tags for the chosen layers, and the protein nodes
  Tg <- build_tags(layers); P <- N[node_type == "protein"]
  # tags that have a symbol (the "+n" rows have none), in the fixed legend order
  hd <- Tg[!is.na(tag)]; hd[, tag := factor(tag, levels = TAG_LEV)]
  # legend text with the number of proteins carrying each tag
  cnt <- hd[, .(np = uniqueN(node)), by = tag]; shown <- TAG_LEV[TAG_LEV %in% as.character(hd$tag)]
  lab <- setNames(sprintf("%s (%d proteins)", shown, cnt$np[match(shown, cnt$tag)]), shown)
  # start from the shared base figure (bubbles, edges, metabolites, class names)
  q <- base_plot() +
    # stalks (under the nodes)
    geom_segment(data = Tg, aes(x0, y0, xend = hx, yend = hy), colour = "grey35", linewidth = 0.25) +
    # proteins: grey circles, as in 15b
    geom_point(data = P, aes(x, y, size = abs_delta, shape = node_type), fill = "grey80", colour = "grey25", stroke = 0.22) +
    # fresh shape and size scales for the tag symbols (the node scales are already taken)
    ggnewscale::new_scale("shape") + ggnewscale::new_scale("size") +
    # the tag symbols at the end of each stalk
    geom_point(data = hd, aes(hx, hy, fill = tag, shape = tag), size = 1.9, colour = "grey15", stroke = 0.3) +
    # fill and shape share one legend (identical name, breaks, labels and guide, so ggplot merges them)
    scale_fill_manual(values = TAG_FILL[shown], labels = lab, breaks = shown, name = "PTM tags",
                      guide = guide_legend(order = 1, ncol = 2, override.aes = list(size = 2.6))) +
    scale_shape_manual(values = TAG_SHAPE[shown], labels = lab, breaks = shown, name = "PTM tags",
                       guide = guide_legend(order = 1, ncol = 2, override.aes = list(size = 2.6))) +
    # the "P" on phosphosite pins, the white dot on O-GlcNAc, "+n" for more sites, and site counts
    geom_text(data = hd[grepl("^phosphosite", tag)], aes(hx, hy, label = "P"), size = 1.05, colour = "white", fontface = "bold") +
    geom_point(data = hd[tag == TAG_LEV[6]], aes(hx, hy), shape = 16, size = 0.45, colour = "white") +
    geom_text(data = Tg[!is.na(more)], aes(hx, hy, label = paste0("+", more)), size = 1.4, colour = "grey10", fontface = "bold") +
    geom_text(data = Tg[!is.na(n) & n > 1], aes(cx, cy, label = n), size = 1.25, colour = "grey20") +
    # the protein labels chosen for this figure
    geom_text_repel(data = labels, aes(x, y, label = lab), size = 2.05, colour = "grey10", min.segment.length = 0.1, segment.size = 0.15,
                    box.padding = 0.45, max.overlaps = Inf, seed = SEED, max.time = 60, max.iter = 1e4) +
    # legend order and look for the edge colour bar and edge types
    guides(colour = guide_colourbar(order = 2, barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5),
           linetype = guide_legend(order = 3, override.aes = list(colour = "grey30", linewidth = 0.6))) +
    # title and caption
    labs(title = title, caption = caption)
  # save and report
  ggsave(file.path(FIG, file), q, width = 11.5, height = 8.4, dpi = 300, bg = "white")
  message("-> ", file.path(FIG, file))
}
# labels: the proteins that matter for each figure
# The drawn proteins.
Pp <- N[node_type == "protein"]
# 16a labels: the proteins with the most responding sites ("responding / measured").
lab_ph <- Pp[ph_resp_any > 0][order(-ph_resp_any, -ph_measured)][seq_len(min(NLAB, .N)), .(x, y, lab = sprintf("%s  %d/%d", node, ph_resp_any, ph_measured))]
# 16b labels: the proteins with the most glycosylation sites.
lab_gl <- Pp[gly_sites > 0][order(-gly_sites, -glycans)][seq_len(min(NLAB, .N)), .(x, y, lab = sprintf("%s  %d sites", node, gly_sites))]
# 16c labels: every drawn protein with a crosstalk residue.
lab_xt <- xp[Pp, on = c(protein = "node"), nomatch = 0][, .(x, y, lab = xt_label)]
# Figure captions (phospho, glycosylation).
CAP_PH <- "MoTrPAC phosphosites (muscle 0.5 / 4 / 24 h, adipose 4 h): one pin per site responding at adj. p < 0.05 vs control at any time point (up to 6, then +n)."
CAP_GL <- "Glycosylation (database knowledge, not measured here): UniProt sites via the team's mnet resource, SNFG symbols, number = sites; site unknown = GlyGen protein-level evidence."
# Draw and save 16a (phospho tags only).
tagged("phospho", lab_ph, "Joint protein-metabolite network, endurance minus resistance edge weights, with MoTrPAC phosphosite tags",
       paste0(CAP_PH, " Labels: proteins with the most responding sites (responding / measured)."), "16a_joint_edge_difference_phospho.png")
# Draw and save 16b (glycosylation tags only).
tagged("glyco", lab_gl, "Joint protein-metabolite network, endurance minus resistance edge weights, with glycosylation tags",
       paste0(CAP_GL, " Labels: proteins with the most glycosylation sites."), "16b_joint_edge_difference_glycosylation.png")
# Draw and save 16c (phospho, glycosylation and crosstalk tags).
tagged(c("phospho", "glyco", "crosstalk"), lab_xt, "Joint protein-metabolite network, endurance minus resistance edge weights, with phosphosite, glycosylation and crosstalk tags",
       paste0(CAP_PH, "\n", CAP_GL, " Diamond: a measured phosphosite is also a known O-glycosylation site (labels list those residues)."),
       "16c_joint_edge_difference_phospho_glyco_crosstalk.png")
# Print a one-line summary of the crosstalk residues.
cat(sprintf("crosstalk residues: %d on %d proteins (%d drawn in the network); responding: %d\n", nrow(XT), uniqueN(XT$protein), nrow(lab_xt), sum(XT$responds != "no")))
