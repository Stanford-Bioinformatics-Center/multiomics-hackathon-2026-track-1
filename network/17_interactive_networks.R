#!/usr/bin/env Rscript
# =====================================================================================================
# 17_interactive_networks.R — STEP 17: INTERACTIVE (NAVIGABLE) VERSIONS OF THE NETWORKS + CYTOSCAPE FILES
# =====================================================================================================
#
# PURPOSE (the question this answers)
#   The static figures (10, 11, 14, 15) show whole networks at once; individual genes, metabolites and
#   edges are hard to find and read. This step makes the three networks navigable, the way Cytoscape's
#   views are: zoom and pan, search, highlight a node's neighbours, filter by edge type or metabolite
#   class, switch between endurance / resistance / difference views, collapse classes, hover for values.
#   It also writes Cytoscape import files so teammates who use Cytoscape get its full toolset on our data.
#
# WHAT THIS SCRIPT DOES (plain language)
#   Nothing is recomputed: nodes, edges, weights and layouts are read from earlier steps.
#     joint network      : steps 14 (edges, nodes) and 15 (class-grouped layout)
#     gene network       : steps 1, 3 (edges) and 10 (layout, as in figures 10a / 11a)
#     metabolite network : steps 1b, 6 (edges) and 10 (layout, as in figures 10b / 11b)
#   For each network, one self-contained HTML page (opens in any browser, no install) with:
#     - FILTERS: omes (RNA / protein / metabolites), tissues, time points, arm (endurance vs control,
#       resistance vs control, endurance minus resistance) and an adj. p threshold. Every filter change
#       recomputes, in the browser: node colour (mean normalised response over the selected dimensions),
#       a black outline where the node is significant in a selected cell (MoTrPAC adj. p), and EDGE WEIGHTS
#       as the same dot products as steps 3 / 6 / 14 restricted to the selected dimensions (with all filters
#       on they equal the pipeline weights; checked to < 1e-6);
#     - MODULES (17_filter_stats.R): structural communities with MoTrPAC run_cameraPR results; the menu marks
#       modules significant in the selection; choosing one highlights it and tabulates its tests;
#     - ANNOTATION LAYERS ("colour nodes by"): MoTrPAC phosphosites responding in the selection (up / down /
#       both), every GlyGen field (phosphosites, kinase sites, glycosylation, glycans, crosstalk residues,
#       mutations, disease, biomarkers, PTM / site / enzyme annotations, pathways, reactions, expression,
#       publications), and GlyGen KINASE -> SUBSTRATE arrows between network proteins (red when a substrate
#       site responds in the selection); tooltips give per-cell logFC / adj. p, responding sites (with kinase
#       and O-GlcNAc flags) and the GlyGen summary;
#     - search, click-to-highlight neighbours, metabolite-class selector, edge-type check boxes, class outlines
#       ("bubbles", as in 15a / 15b), collapse / expand classes. Proteins are circles, metabolites triangles.
#   And for each network a Cytoscape.js JSON file (.cyjs, with node positions) plus one Cytoscape style
#   file with three styles (EE, RE, difference), and plain node / edge tables.
#
# HOW TO RUN
#   After steps 10, 14, 15, network/inventory and 17_filter_stats.R:   Rscript network/17_interactive_networks.R   (about 30 seconds; needs pandoc,
#   which ships with RStudio / Positron / Quarto, to make the pages self-contained)
#   Open:  $HACK_FIG/17_interactive/17a_joint_network.html (and 17b, 17c) in a browser.
#   Cytoscape desktop: File > Import > Network from File > 17_*_network.cyjs; File > Import > Styles from
#   File > 17_cytoscape_styles.xml; then pick "hackathon EE", "hackathon RE" or "hackathon difference" in
#   the Style panel. Positions come with the .cyjs file (no layout needed).
#
# DATA AND PROVENANCE
#   MoTrPAC human pre-suspension results (MotrpacHumanPreSuspensionAnalysis v0.2.4), normalised per ome
#   (steps 1 / 1b); STRING combined score >= 700 (step 2); Rhea release 142 (step 5); RefMet classes
#   (step 1c). Upstream QC is documented in those steps and in the README ("Critical QC step").
#
# TECH STACK
#   R 4.4; data.table, visNetwork 2.1 (vis-network JavaScript library), htmlwidgets (+ pandoc for
#   self-contained pages), htmltools, jsonlite. Plain JavaScript (embedded below) adds the controls.
#
# INPUTS (files and columns)
#   $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 15_class_layout.csv       joint network
#   $HACK_OUT/01_nodes_{EE,RE}.csv, 03_weighted_edges.csv, 10_layout_genes.csv   gene network
#   $HACK_OUT/01b_metab_nodes_{EE,RE}.csv, 06_metabolite_edges.csv, 10_layout_metabolites.csv  metabolites
#   $HACK_OUT/05_metabolite_protein_links.csv (Rhea reactions), 01c_metabolite_ids.csv (classes)
#   $HACK_OUT/17_node_cell_stats.csv, 17_modules.csv, 17_module_camera.csv, 17_phospho_site_stats.csv,
#     17_kinase_edges.csv, 17_glygen_protein_annotation.csv   (from 17_filter_stats.R)
#
# OUTPUTS (files, locations)
#   $HACK_FIG/17_interactive/17a_joint_network.html, 17b_gene_network.html, 17c_metabolite_network.html
#   $HACK_OUT/17_cytoscape/17_{joint,gene,metabolite}_network.cyjs    Cytoscape.js JSON with positions
#   $HACK_OUT/17_cytoscape/17_{joint,gene,metabolite}_{nodes,edges}.csv  the same data as plain tables
#   $HACK_OUT/17_cytoscape/17_cytoscape_styles.xml                     three Cytoscape styles
#   (never in the repo)
#
# EXPECTED OUTPUT (2026-09-26) AND VALIDATION
#   joint 364 nodes / 764 edges; genes 286 / 431; metabolites 44 / 147. The script stops if a network's
#   node or edge count differs from its source table, if a layout is missing a node, or if an edge weight
#   differs from the source. 99_validate_outputs.R re-checks the .cyjs files against the source tables.
#
# KNOWN LIMITS
#   The Cytoscape files were written to the documented formats but not opened in Cytoscape desktop here
#   (it is not installed on the machine that produced them). Colours in the pages are relative to each
#   network's own 95th percentile (as in the static figures); differences between arms are not tested.
# =====================================================================================================

# Load packages quietly.
suppressMessages({ library(data.table); library(visNetwork); library(htmlwidgets); library(htmltools); library(jsonlite) })

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv("HACK_OUT", unset = path.expand("~/Desktop/output/hackathon-2026-track1/network"))
# Where figures go (override with HACK_FIG); outside the repo on purpose.
FIG <- Sys.getenv("HACK_FIG", unset = path.expand("~/Desktop/output/hackathon"))
# The interactive pages and the Cytoscape files get their own folders.
HTML_DIR <- file.path(FIG, "17_interactive"); CY_DIR <- file.path(OUT, "17_cytoscape")
dir.create(HTML_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(CY_DIR, recursive = TRUE, showWarnings = FALSE)
# Canvas size in pixels that the 0..1 layouts are stretched to.
W <- 1600; H <- 1000
# Edge-type colours (as in step 14) and the colours of the two diverging scales (as in steps 10 / 11 / 14).
TYPE_COL <- c("protein - protein" = "#8C8C8C", "metabolite - metabolite" = "#1B7837", "metabolite - protein" = "#8C510A")
RESP_PAL <- c("#6A3D9A", "#FFFFFF", "#E66100")   # mean response: down - none - up
DIFF_PAL <- c("#2166AC", "#D9D9D9", "#B2182B")   # difference: higher in RE - same - higher in EE

# ---- helpers ------------------------------------------------------------------------------------------
# Map values onto a three-colour scale symmetric around 0, squashing anything beyond +/- lim.
ramp3 <- function(v, lim, pal) { t <- (pmin(pmax(v / lim, -1), 1) + 1) / 2; rgb(colorRamp(pal)(t), maxColorValue = 255) }
# Short number formatting for tooltips.
f3 <- function(v) formatC(v, digits = 3, format = "g")
# HTML-escape text for tooltips.
esc <- function(v) htmlEscape(as.character(v))
# Node size in pixels from a non-negative quantity (square-root scaling, relative to the network's largest).
nsize <- function(v) 6 + 16 * sqrt(v / max(v, na.rm = TRUE))
# Edge width in pixels from a non-negative quantity, capped at its 95th percentile.
ewidth <- function(v, cap) 0.6 + 5 * pmin(v / cap, 1)

# Build one network's node and edge tables with every view's colours, sizes and tooltips.
#   N: node, node_type (protein / metabolite), class, resp_EE, resp_RE, x, y (0..1), hull (outline group)
#   E: a, b, edge_type, w_EE, w_RE, w_diff, info (tooltip HTML with the edge's evidence)
prepare <- function(N, E) {
  # node strength per arm = sum of |w| over its edges; difference = sum|w_EE| - sum|w_RE|
  both <- rbind(E[, .(node = a, w_EE, w_RE)], E[, .(node = b, w_EE, w_RE)])
  st <- both[, .(strength_EE = sum(abs(w_EE)), strength_RE = sum(abs(w_RE)), degree = .N), by = node]
  N <- st[N, on = "node"][, strength_diff := strength_EE - strength_RE]
  # colour limits: 95th percentile of |mean response| (both arms) and of |w_diff| (as in the figures)
  lim_r <- as.numeric(quantile(abs(c(N$resp_EE, N$resp_RE)), 0.95)); lim_d <- as.numeric(quantile(abs(E$w_diff), 0.95))
  # node colours per view (difference view: plain grey, as in figures 11 / 14b)
  N[, `:=`(colEE = ramp3(resp_EE, lim_r, RESP_PAL), colRE = ramp3(resp_RE, lim_r, RESP_PAL), colD = "#CCCCCC")]
  # node sizes per view (EE / RE on one shared scale so the arms are comparable; difference = |strength difference|)
  smax <- max(N$strength_EE, N$strength_RE)
  N[, `:=`(sizeEE = 6 + 16 * sqrt(strength_EE / smax), sizeRE = 6 + 16 * sqrt(strength_RE / smax), sizeD = nsize(abs(strength_diff)))]
  # node tooltips
  N[, title := sprintf("<b>%s</b><br>%s%s<br>mean response: EE %s · RE %s<br>strength (sum |w|): EE %s · RE %s<br>edges: %d",
                       esc(node), node_type, fifelse(node_type == "metabolite", paste0(" · ", esc(class)), ""),
                       f3(resp_EE), f3(resp_RE), f3(strength_EE), f3(strength_RE), degree)]
  # edge colours, widths and line styles per view
  cap <- as.numeric(quantile(abs(c(E$w_EE, E$w_RE)), 0.95))
  E[, `:=`(colEE = TYPE_COL[as.character(edge_type)], colRE = TYPE_COL[as.character(edge_type)], colD = ramp3(w_diff, lim_d, DIFF_PAL),
           wdEE = ewidth(abs(w_EE), cap), wdRE = ewidth(abs(w_RE), cap), wdD = 0.6 + 6 * pmin(abs(w_diff) / lim_d, 1),
           # EE / RE: dashed when the weight is negative; difference view: line style = edge type
           dashEE = fifelse(w_EE < 0, "neg", "solid"), dashRE = fifelse(w_RE < 0, "neg", "solid"),
           dashD = c("protein - protein" = "solid", "metabolite - metabolite" = "mm", "metabolite - protein" = "mp")[as.character(edge_type)])]
  # edge tooltips
  E[, title := sprintf("<b>%s — %s</b><br>%s<br>w_EE %s · w_RE %s · w_EE − w_RE %s%s",
                       esc(a), esc(b), edge_type, f3(w_EE), f3(w_RE), f3(w_diff), fifelse(is.na(info) | info == "", "", paste0("<br>", info)))]
  list(N = N, E = E, lim_r = lim_r, lim_d = lim_d)
}

# ---- page data: everything the filters, modules and annotation layers need (from 17_filter_stats.R) --------
# Normalised node vectors (steps 1 / 1b) keyed "tissue|ome|time", per arm; used to recompute edge weights over
# the selected dimensions only (the same dot products as steps 3, 6 and 14, restricted to the selection).
vec_tab <- function(fE, fR, id) {
  e <- fread(file.path(OUT, fE)); r <- fread(file.path(OUT, fR)); d <- setdiff(names(e), c("entrez_gene", "gene_symbol", "metabolite"))
  # the always-empty adipose protein columns are read as logical; make every value column numeric
  for (x in list(e, r)) for (k in d) set(x, j = k, value = as.numeric(x[[k]]))
  m <- function(x, arm) { l <- melt(x[, c(id, d), with = FALSE], id.vars = id, variable.name = "dim", value.name = "v", na.rm = TRUE)
    l[, `:=`(key = gsub("_", "|", as.character(dim)), arm = arm)]; setnames(l, id, "node"); l[, .(node, arm, key, v)] }
  rbind(m(e, "EE"), m(r, "RE"))
}
VEC <- rbind(vec_tab("01_nodes_EE.csv", "01_nodes_RE.csv", "gene_symbol"), vec_tab("01b_metab_nodes_EE.csv", "01b_metab_nodes_RE.csv", "metabolite"))
NSTAT <- fread(file.path(OUT, "17_node_cell_stats.csv"))
PHS <- fread(file.path(OUT, "17_phospho_site_stats.csv"))
MODS <- fread(file.path(OUT, "17_modules.csv")); CAM <- fread(file.path(OUT, "17_module_camera.csv"))
MNAME <- fread(file.path(OUT, "17_module_names.csv")); MORA <- fread(file.path(OUT, "17_module_ora.csv"))
KIN <- fread(file.path(OUT, "17_kinase_edges.csv")); ANN <- fread(file.path(OUT, "17_glygen_protein_annotation.csv"))
for (f in c("17_node_cell_stats.csv", "17_phospho_site_stats.csv", "17_modules.csv", "17_module_camera.csv", "17_kinase_edges.csv"))
  if (!file.exists(file.path(OUT, f))) stop("run network/17_filter_stats.R first (", f, " missing)")
# Annotation fields offered in the "colour nodes by" menu: column, label, kind (count / binary).
# Node-colour menu: non-PTM annotations only (phosphorylation and glycosylation are drawn as tags on the nodes).
ANN_FIELDS <- list(list(col = "is_kinase", label = "mnet: is a kinase", binary = TRUE), list(col = "substrate_sites", label = "mnet: substrate sites it phosphorylates (as kinase)"),
                   list(col = "mutations", label = "GlyGen: mutations / SNVs"), list(col = "disease", label = "GlyGen: disease associations"),
                   list(col = "biomarkers", label = "GlyGen: biomarkers"), list(col = "ptm_annotation", label = "GlyGen: PTM annotations"),
                   list(col = "site_annotation", label = "GlyGen: active / binding sites"), list(col = "enzyme", label = "GlyGen: enzyme (EC) annotations"),
                   list(col = "pathways", label = "GlyGen: pathways"), list(col = "reactions", label = "GlyGen: reactions"),
                   list(col = "expression_tissues", label = "GlyGen: normal tissues expressed"), list(col = "publications", label = "GlyGen: publications"))

# Assemble the page data for one network as a JSON string (parsed by the page's JavaScript).
page_data <- function(P, net, types) {
  nodes <- P$N$node
  # vectors: node -> {EE: {key: v}, RE: {key: v}}
  v <- VEC[node %in% nodes]
  vecs <- lapply(split(v, v$node), function(x) lapply(split(x, x$arm), function(y) as.list(setNames(round(y$v, 6), y$key))))
  # node statistics: node -> {"tissue|ome|arm|time": [logFC, adj_p]}
  s <- NSTAT[node %in% nodes]
  stats <- lapply(split(s, s$node), function(x) setNames(lapply(seq_len(nrow(x)), function(i) c(signif(x$logFC[i], 4), signif(x$adj_p[i], 3))),
                                                          paste(x$tissue, x$ome, x$arm, x$time, sep = "|")))
  # phosphosites: protein -> [{site, f, kin, xt, known, c: {"tissue|arm|time": [logFC, adj_p]}}]
  ph <- PHS[protein %in% nodes]
  phos <- lapply(split(ph, ph$protein), function(x) unname(lapply(split(x, x$feature_id), function(y)
    list(site = y$site[1], f = y$feature_id[1], kin = ifelse(is.na(y$kinases[1]), "", y$kinases[1]), xt = isTRUE(y$crosstalk[1]), known = isTRUE(y$known_in_glygen[1]),
         c = setNames(lapply(seq_len(nrow(y)), function(i) c(signif(y$logFC[i], 4), signif(y$adj_p[i], 3))), paste(y$tissue, y$arm, y$time, sep = "|"))))))
  # modules of this network with their CAMERA-PR results: [{id, members, cam: {"tissue|ome|arm|time": [z, dir, fdr, n]}}]
  md <- MODS[network == net & node %in% nodes]; cm <- CAM[network == net]
  mods <- unname(lapply(split(md, md$module), function(x) { cc <- cm[module == x$module[1]]
    top <- MORA[module == x$module[1] & adj_p_value < 0.05 & database != "CELLMARKER"][order(adj_p_value)][seq_len(min(.N, 5))]
    list(id = x$module[1], name = MNAME[module == x$module[1], name], members = I(x$node), n_prot = sum(x$node_type == "protein"), n_met = sum(x$node_type == "metabolite"),
         ora = unname(lapply(seq_len(nrow(top)), function(i) list(top$set_label[i], signif(top$adj_p_value[i], 2), top$overlap_n[i], gsub(";", ", ", top$overlap[i])))),
         cam = setNames(lapply(seq_len(nrow(cc)), function(i) list(signif(cc$z[i], 3), cc$direction[i], signif(cc$fdr[i], 3), cc$n_members_tested[i])),
                        paste(cc$tissue, cc$ome, cc$arm, cc$time, sep = "|"))) }))
  # GlyGen annotations: protein -> {field: count}
  an <- ANN[protein %in% nodes]
  ann <- lapply(split(an, an$protein), function(x) as.list(x[, !"protein"]))
  # kinase -> substrate edges inside this network (both ends drawn)
  ke <- KIN[kinase %in% nodes & substrate %in% nodes]
  kin <- unname(lapply(seq_len(nrow(ke)), function(i) as.list(ke[i])))
  omes <- intersect(c("rna", "prot", "metab"), unique(sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", v$key)))
  jsonlite::toJSON(list(net = net, vecs = vecs, stats = stats, phos = phos, mods = mods, ann = ann, kin = kin, omes = I(omes),
                        annFields = ANN_FIELDS, types = I(types)), auto_unbox = TRUE, digits = NA, na = "null")
}


# The controls (plain JavaScript, run once the page has drawn the network). cfg.json comes from page_data().
JS <- r"---(
function(el, x, cfg) {
  var D = JSON.parse(cfg.json), S = cfg.static;
  var net = document.getElementById("graph" + el.id).chart;
  var nodes = net.body.data.nodes, edges = net.body.data.edges;
  var TIS = ["adipose", "blood", "muscle"], TIMES = ["0.5h", "4h", "24h"], OMES = {rna: "RNA", prot: "protein", metab: "metabolites"};
  var st = { tags: {mp: 1, kp: 0, N: 1, O: 1, OG: 1, unk: 0, xt: 1}, pins: {}, omes: {}, tis: {adipose: 1, blood: 1, muscle: 1}, times: {"0.5h": 1, "4h": 1, "24h": 1}, arm: "EE", thr: 0.05, colour: "response",
             focus: null, collapsed: false, hulls: S.hulls.length > 0, typeOn: {}, kinase: false, module: "" };
  D.omes.forEach(function (o) { st.omes[o] = 1; }); D.types.forEach(function (t) { st.typeOn[t] = true; });
  var DASH = { solid: false, neg: [5, 5], mm: [10, 5], mp: [2, 4], kin: [3, 3] };
  var TYPE_COL = {"protein - protein": "#8C8C8C", "metabolite - metabolite": "#1B7837", "metabolite - protein": "#8C510A"};
  var RESP = ["#6A3D9A", "#FFFFFF", "#E66100"], DIFF = ["#2166AC", "#F2F2F2", "#B2182B"];
  // discrete bins for GlyGen counts: [from, to, label, colour] (ColorBrewer Greens)
  var BINS = [[0, 0, "0", "#FFFFFF"], [1, 1, "1", "#E5F5E0"], [2, 4, "2–4", "#C7E9C0"], [5, 9, "5–9", "#A1D99B"], [10, 24, "10–24", "#74C476"],
              [25, 99, "25–99", "#31A354"], [100, 1e12, "100 or more", "#006D2C"]];
  var ARMLAB = {EE: "endurance vs control", RE: "resistance vs control", ER: "endurance minus resistance"};
  // ---------- small helpers ----------
  function hex(c) { return [parseInt(c.substr(1, 2), 16), parseInt(c.substr(3, 2), 16), parseInt(c.substr(5, 2), 16)]; }
  function mix(a, b, t) { a = hex(a); b = hex(b); return "rgb(" + [0, 1, 2].map(function (i) { return Math.round(a[i] + (b[i] - a[i]) * t); }).join(",") + ")"; }
  function div3(v, lim, pal) { var t = Math.max(-1, Math.min(1, v / lim)); return t < 0 ? mix(pal[1], pal[0], -t) : mix(pal[1], pal[2], t); }
  function seq3(t, pal) { t = Math.max(0, Math.min(1, t)); return t < 0.5 ? mix(pal[0], pal[1], t * 2) : mix(pal[1], pal[2], (t - 0.5) * 2); }
  function q95(a) { a = a.filter(function (v) { return isFinite(v); }).map(Math.abs).sort(function (x, y) { return x - y; }); return a.length ? (a[Math.floor(0.95 * (a.length - 1))] || a[a.length - 1] || 1) : 1; }
  function f3(v) { return (v === null || v === undefined) ? "NA" : (Math.abs(v) >= 0.001 && Math.abs(v) < 1000 ? (+v).toPrecision(3) : (+v).toExponential(1)); }
  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;"); }
  function on(obj) { return Object.keys(obj).filter(function (k) { return obj[k]; }); }
  function grad(pal) { return "<span class='hk-grad' style='background:linear-gradient(90deg," + pal.join(",") + ")'></span>"; }
  // selected dimension keys "tissue|ome|time"
  function keys(omes) { var k = []; omes.forEach(function (o) { on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) { k.push(t + "|" + o + "|" + h); }); }); }); return k; }
  function vec(n, arm) { return (D.vecs[n] || {})[arm] || {}; }
  // edge weight over the selected dimensions for one arm; returns [w, number of contributing terms]
  function weight(e, arm) {
    var a = vec(e.from, arm), b = vec(e.to, arm), w = 0, n = 0;
    if (e.etype === "metabolite - protein") {
      if (!st.omes.metab) return [0, 0];
      var go = ["rna", "prot"].filter(function (o) { return st.omes[o]; });
      on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) {
        var m = a[t + "|metab|" + h]; if (m === undefined) return;
        go.forEach(function (o) { var g = b[t + "|" + o + "|" + h]; if (g !== undefined) { w += m * g; n++; } }); }); });
      return [w, n];
    }
    var om = (e.etype === "metabolite - metabolite") ? ["metab"] : ["rna", "prot"];
    keys(om.filter(function (o) { return st.omes[o]; })).forEach(function (k) { if (a[k] !== undefined && b[k] !== undefined) { w += a[k] * b[k]; n++; } });
    return [w, n];
  }
  // node mean response over the selected dimensions (arm EE / RE, or EE - RE)
  function response(n) {
    var ks = keys(on(st.omes)), vE = vec(n, "EE"), vR = vec(n, "RE"), s = 0, c = 0;
    ks.forEach(function (k) { if (vE[k] === undefined) return; s += (st.arm === "EE" ? vE[k] : st.arm === "RE" ? vR[k] : vE[k] - vR[k]); c++; });
    return c ? s / c : null;
  }
  // significant cells of a node for the current arm and selection
  function sigCells(n) {
    var out = [], S2 = D.stats[n] || {};
    on(st.omes).forEach(function (o) { on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) {
      var r = S2[t + "|" + o + "|" + st.arm + "|" + h]; if (r) out.push([t, o, h, r[0], r[1]]); }); }); });
    return out;
  }
  // phosphosites responding under the current filter (tissues muscle / adipose, arm, times)
  function phosCells(n) {
    var out = [];
    (D.phos[n] || []).forEach(function (s) { var best = null;
      on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) { var r = s.c[t + "|" + st.arm + "|" + h];
        if (r && (best === null || r[1] < best[4])) best = [t, h, r[0], s, r[1]]; }); });
      if (best) out.push(best); });
    return out;
  }
  // ---------- toolbar ----------
  var bar = document.createElement("div"); bar.className = "hk-bar";
  var cb = function (cls, val, lab, chk) { return "<label><input type='checkbox' class='" + cls + "' value='" + val + "'" + (chk ? " checked" : "") + "> " + lab + "</label>"; };
  var h = "<div><b>Omes</b> " + D.omes.map(function (o) { return cb("hk-ome", o, OMES[o], 1); }).join(" ") +
          "<span class='hk-sep'></span><b>Tissues</b> " + TIS.map(function (t) { return cb("hk-tis", t, t, 1); }).join(" ") +
          "<span class='hk-sep'></span><b>Time</b> " + TIMES.map(function (t) { return cb("hk-time", t, t, 1); }).join(" ") +
          "<span class='hk-sep'></span><b>Arm</b> " + ["EE", "RE", "ER"].map(function (a) { return "<label><input type='radio' name='" + el.id + "-arm' class='hk-arm' value='" + a + "'" + (a === "EE" ? " checked" : "") + "> " + ARMLAB[a] + "</label>"; }).join(" ") +
          "<span class='hk-sep'></span><b>adj. p &lt;</b> <input class='hk-thr' type='number' step='0.01' min='0' max='1' value='0.05' style='width:55px'></div>";
  h += "<div><b>Colour nodes by</b> <select class='hk-colour'><option value='response'>exercise response (filtered)</option>" +
       D.annFields.map(function (f) { return "<option value='ann:" + f.col + "'>" + esc(f.label) + "</option>"; }).join("") + "</select>" +
       "<span class='hk-sep'></span><b>Module</b> <select class='hk-mod'><option value=''>none</option></select>" +
       "<span class='hk-sep'></span><b>Find</b> <input class='hk-find' list='" + el.id + "-dl' placeholder='gene or metabolite'>" +
       "<datalist id='" + el.id + "-dl'>" + nodes.getIds().sort().map(function (i) { return "<option value=\"" + String(i).replace(/"/g, "&quot;") + "\">"; }).join("") + "</datalist>";
  if (S.classes.length) h += "<span class='hk-sep'></span><b>Class</b> <select class='hk-cls'><option value=''>all</option>" + S.classes.map(function (c) { return "<option>" + esc(c) + "</option>"; }).join("") + "</select>";
  h += "</div><div><b>PTM tags</b> (any combination) " + cb("hk-tag", "mp", "MoTrPAC phospho (red EE / blue RE / purple both)", 1) + " " +
       cb("hk-tag", "kp", "known phosphosites", 0) + " " + cb("hk-tag", "N", "N-linked glyco", 1) + " " + cb("hk-tag", "O", "O-linked glyco (GalNAc)", 1) + " " +
       cb("hk-tag", "OG", "O-GlcNAc", 1) + " " + cb("hk-tag", "unk", "glycosylated, site unknown", 0) + " " + cb("hk-tag", "xt", "phospho = O-glyco residue", 1);
  h += "</div><div><b>Edges</b> " + D.types.map(function (t) { return cb("hk-type", t, t, 1); }).join(" ") +
       (D.kin.length ? " " + cb("hk-kin", "1", "kinase → substrate (OmniPath via mnet; " + D.kin.length + ")", 0) : "") +
       (S.hulls.length ? "<span class='hk-sep'></span>" + cb("hk-hull", "1", "class outlines", 1) + " <button class='hk-col'>Collapse classes</button>" : "") +
       "<span class='hk-sep'></span><button class='hk-reset'>Reset</button></div>";
  bar.innerHTML = h;
  var legend = document.createElement("div"); legend.className = "hk-legend";
  var legendBox = document.createElement("div"); legendBox.className = "hk-legbox";
  var panel = document.createElement("div"); panel.className = "hk-panel";
  var help = document.createElement("div"); help.className = "hk-help";
  help.innerHTML = "Filters recompute node colours, significance outlines and edge weights (dot products over the selected dimensions only) · click a node to highlight it and its neighbours, empty space to clear · hover for values" + (S.classes.length ? " · double-click a collapsed class to open it" : "");
  el.parentNode.insertBefore(bar, el); el.parentNode.appendChild(legend); el.parentNode.appendChild(panel); el.parentNode.appendChild(help);
  var q = function (s) { return bar.querySelector(s); };
  // the network and its legend side by side: the legend panel sits to the right of the canvas, never over it
  var wrap = document.createElement("div"); wrap.className = "hk-wrap"; el.parentNode.insertBefore(wrap, el); wrap.appendChild(el); wrap.appendChild(legendBox);
  el.style.flex = "1 1 auto"; el.style.minWidth = "0"; el.style.width = "auto";
  wrap.parentNode.insertBefore(legend, wrap);   // the selection summary goes above the network, clear of the zoom buttons
  setTimeout(function () { net.setSize("100%", "760px"); net.redraw(); net.fit(); }, 60);
  // kinase edges are added once (hidden until switched on)
  D.kin.forEach(function (k, i) { edges.add({ id: "kin" + i, from: k.kinase, to: k.substrate, etype: "kinase", arrows: { to: { enabled: true, scaleFactor: 0.6 } },
    hidden: true, width: 1.6, dashes: DASH.kin, color: { color: "#111111", highlight: "#111111" }, smooth: { enabled: true, type: "curvedCW", roundness: 0.2 }, kin: k }); });
  // ---------- classes: collapse / outlines ----------
  function openAll() { S.classes.forEach(function (c) { var id = "class: " + c; if (net.isCluster(id)) net.openCluster(id); }); }
  function clusterAll() {
    S.classes.forEach(function (c) { var ids = nodes.getIds({ filter: function (n) { return n.group === c; } }); if (ids.length < 2) return;
      var p = net.getPositions(ids), cx = 0, cy = 0; ids.forEach(function (i) { cx += p[i].x; cy += p[i].y; });
      net.cluster({ joinCondition: function (o) { return o.group === c; }, clusterNodeProperties: { id: "class: " + c, label: c + " (" + ids.length + ")", shape: "triangle", size: 20,
        x: cx / ids.length, y: cy / ids.length, physics: false, font: { size: 18 }, color: { background: "#D9D9D9", border: "#404040" }, title: "<b>" + esc(c) + "</b><br>" + ids.length + " metabolites (double-click to open)" } }); });
  }
  function hull(pts) { pts.sort(function (a, b) { return a.x - b.x || a.y - b.y; }); if (pts.length < 3) return pts;
    var cr = function (o, a, b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }, lo = [], up = [];
    pts.forEach(function (p) { while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], p) <= 0) lo.pop(); lo.push(p); });
    pts.slice().reverse().forEach(function (p) { while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], p) <= 0) up.pop(); up.push(p); });
    lo.pop(); up.pop(); return lo.concat(up); }
  // PTM tags: a short stalk from the node edge ending in a symbol, fanned out clockwise from the upper right
  function star(ctx, x, y, r) { ctx.beginPath(); for (var i = 0; i < 10; i++) { var a = -Math.PI / 2 + i * Math.PI / 5, rr = i % 2 ? r * 0.45 : r; ctx.lineTo(x + rr * Math.cos(a), y + rr * Math.sin(a)); } ctx.closePath(); }
  net.on("afterDrawing", function (ctx) {
    var ids = Object.keys(st.pins); if (!ids.length) return;
    var pos = net.getPositions(ids);
    ids.forEach(function (id) { if (!pos[id] || net.isCluster(id)) return; var nd = nodes.get(id); if (!nd || nd.hidden) return;
      var r = (nd.size || 10), P = st.pins[id];
      ctx.save(); ctx.globalAlpha = (st.keep && !st.keep[id]) ? 0.12 : 1;
      P.forEach(function (p, k) { var ang = (-70 + k * 32) * Math.PI / 180, c = Math.cos(ang), sn = Math.sin(ang);
        var x0 = pos[id].x + c * r, y0 = pos[id].y + sn * r, x1 = pos[id].x + c * (r + 13), y1 = pos[id].y + sn * (r + 13), hx = pos[id].x + c * (r + 18.5), hy = pos[id].y + sn * (r + 18.5);
        ctx.strokeStyle = "#555"; ctx.lineWidth = 1.1; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
        ctx.lineWidth = 0.9; ctx.strokeStyle = "#222";
        if (p.kind === "P" || p.kind === "Pdb") { ctx.fillStyle = p.kind === "P" ? p.fill : "#D0D0D0"; ctx.beginPath(); ctx.arc(hx, hy, 5.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
          ctx.fillStyle = p.kind === "P" ? "#FFFFFF" : "#222"; ctx.font = "bold 7.5px Helvetica"; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText("P", hx, hy + 0.4); }
        else if (p.kind === "sq" || p.kind === "sqo" || p.kind === "sqh") { ctx.fillStyle = p.kind === "sqh" ? "#FFFFFF" : p.fill; ctx.strokeStyle = p.kind === "sqh" ? "#0072BC" : "#222";
          if (p.kind === "sqh") ctx.setLineDash([2, 1.5]); ctx.fillRect(hx - 4.5, hy - 4.5, 9, 9); ctx.strokeRect(hx - 4.5, hy - 4.5, 9, 9); ctx.setLineDash([]);
          if (p.kind === "sqo") { ctx.fillStyle = "#FFFFFF"; ctx.beginPath(); ctx.arc(hx, hy, 2, 0, 2 * Math.PI); ctx.fill(); } }
        else if (p.kind === "star") { ctx.fillStyle = "#FFD700"; star(ctx, hx, hy, 6.5); ctx.fill(); ctx.stroke(); }
        else if (p.kind === "more") { ctx.fillStyle = "#222"; ctx.font = "bold 8px Helvetica"; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText("+" + p.n, hx, hy); }
        if (p.n > 1 && p.kind !== "more") { ctx.fillStyle = "#222"; ctx.font = "7px Helvetica"; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText(p.n, pos[id].x + c * (r + 27), pos[id].y + sn * (r + 27)); }
      });
      ctx.restore(); });
  });
  net.on("beforeDrawing", function (ctx) { if (!st.hulls || st.collapsed) return;
    S.hulls.forEach(function (g) { var ids = nodes.getIds({ filter: function (n) { return n.hull === g.key; } }); if (!ids.length) return;
      var p = net.getPositions(ids), pts = ids.map(function (i) { return { x: p[i].x, y: p[i].y }; }), hp = hull(pts);
      var path = function () { ctx.beginPath(); ctx.moveTo(hp[0].x, hp[0].y); hp.forEach(function (v) { ctx.lineTo(v.x, v.y); }); ctx.closePath(); };
      ctx.save(); ctx.lineJoin = "round"; ctx.lineCap = "round";
      path(); ctx.lineWidth = 46; ctx.strokeStyle = "rgba(90,90,90,0.55)"; ctx.stroke(); path(); ctx.lineWidth = 43; ctx.strokeStyle = "#F3F3F3"; ctx.stroke(); ctx.fillStyle = "#F3F3F3"; ctx.fill();
      var top = Math.min.apply(null, pts.map(function (v) { return v.y; })), mx = pts.reduce(function (s, v) { return s + v.x; }, 0) / pts.length;
      ctx.font = "italic bold 15px Helvetica"; ctx.fillStyle = "#4D4D4D"; ctx.textAlign = "center"; ctx.fillText(g.label, mx, top - 28); ctx.restore(); }); });
  // ---------- the main redraw ----------
  function apply() {
    if (st.collapsed) openAll();
    if (st.colour.indexOf("ann:") === 0 && !D.annFields.some(function (f) { return "ann:" + f.col === st.colour; })) st.colour = "response";
    // 1. edge weights for the selection
    var eAll = edges.get({ filter: function (e) { return e.etype !== "kinase"; } }), W = {}, strength = {};
    eAll.forEach(function (e) { var a = weight(e, "EE"), b = weight(e, "RE"), w = st.arm === "EE" ? a[0] : st.arm === "RE" ? b[0] : a[0] - b[0];
      W[e.id] = { w: w, wE: a[0], wR: b[0], n: a[1] }; });
    var visible = function (e) { return st.typeOn[e.etype] && W[e.id].n > 0; };
    var wl = eAll.filter(visible).map(function (e) { return W[e.id].w; }), lim = q95(wl);
    eAll.forEach(function (e) { if (visible(e)) { var a = Math.abs(W[e.id].w); strength[e.from] = (strength[e.from] || 0) + a; strength[e.to] = (strength[e.to] || 0) + a; } });
    // 2. node values for the colour mode
    var nAll = nodes.get(), val = {}, vals = [];
    nAll.forEach(function (n) { var v = null;
      if (st.colour === "response") v = response(n.id);
      else if (st.colour === "phospho") { var pc = phosCells(n.id); v = (D.phos[n.id] ? pc : null); }
      else { var f = st.colour.substr(4), a = D.ann[n.id]; v = a ? a[f] : null; }
      val[n.id] = v; if (typeof v === "number") vals.push(v); });
    var vlim = q95(vals), vmax = Math.max.apply(null, vals.concat([1]));
    var smax = Math.max.apply(null, Object.keys(strength).map(function (k) { return strength[k]; }).concat([1e-9]));
    var keep = null; if (st.focus) { keep = {}; st.focus.forEach(function (i) { keep[i] = true; net.getConnectedNodes(i).forEach(function (j) { keep[j] = true; }); }); }
    var thr = st.thr, fset = {}, tally = {}; if (st.focus) st.focus.forEach(function (i) { fset[i] = true; });
    // 3. node styles + tooltips
    nodes.update(nAll.map(function (n) {
      var v = val[n.id], col = "#E6E6E6", border = "#555555", bw = 0.8, sc = sigCells(n.id), sig = sc.filter(function (r) { return r[4] < thr; });
      if (st.colour === "response") { if (v !== null) col = div3(v, vlim, st.arm === "ER" ? DIFF : RESP); if (sig.length) { border = "#000000"; bw = 3; } }
      else if (st.colour === "phospho") { if (v === null) col = "#E6E6E6"; else { var up = v.filter(function (r) { return r[4] < thr && r[2] > 0; }).length, dn = v.filter(function (r) { return r[4] < thr && r[2] < 0; }).length;
          col = up && dn ? "#8C510A" : up ? "#E66100" : dn ? "#5E3C99" : "#FFFFFF"; if (up || dn) { border = "#000000"; bw = 2; } } }
      else { var fld = D.annFields.filter(function (f) { return "ann:" + f.col === st.colour; })[0];
        if (v === null || v === undefined) { col = "#E6E6E6"; tally["no record"] = (tally["no record"] || 0) + 1; }
        else if (fld.binary) { col = v > 0 ? "#1B9E77" : "#FFFFFF"; var lb = v > 0 ? "yes" : "no"; tally[lb] = (tally[lb] || 0) + 1; }
        else { var bn = BINS.filter(function (b) { return v >= b[0] && v <= b[1]; })[0]; col = bn[3]; tally[bn[2]] = (tally[bn[2]] || 0) + 1; } }
      if (st.colour === "phospho") { var pl = v === null ? "not measured" : col === "#E66100" ? "up" : col === "#5E3C99" ? "down" : col === "#8C510A" ? "up and down" : "measured, none respond"; tally[pl] = (tally[pl] || 0) + 1; }
      if (st.colour === "response") { var rl = v === null ? "no data in selection" : sig.length ? "significant" : "not significant"; tally[rl] = (tally[rl] || 0) + 1; }
      var onF = !keep || keep[n.id], size = 6 + 16 * Math.sqrt((strength[n.id] || 0) / smax);
      // tooltip
      var t = "<b>" + esc(n.id) + "</b> (" + (n.shape === "triangle" ? "metabolite · " + esc(n.group) : "protein") + ")<br><i>" + ARMLAB[st.arm] + ", selection</i>: mean normalised response " + f3(response(n.id));
      if (sc.length) t += "<br>" + sc.slice(0, 9).map(function (r) { return (r[4] < thr ? "<b>" : "") + r[0] + " " + r[1] + " " + r[2] + ": logFC " + f3(r[3]) + ", adj. p " + f3(r[4]) + (r[4] < thr ? "</b>" : ""); }).join("<br>") + (sc.length > 9 ? "<br>… " + (sc.length - 9) + " more" : "");
      if (D.phos[n.id]) { var pcs = phosCells(n.id).filter(function (r) { return r[4] < thr; });
        t += "<br><u>MoTrPAC phosphosites</u>: " + D.phos[n.id].length + " measured, " + pcs.length + " respond in selection" + (pcs.length ? "<br>" + pcs.slice(0, 8).map(function (r) { return r[3].site + " (" + r[0] + " " + r[1] + ") logFC " + f3(r[2]) + ", adj. p " + f3(r[4]) + (r[3].kin ? " · kinase " + esc(r[3].kin) : "") + (r[3].xt ? " · also an O-glycosylation site" : ""); }).join("<br>") : ""); }
      var a = D.ann[n.id]; if (a) t += "<br><u>PTM (mnet)</u>: " + a.phosphosites + " phosphosites (" + a.kinase_sites + " with kinase)" + (a.is_kinase ? " · kinase (" + a.substrate_sites + " substrate sites)" : "") + " · " + a.glyco_sites + " glycosylation sites (N " + a.glyco_N_sites + ", O " + a.glyco_O_sites + ")" + (a.crosstalk_residues ? " · " + a.crosstalk_residues + " phospho = O-glyco residues" : "") +
        "<br><u>GlyGen</u>: " + a.glycans + " glycan structures" + (a.glyco_protein_level ? " · glycosylated (protein-level evidence)" : "") + " · " + a.mutations + " mutations · " + a.disease + " diseases · " + a.pathways + " pathways · " + a.publications + " publications";
      var mm = D.mods.filter(function (m) { return m.members.indexOf(n.id) >= 0; })[0]; if (mm) t += "<br>module " + mm.id;
      return { id: n.id, size: size, title: t, borderWidth: onF ? bw : 0.5, font: { color: onF ? "#1A1A1A" : "rgba(0,0,0,0.06)" },
               color: { background: onF ? col : "rgba(230,230,230,0.25)", border: onF ? border : "rgba(170,170,170,0.25)", highlight: { background: col, border: "#000000" }, hover: { background: col, border: "#000000" } } };
    }));
    // 3b. PTM tags per protein for the current selection
    var PIN = {EE: "#E41A1C", RE: "#377EB8", both: "#984EA3"}, pinTally = {EE: 0, RE: 0, both: 0, kp: 0, N: 0, O: 0, OG: 0, unk: 0, xt: 0};
    st.pins = {}; st.keep = keep;
    nAll.forEach(function (n) { var P = [], a = D.ann[n.id] || null;
      if (st.tags.mp && D.phos[n.id]) { var rs = [];
        D.phos[n.id].forEach(function (s) { var ee = false, re = false, mp = 1;
          on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) { var e = s.c[t + "|EE|" + h], r = s.c[t + "|RE|" + h];
            if (e && e[1] < thr) { ee = true; mp = Math.min(mp, e[1]); } if (r && r[1] < thr) { re = true; mp = Math.min(mp, r[1]); } }); });
          if (ee || re) rs.push({ kind: "P", arm: ee && re ? "both" : ee ? "EE" : "RE", p: mp, site: s.site }); });
        rs.sort(function (x, y) { return x.p - y.p; });
        rs.forEach(function (r) { pinTally[r.arm]++; });
        rs.slice(0, 6).forEach(function (r) { P.push({ kind: "P", fill: PIN[r.arm] }); });
        if (rs.length > 6) P.push({ kind: "more", n: rs.length - 6 }); }
      if (a && st.tags.kp && a.phosphosites > 0) { P.push({ kind: "Pdb", n: a.phosphosites }); pinTally.kp++; }
      if (a) {
        if (st.tags.N && a.glyco_N_sites > 0) { P.push({ kind: "sq", fill: "#0072BC", n: a.glyco_N_sites }); pinTally.N++; }
        if (st.tags.O && a.glyco_O_sites > 0) { P.push({ kind: "sq", fill: "#FFD400", n: a.glyco_O_sites }); pinTally.O++; }
        if (st.tags.OG && a.glyco_OGlcNAc_sites > 0) { P.push({ kind: "sqo", fill: "#0072BC", n: a.glyco_OGlcNAc_sites }); pinTally.OG++; }
        if (st.tags.unk && !(a.glyco_sites > 0) && a.glyco_protein_level > 0) { P.push({ kind: "sqh" }); pinTally.unk++; } }
      if (a && st.tags.xt && a.crosstalk_residues > 0) { P.push({ kind: "star", n: a.crosstalk_residues }); pinTally.xt++; }
      if (P.length) st.pins[n.id] = P; });
    st.pinTally = pinTally;
    // 4. edge styles + tooltips
    edges.update(eAll.map(function (e) { var r = W[e.id], vis = visible(e), onF = !st.focus || fset[e.from] || fset[e.to];
      var c = st.arm === "ER" ? div3(r.w, lim, ["#2166AC", "#D9D9D9", "#B2182B"]) : TYPE_COL[e.etype];
      var d = st.arm === "ER" ? (e.etype === "protein - protein" ? "solid" : e.etype === "metabolite - metabolite" ? "mm" : "mp") : (r.w < 0 ? "neg" : "solid");
      return { id: e.id, hidden: !vis, width: 0.6 + 5 * Math.min(Math.abs(r.w) / lim, 1), dashes: DASH[d],
               color: { color: onF ? c : "rgba(200,200,200,0.12)", highlight: c, hover: c, opacity: onF ? 0.85 : 1 },
               title: "<b>" + esc(e.from) + " — " + esc(e.to) + "</b><br>" + e.etype + "<br>selection (" + r.n + " terms): w_EE " + f3(r.wE) + " · w_RE " + f3(r.wR) + " · w_EE − w_RE " + f3(r.wE - r.wR) + (e.info ? "<br>" + e.info : "") }; }));
    // kinase edges: shown if switched on; tooltip lists the substrate's MoTrPAC response at those sites
    edges.update(edges.get({ filter: function (e) { return e.etype === "kinase"; } }).map(function (e) {
      var k = e.kin, sites = k.sites.split(","), resp = (D.phos[k.substrate] || []).filter(function (s) { return sites.indexOf(s.site) >= 0; });
      var rtxt = resp.map(function (s) { var best = null; on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) { var r = s.c[t + "|" + st.arm + "|" + h]; if (r && (!best || r[1] < best[1])) best = [r[0], r[1], t, h]; }); });
        return s.site + ": " + (best ? "logFC " + f3(best[0]) + ", adj. p " + f3(best[1]) + " (" + best[2] + " " + best[3] + ")" : "not measured in selection"); });
      var any = resp.some(function (s) { return on(st.tis).some(function (t) { return on(st.times).some(function (h) { var r = s.c[t + "|" + st.arm + "|" + h]; return r && r[1] < st.thr; }); }); });
      return { id: e.id, hidden: !st.kinase, width: any ? 3 : 1.6, color: { color: any ? "#D7301F" : "#111111", highlight: "#D7301F" },
               title: "<b>" + esc(k.kinase) + " → " + esc(k.substrate) + "</b> (kinase → substrate, OmniPath via mnet: " + esc(k.source) + ")<br>sites: " + esc(k.sites) +
                      "<br>MoTrPAC at these sites (" + ARMLAB[st.arm] + "): " + (rtxt.length ? "<br>" + rtxt.join("<br>") : "not measured") + (any ? "<br><b>a substrate site responds (red edge)</b>" : "") }; }));
    // 5. legend: rebuilt for whatever is shown (node colour mode, outline, edges, sizes, arrows), with counts
    var sw = function (c, lab, n, shape) { return "<div class='hk-row'><span class='hk-swatch' style='background:" + c + (shape === "ring" ? ";border:3px solid #000;background:#FFF" : "") + "'></span>" + lab + (n !== undefined ? " <span class='hk-n'>(" + n + ")</span>" : "") + "</div>"; };
    var bar3 = function (pal, lo, mid, hi, loLab, hiLab) { return "<div class='hk-bar3' style='background:linear-gradient(90deg," + pal.join(",") + ")'></div><div class='hk-ticks'><span>" + lo + "</span><span>" + mid + "</span><span>" + hi + "</span></div><div class='hk-ticks hk-sub'><span>" + loLab + "</span><span>" + hiLab + "</span></div>"; };
    var G = "<div class='hk-lt'>Nodes</div>";
    if (st.colour === "response") {
      G += "<div class='hk-ls'>fill = mean normalised response (" + ARMLAB[st.arm] + ", selection)</div>" + (st.arm === "ER" ? bar3(DIFF, "−" + f3(vlim), "0", "+" + f3(vlim), "higher in resistance", "higher in endurance") : bar3(RESP, "−" + f3(vlim), "0", "+" + f3(vlim), "down vs control", "up vs control")) +
           sw("#FFF", "black outline: adj. p &lt; " + thr + " in a selected cell", tally["significant"] || 0, "ring") + sw("#E6E6E6", "grey: no data in the selection", tally["no data in selection"] || 0);
    } else if (st.colour === "phospho") {
      G += "<div class='hk-ls'>fill = MoTrPAC phosphosites responding (adj. p &lt; " + thr + ", " + ARMLAB[st.arm] + ", selected muscle / adipose times)</div>" +
           [["#E66100", "up"], ["#5E3C99", "down"], ["#8C510A", "up and down"], ["#FFFFFF", "measured, none respond"], ["#E6E6E6", "not measured"]].map(function (x) { return sw(x[0], x[1], tally[x[1]] || 0); }).join("");
    } else { var fl = D.annFields.filter(function (f) { return "ann:" + f.col === st.colour; })[0];
      G += "<div class='hk-ls'>fill = " + esc(fl.label) + " (per protein)</div>" + (fl.binary ? sw("#1B9E77", "yes", tally["yes"] || 0) + sw("#FFFFFF", "no", tally["no"] || 0)
           : BINS.map(function (b) { return sw(b[3], b[2], tally[b[2]] || 0); }).join("")) + sw("#E6E6E6", "no record (incl. metabolites)", tally["no record"] || 0);
    }
    G += "<div class='hk-ls'>circle = protein, triangle = metabolite; size = strength (sum |w|) over the shown edges</div>";
    var T = st.pinTally || {}, pin = function (bg, fg, txt, lab, n, cls) { return "<div class='hk-row'><span class='hk-pin " + (cls || "") + "' style='background:" + bg + ";color:" + fg + "'>" + txt + "</span>" + lab + (n !== undefined ? " <span class='hk-n'>(" + n + ")</span>" : "") + "</div>"; };
    if (on(st.tags).length) {
      G += "<div class='hk-lt'>PTM tags (stalks on proteins)</div>";
      if (st.tags.mp) G += "<div class='hk-ls'>MoTrPAC phosphosites responding (adj. p &lt; " + thr + ", selected tissues / times; one pin per site, up to 6, then +n)</div>" +
        pin("#E41A1C", "#FFF", "P", "after endurance only", T.EE + " sites") + pin("#377EB8", "#FFF", "P", "after resistance only", T.RE + " sites") + pin("#984EA3", "#FFF", "P", "after both", T.both + " sites");
      if (st.tags.kp) G += pin("#D0D0D0", "#222", "P", "known phosphosites (mnet: UniProt + OmniPath; number = sites)", T.kp + " proteins");
      if (st.tags.N || st.tags.O || st.tags.OG || st.tags.unk) G += "<div class='hk-ls'>glycosylation, SNFG symbols (UniProt via mnet; number = sites)</div>" +
        (st.tags.N ? pin("#0072BC", "#FFF", "", "N-linked (GlcNAc)", T.N + " proteins", "hk-sq") : "") +
        (st.tags.O ? pin("#FFD400", "#222", "", "O-linked, mucin-type (GalNAc)", T.O + " proteins", "hk-sq") : "") +
        (st.tags.OG ? pin("#0072BC", "#FFF", "●", "O-GlcNAc", T.OG + " proteins", "hk-sq") : "") +
        (st.tags.unk ? pin("#FFFFFF", "#0072BC", "", "glycosylated, site unknown (GlyGen protein-level)", T.unk + " proteins", "hk-sq hk-dash") : "");
      if (st.tags.xt) G += "<div class='hk-row'><span class='hk-star'>★</span>a MoTrPAC phosphosite that is also an O-glycosylation site <span class='hk-n'>(" + T.xt + " proteins)</span></div>";
    }
    G += "<div class='hk-lt'>Edges</div>" + (st.arm === "ER" ? "<div class='hk-ls'>colour = w_EE − w_RE; width = |difference|</div>" + bar3(["#2166AC", "#D9D9D9", "#B2182B"], "−" + f3(lim), "0", "+" + f3(lim), "higher in resistance", "higher in endurance") + "<div class='hk-ls'>solid = protein–protein, long dash = metabolite–metabolite, dotted = metabolite–protein</div>"
         : "<div class='hk-ls'>width = |w| (up to " + f3(lim) + "); dashed = negative weight</div>" + D.types.map(function (t) { return "<div class='hk-row'><span class='hk-line' style='background:" + TYPE_COL[t] + "'></span>" + t + "</div>"; }).join(""));
    if (st.kinase) G += "<div class='hk-row'><span class='hk-line' style='background:#111'></span>→ kinase → substrate (OmniPath via mnet)</div><div class='hk-row'><span class='hk-line' style='background:#D7301F'></span>→ a substrate site responds in the selection</div>";
    legendBox.innerHTML = G;
    legend.innerHTML = "<b>Selection</b>: " + ARMLAB[st.arm] + " · " + on(st.omes).map(function (o) { return OMES[o]; }).join(", ") + " · " + on(st.tis).join(", ") + " · " + on(st.times).join(", ") + " · adj. p &lt; " + thr;
    // 6. modules: significance for the selection, and the selected module's table
    var ms = q(".hk-mod"), cur = ms.value;
    ms.innerHTML = "<option value=''>none</option>" + D.mods.map(function (m) { var best = modCells(m).reduce(function (b, r) { return (!b || r[5] < b[5]) ? r : b; }, null);
      return "<option value='" + m.id + "'" + (m.id === cur ? " selected" : "") + ">" + esc(m.name) + " (" + m.n_prot + " prot, " + m.n_met + " met)" + (best && best[5] < thr ? " ★ FDR " + f3(best[5]) : "") + "</option>"; }).join("");
    var M = D.mods.filter(function (m) { return m.id === cur; })[0];
    if (M) { var rows = modCells(M);
      panel.innerHTML = "<b>" + esc(M.name) + "</b> (" + M.id + "): " + M.members.length + " nodes<br><u>Pathways over-represented among its members</u> (ORA vs our 471 genes / 450 metabolites, BH adj. p): " +
        (M.ora.length ? "<ul class='hk-ul'>" + M.ora.map(function (o) { return "<li>" + esc(o[0]) + " — adj. p " + f3(o[1]) + ", " + o[2] + " members (" + esc(o[3]) + ")</li>"; }).join("") + "</ul>" : " none significant.<br>") +
        "<u>Exercise response of the module</u> — CAMERA-PR (MoTrPAC run_cameraPR; competitive, within each ome) for " + ARMLAB[st.arm] + " in the selection:" +
        (rows.length ? "<table><tr><th>tissue</th><th>ome</th><th>time</th><th>members tested</th><th>direction</th><th>z</th><th>FDR</th></tr>" + rows.map(function (r) { return "<tr" + (r[5] < thr ? " class='hk-sig'" : "") + "><td>" + r[0] + "</td><td>" + OMES[r[1]] + "</td><td>" + r[2] + "</td><td>" + r[6] + "</td><td>" + r[4] + "</td><td>" + f3(r[3]) + "</td><td>" + f3(r[5]) + "</td></tr>"; }).join("") + "</table>" : " no tested cell (fewer than 5 members measured in the selected omes / tissues).");
    } else panel.innerHTML = "";
    bar.querySelectorAll(".hk-arm").forEach(function (b) { b.checked = b.value === st.arm; });
    if (st.collapsed) clusterAll();
  }
  function modCells(m) { var out = []; on(st.omes).forEach(function (o) { on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) {
    var r = m.cam[t + "|" + o + "|" + st.arm + "|" + h]; if (r) out.push([t, o, h, r[0], r[1], r[2], r[3]]); }); }); }); return out; }
  // ---------- wiring ----------
  var chk = function (cls, obj) { bar.querySelectorAll(cls).forEach(function (c) { c.onchange = function () { obj[c.value] = c.checked ? 1 : 0; apply(); }; }); };
  chk(".hk-ome", st.omes); chk(".hk-tis", st.tis); chk(".hk-time", st.times); chk(".hk-tag", st.tags);
  bar.querySelectorAll(".hk-arm").forEach(function (r) { r.onchange = function () { st.arm = r.value; apply(); }; });
  q(".hk-thr").onchange = function () { var v = parseFloat(this.value); if (v > 0 && v <= 1) { st.thr = v; apply(); } };
  q(".hk-colour").onchange = function () { st.colour = this.value; apply(); };
  q(".hk-mod").onchange = function () { var m = D.mods.filter(function (x) { return x.id === q(".hk-mod").value; })[0];
    st.focus = m ? m.members.slice() : null; apply(); if (m) net.fit({ nodes: m.members, animation: { duration: 600 } }); };
  bar.querySelectorAll(".hk-type").forEach(function (c) { c.onchange = function () { st.typeOn[c.value] = c.checked; apply(); }; });
  if (q(".hk-kin")) q(".hk-kin").onchange = function () { st.kinase = this.checked; apply(); };
  q(".hk-find").onchange = function () { var id = this.value; if (nodes.get(id) === null) return;
    if (st.collapsed) { st.collapsed = false; openAll(); var b = q(".hk-col"); if (b) b.textContent = "Collapse classes"; }
    st.focus = [id]; apply(); net.selectNodes([id]); net.focus(id, { scale: 1.3, animation: { duration: 600 } }); };
  if (q(".hk-cls")) q(".hk-cls").onchange = function () { var c = this.value; if (!c) { st.focus = null; apply(); net.fit({ animation: true }); return; }
    st.focus = nodes.getIds({ filter: function (n) { return n.group === c; } }); apply(); net.fit({ nodes: st.focus, animation: { duration: 600 } }); };
  if (q(".hk-hull")) q(".hk-hull").onchange = function () { st.hulls = this.checked; net.redraw(); };
  if (q(".hk-col")) q(".hk-col").onclick = function () { if (st.collapsed) { st.collapsed = false; openAll(); this.textContent = "Collapse classes"; } else { st.collapsed = true; clusterAll(); this.textContent = "Expand classes"; } net.redraw(); };
  q(".hk-reset").onclick = function () { location.reload(); };
  net.on("click", function (p) { if (p.nodes.length && !net.isCluster(p.nodes[0])) { st.focus = [p.nodes[0]]; apply(); }
    else if (!p.nodes.length && !p.edges.length && st.focus) { st.focus = null; q(".hk-mod").value = ""; apply(); } });
  net.on("doubleClick", function (p) { if (p.nodes.length && net.isCluster(p.nodes[0])) net.openCluster(p.nodes[0]); });
  window.hkState = st; window.hkWeight = weight; window.hkApply = apply;
  apply();
}
)---"

# Page styling for the toolbar, legend and module panel.
CSS <- tags$style(HTML("
  body { font-family: Helvetica, Arial, sans-serif; color: #1A1A1A; margin: 10px 16px; }
  .hk-bar { font-size: 12.5px; padding: 4px 0 6px; line-height: 2.0; }
  .hk-bar button { font-size: 12px; margin-right: 3px; padding: 1px 8px; border: 1px solid #999; background: #F7F7F7; border-radius: 3px; cursor: pointer; }
  .hk-bar input.hk-find { width: 170px; font-size: 12px; } .hk-bar select { font-size: 12px; max-width: 330px; } .hk-bar label { margin-right: 5px; }
  .hk-sep { display: inline-block; width: 12px; }
  .hk-legend { font-size: 11.5px; padding: 6px 0 2px; line-height: 1.6; } .hk-help { font-size: 11px; color: #666; padding-top: 4px; }
  .hk-panel { font-size: 12px; padding: 6px 0; } .hk-panel table { border-collapse: collapse; margin-top: 4px; }
  .hk-panel td, .hk-panel th { border: 1px solid #DDD; padding: 2px 8px; text-align: left; } .hk-panel tr.hk-sig td { font-weight: bold; background: #FFF4D6; }
  .hk-grad { display: inline-block; width: 80px; height: 9px; vertical-align: middle; border: 1px solid #999; margin: 0 4px; }
  .hk-sw2 { display: inline-block; width: 11px; height: 11px; vertical-align: middle; border: 1px solid #777; margin: 0 3px 0 6px; border-radius: 6px; }
  .hk-wrap { display: flex; gap: 12px; align-items: flex-start; }
  .hk-legbox { flex: 0 0 290px; max-height: 800px; overflow-y: auto; background: #FFFFFF; border: 1px solid #BBB; border-radius: 4px; padding: 6px 10px; font-size: 11px; margin-top: 50px; }
  .hk-lt { font-weight: bold; font-size: 12px; margin-top: 4px; } .hk-ls { color: #444; margin: 2px 0; } .hk-row { margin: 1px 0; } .hk-n { color: #777; }
  .hk-swatch { display: inline-block; width: 12px; height: 12px; border: 1px solid #777; border-radius: 7px; vertical-align: middle; margin-right: 6px; box-sizing: border-box; }
  .hk-line { display: inline-block; width: 22px; height: 3px; vertical-align: middle; margin-right: 6px; }
  .hk-bar3 { height: 10px; border: 1px solid #999; margin-top: 3px; } .hk-ticks { display: flex; justify-content: space-between; font-size: 10px; color: #333; } .hk-sub { color: #777; }
  .hk-ul { margin: 2px 0 4px 18px; padding: 0; }
  .hk-pin { display: inline-block; width: 13px; height: 13px; border-radius: 7px; border: 1px solid #222; font-size: 8.5px; font-weight: bold; text-align: center; line-height: 13px; vertical-align: middle; margin-right: 6px; box-sizing: border-box; }
  .hk-pin.hk-sq { border-radius: 0; font-size: 7px; } .hk-pin.hk-dash { border: 1.5px dashed #0072BC; }
  .hk-star { color: #FFD700; font-size: 14px; -webkit-text-stroke: 0.6px #222; margin-right: 5px; vertical-align: middle; }
  div.vis-tooltip { font-family: Helvetica, Arial, sans-serif; font-size: 11.5px; line-height: 1.4; max-width: 460px; white-space: normal; }
"))

# Build and save one interactive page.
page <- function(P, net, types, classes, hulls, file, title, subtitle) {
  vn <- P$N[, .(id = node, label = node, group = class, shape = fifelse(node_type == "metabolite", "triangle", "dot"),
                x = x * W, y = (1 - y) * H, physics = FALSE, hull, size = sizeEE, color = colEE)]
  ve <- P$E[, .(id = paste0("e", seq_len(.N)), from = a, to = b, etype = as.character(edge_type), info = gsub("\n", " ", fifelse(is.na(info), "", info)), w_EE_ref = w_EE)]
  cfg <- list(json = as.character(page_data(P, net, types)), static = list(classes = I(classes), hulls = hulls))
  w <- visNetwork(as.data.frame(vn), as.data.frame(ve), width = "100%", height = "760px",
                  main = list(text = title, style = "font-family:Helvetica;font-weight:bold;font-size:17px;text-align:left"),
                  submain = list(text = subtitle, style = "font-family:Helvetica;font-size:12px;color:#555;text-align:left")) |>
    visPhysics(enabled = FALSE) |> visEdges(smooth = FALSE, selectionWidth = 1.5) |>
    visNodes(borderWidth = 0.8, font = list(size = 13, face = "Helvetica"), scaling = list(label = list(enabled = TRUE, min = 13, max = 13, drawThreshold = 9))) |>
    visInteraction(navigationButtons = TRUE, hover = TRUE, tooltipDelay = 80, hideEdgesOnDrag = TRUE, multiselect = TRUE) |>
    onRender(JS, data = cfg)
  w <- prependContent(w, CSS)
  # a fixed widget ID (visNetwork ignores elementId as an argument), so reruns write byte-identical pages
  w$elementId <- paste0("hk-", net)
  saveWidget(w, file.path(HTML_DIR, file), selfcontained = TRUE, title = title)
  unlink(file.path(HTML_DIR, sub("\\.html$", "_files", file)), recursive = TRUE)
  message("-> ", file.path(HTML_DIR, file))
}

# Cytoscape export: .cyjs (Cytoscape.js JSON with positions) and plain tables for one network.
cyto <- function(P, stem, name) {
  # node data: identifiers, attributes, and ready-made per-view colours / sizes / shape for the styles
  nd <- P$N[, .(id = node, name = node, shared_name = node, node_type, class, resp_EE, resp_RE, strength_EE, strength_RE, strength_diff,
                degree, col_EE = colEE, col_RE = colRE, col_diff = colD, size_EE = 2 * sizeEE, size_RE = 2 * sizeRE, size_diff = 2 * sizeD,
                cy_shape = fifelse(node_type == "metabolite", "TRIANGLE", "ELLIPSE"))]
  # edge data: identifiers, weights and ready-made per-view colours / widths / line types
  line <- c(solid = "SOLID", neg = "EQUAL_DASH", mm = "LONG_DASH", mp = "DOT")
  ed <- P$E[, .(id = paste0("e", seq_len(.N)), source = a, target = b, name = paste(a, "(interacts with)", b), interaction = "interacts with",
                edge_type = as.character(edge_type), w_EE, w_RE, w_diff, evidence = gsub("<[^>]+>", " ", info),
                col_EE = colEE, col_RE = colRE, col_diff = colD, width_EE = wdEE, width_RE = wdRE, width_diff = wdD,
                line_EE = line[dashEE], line_RE = line[dashRE], line_diff = line[dashD])]
  # plain tables
  fwrite(nd, file.path(CY_DIR, sprintf("17_%s_nodes.csv", stem))); fwrite(ed, file.path(CY_DIR, sprintf("17_%s_edges.csv", stem)))
  # Cytoscape.js JSON: one element per node (with position) and per edge
  nodes_js <- lapply(seq_len(nrow(nd)), function(i) list(data = as.list(nd[i]), position = list(x = P$N$x[i] * W, y = (1 - P$N$y[i]) * H)))
  edges_js <- lapply(seq_len(nrow(ed)), function(i) list(data = as.list(ed[i])))
  js <- list(format_version = "1.0", generated_by = "17_interactive_networks.R", target_cytoscapejs_version = "~2.1",
             data = list(name = name, shared_name = name), elements = list(nodes = nodes_js, edges = edges_js))
  write_json(js, file.path(CY_DIR, sprintf("17_%s_network.cyjs", stem)), auto_unbox = TRUE, digits = NA, pretty = FALSE, na = "null")
}

# ---- joint network (steps 14, 15) -----------------------------------------------------------------------
jn <- fread(file.path(OUT, "14_joint_nodes.csv")); je <- fread(file.path(OUT, "14_joint_edges.csv"))
cl <- fread(file.path(OUT, "15_class_layout.csv"))
# evidence for each edge: STRING score (protein - protein), link rule (metabolite - metabolite), Rhea reactions (metabolite - protein)
w3 <- fread(file.path(OUT, "03_weighted_edges.csv")); w6 <- fread(file.path(OUT, "06_metabolite_edges.csv"))
rh <- unique(fread(file.path(OUT, "05_metabolite_protein_links.csv"))[, .(metabolite, gene_symbol, n_reactions, example_reactions)], by = c("metabolite", "gene_symbol"))
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "~")
ev <- rbind(w3[, .(k = key(symbol_a, symbol_b), info = paste0("STRING combined score ", combined_score))],
            w6[, .(k = key(metabolite_a, metabolite_b), info = paste0(esc(class), " · ", link_type,
                   fifelse(shared_proteins != "", paste0("<br>shared enzymes: ", esc(gsub(";", ", ", shared_proteins))), ""),
                   fifelse(!is.na(string_protein_pairs) & string_protein_pairs != "", paste0("<br>STRING-linked enzymes: ", esc(gsub(";", ", ", string_protein_pairs))), "")))],
            rh[, .(k = key(metabolite, gene_symbol), info = sprintf("Rhea: %d reaction(s), e.g. %s", n_reactions, esc(gsub(";", ", ", example_reactions))))])
# joint nodes with the class-grouped layout; outlines group metabolites by class
JN <- jn[, .(node, node_type, class, resp_EE, resp_RE)][cl[, .(node, x, y)], on = "node"]
JN[, hull := fifelse(node_type == "metabolite", class, "")]
JE <- je[, .(a = node_a, b = node_b, edge_type, w_EE, w_RE, w_diff)][, info := ev$info[match(key(a, b), ev$k)]]
# safety checks: same sizes as step 14, every node placed, weights unchanged
stopifnot(nrow(JN) == nrow(jn), !anyNA(JN$x), nrow(JE) == nrow(je), !anyNA(JE$info))
PJ <- prepare(JN, JE)
jcls <- JN[node_type == "metabolite", .N, by = class][order(-N), class]
page(PJ, "joint", names(TYPE_COL), jcls, lapply(jcls, function(c) list(key = c, label = c)), "17a_joint_network.html",
     "Joint protein-metabolite network: endurance, resistance and their difference",
     sprintf("%d nodes (%d proteins, %d metabolites) · %d edges: STRING v12 >= 700 and Rhea catalysis / transport (team mnet resource), shared/STRING-linked enzyme + same RefMet super class · layout as figures 15a / 15b", nrow(PJ$N), sum(PJ$N$node_type == "protein"), sum(PJ$N$node_type == "metabolite"), nrow(PJ$E)))
cyto(PJ, "joint", "hackathon joint protein-metabolite network")

# ---- gene network (steps 1, 3, 10) ------------------------------------------------------------------------
ge <- fread(file.path(OUT, "01_nodes_EE.csv")); gr <- fread(file.path(OUT, "01_nodes_RE.csv")); gl <- fread(file.path(OUT, "10_layout_genes.csv"))
GN <- data.table(node = ge$gene_symbol, node_type = "protein", class = "protein",
                 resp_EE = rowMeans(as.matrix(ge[, -(1:2)]), na.rm = TRUE), resp_RE = rowMeans(as.matrix(gr[, -(1:2)]), na.rm = TRUE))
# only genes with an edge (as in figure 10a), placed as in 10a / 11a
GN <- GN[gl, on = "node"][, hull := ""]
GE <- w3[, .(a = symbol_a, b = symbol_b, edge_type = "protein - protein", w_EE, w_RE, w_diff = w_EE - w_RE, info = paste0("STRING combined score ", combined_score))]
stopifnot(nrow(GE) == nrow(w3), all(c(GE$a, GE$b) %in% GN$node), !anyNA(GN$resp_EE))
PG <- prepare(GN, GE)
page(PG, "gene", "protein - protein", character(0), list(), "17b_gene_network.html",
     "Gene network: endurance, resistance and their difference",
     sprintf("%d genes with an edge (of 471) · %d STRING edges (combined score >= 700) · layout as figures 10a / 11a", nrow(GN), nrow(GE)))
cyto(PG, "gene", "hackathon gene network")

# ---- metabolite network (steps 1b, 6, 10) -----------------------------------------------------------------
me <- fread(file.path(OUT, "01b_metab_nodes_EE.csv")); mr <- fread(file.path(OUT, "01b_metab_nodes_RE.csv")); ml <- fread(file.path(OUT, "10_layout_metabolites.csv"))
ids <- fread(file.path(OUT, "01c_metabolite_ids.csv"))
MN <- data.table(node = me$metabolite, node_type = "metabolite", class = ids$super_class[match(me$metabolite, ids$metabolite)],
                 resp_EE = rowMeans(as.matrix(me[, -1])), resp_RE = rowMeans(as.matrix(mr[, -1])))[ml, on = "node"]
ME <- w6[, .(a = metabolite_a, b = metabolite_b, edge_type = "metabolite - metabolite", w_EE, w_RE, w_diff = w_EE - w_RE)][, info := ev$info[match(key(a, b), ev$k)]]
# outlines: one per connected group (a class can form several groups in this layout, as labelled in figure 10b)
comp <- igraph::components(igraph::graph_from_data_frame(ME[, .(a, b)], directed = FALSE, vertices = MN[, .(node)]))$membership
MN[, hull := paste(class, comp[node], sep = " | ")]
stopifnot(nrow(ME) == nrow(w6), all(c(ME$a, ME$b) %in% MN$node), !anyNA(MN$class))
PM <- prepare(MN, ME)
mh <- unique(MN[, .(hull, class)])
mcls <- MN[, .N, by = class][order(-N), class]
page(PM, "metabolite", "metabolite - metabolite", mcls, lapply(seq_len(nrow(mh)), function(i) list(key = mh$hull[i], label = mh$class[i])),
     "17c_metabolite_network.html", "Metabolite network: endurance, resistance and their difference",
     sprintf("%d metabolites with an edge (of 450) · %d edges: shared or STRING-linked Rhea enzyme + same RefMet super class · layout as figures 10b / 11b", nrow(MN), nrow(ME)))
cyto(PM, "metabolite", "hackathon metabolite network")

# ---- Cytoscape styles (one file, three styles; passthrough mappings of the ready-made columns) ------------
vp <- function(name, default, attr, type) sprintf('      <visualProperty name="%s" default="%s"><passthroughMapping attributeName="%s" attributeType="%s"/></visualProperty>', name, default, attr, type)
style <- function(label, v) paste0('  <visualStyle name="hackathon ', label, '">\n    <network>\n      <visualProperty name="NETWORK_BACKGROUND_PAINT" default="#FFFFFF"/>\n    </network>\n',
  '    <node>\n      <dependency name="nodeSizeLocked" value="true"/>\n',
  paste(c(vp("NODE_FILL_COLOR", "#CCCCCC", paste0("col_", v), "string"), vp("NODE_SIZE", "20", paste0("size_", v), "float"),
          vp("NODE_SHAPE", "ELLIPSE", "cy_shape", "string"), vp("NODE_LABEL", "", "name", "string"),
          '      <visualProperty name="NODE_BORDER_PAINT" default="#404040"/>', '      <visualProperty name="NODE_BORDER_WIDTH" default="0.8"/>',
          '      <visualProperty name="NODE_LABEL_FONT_SIZE" default="9"/>'), collapse = "\n"),
  '\n    </node>\n    <edge>\n      <dependency name="arrowColorMatchesEdge" value="true"/>\n',
  paste(c(vp("EDGE_STROKE_UNSELECTED_PAINT", "#8C8C8C", paste0("col_", v), "string"), vp("EDGE_WIDTH", "1", paste0("width_", v), "float"),
          vp("EDGE_LINE_TYPE", "SOLID", paste0("line_", v), "string"), '      <visualProperty name="EDGE_TRANSPARENCY" default="200"/>'), collapse = "\n"),
  '\n    </edge>\n  </visualStyle>')
writeLines(c('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<vizmap id="VizMap-hackathon-17" documentVersion="3.1">',
             style("EE", "EE"), style("RE", "RE"), style("difference", "diff"), '</vizmap>'), file.path(CY_DIR, "17_cytoscape_styles.xml"))
message("-> ", CY_DIR)

# Show the network sizes.
print(data.table(network = c("joint", "genes", "metabolites"), nodes = c(nrow(PJ$N), nrow(PG$N), nrow(PM$N)), edges = c(nrow(PJ$E), nrow(PG$E), nrow(PM$E))))
