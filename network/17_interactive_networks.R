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
KIN <- fread(file.path(OUT, "17_kinase_edges.csv")); ANN <- fread(file.path(OUT, "17_glygen_protein_annotation.csv"))
for (f in c("17_node_cell_stats.csv", "17_phospho_site_stats.csv", "17_modules.csv", "17_module_camera.csv", "17_kinase_edges.csv"))
  if (!file.exists(file.path(OUT, f))) stop("run network/17_filter_stats.R first (", f, " missing)")
# Annotation fields offered in the "colour nodes by" menu: column, label, kind (count / binary).
ANN_FIELDS <- list(list(col = "glygen_phosphosites", label = "GlyGen: phosphosites"), list(col = "glygen_kinase_sites", label = "GlyGen: phosphosites with a known kinase"),
                   list(col = "glycosylated", label = "GlyGen: glycosylated (any)", binary = TRUE), list(col = "glyco_sites", label = "GlyGen: glycosylation sites"),
                   list(col = "glyco_N_sites", label = "GlyGen: N-linked sites"), list(col = "glyco_O_sites", label = "GlyGen: O-linked sites"),
                   list(col = "glycans", label = "GlyGen: glycan structures at sites"), list(col = "crosstalk_residues", label = "phospho = O-glyco residues (crosstalk)"),
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
    list(id = x$module[1], members = I(x$node), n_prot = sum(x$node_type == "protein"), n_met = sum(x$node_type == "metabolite"),
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
  var st = { omes: {}, tis: {adipose: 1, blood: 1, muscle: 1}, times: {"0.5h": 1, "4h": 1, "24h": 1}, arm: "EE", thr: 0.05, colour: "response",
             focus: null, collapsed: false, hulls: S.hulls.length > 0, typeOn: {}, kinase: false, module: "" };
  D.omes.forEach(function (o) { st.omes[o] = 1; }); D.types.forEach(function (t) { st.typeOn[t] = true; });
  var DASH = { solid: false, neg: [5, 5], mm: [10, 5], mp: [2, 4], kin: [3, 3] };
  var TYPE_COL = {"protein - protein": "#8C8C8C", "metabolite - metabolite": "#1B7837", "metabolite - protein": "#8C510A"};
  var RESP = ["#6A3D9A", "#FFFFFF", "#E66100"], DIFF = ["#2166AC", "#F2F2F2", "#B2182B"], SEQ = ["#F7FCF5", "#74C476", "#00441B"];
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
       "<option value='phospho'>MoTrPAC phosphosites (filtered)</option>" +
       D.annFields.map(function (f) { return "<option value='ann:" + f.col + "'>" + esc(f.label) + "</option>"; }).join("") + "</select>" +
       "<span class='hk-sep'></span><b>Module</b> <select class='hk-mod'><option value=''>none</option></select>" +
       "<span class='hk-sep'></span><b>Find</b> <input class='hk-find' list='" + el.id + "-dl' placeholder='gene or metabolite'>" +
       "<datalist id='" + el.id + "-dl'>" + nodes.getIds().sort().map(function (i) { return "<option value=\"" + String(i).replace(/"/g, "&quot;") + "\">"; }).join("") + "</datalist>";
  if (S.classes.length) h += "<span class='hk-sep'></span><b>Class</b> <select class='hk-cls'><option value=''>all</option>" + S.classes.map(function (c) { return "<option>" + esc(c) + "</option>"; }).join("") + "</select>";
  h += "</div><div><b>Edges</b> " + D.types.map(function (t) { return cb("hk-type", t, t, 1); }).join(" ") +
       (D.kin.length ? " " + cb("hk-kin", "1", "kinase → substrate (GlyGen; " + D.kin.length + ")", 0) : "") +
       (S.hulls.length ? "<span class='hk-sep'></span>" + cb("hk-hull", "1", "class outlines", 1) + " <button class='hk-col'>Collapse classes</button>" : "") +
       "<span class='hk-sep'></span><button class='hk-reset'>Reset</button></div>";
  bar.innerHTML = h;
  var legend = document.createElement("div"); legend.className = "hk-legend";
  var panel = document.createElement("div"); panel.className = "hk-panel";
  var help = document.createElement("div"); help.className = "hk-help";
  help.innerHTML = "Filters recompute node colours, significance outlines and edge weights (dot products over the selected dimensions only) · click a node to highlight it and its neighbours, empty space to clear · hover for values" + (S.classes.length ? " · double-click a collapsed class to open it" : "");
  el.parentNode.insertBefore(bar, el); el.parentNode.appendChild(legend); el.parentNode.appendChild(panel); el.parentNode.appendChild(help);
  var q = function (s) { return bar.querySelector(s); };
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
    var thr = st.thr, fset = {}; if (st.focus) st.focus.forEach(function (i) { fset[i] = true; });
    // 3. node styles + tooltips
    nodes.update(nAll.map(function (n) {
      var v = val[n.id], col = "#E6E6E6", border = "#555555", bw = 0.8, sc = sigCells(n.id), sig = sc.filter(function (r) { return r[4] < thr; });
      if (st.colour === "response") { if (v !== null) col = div3(v, vlim, st.arm === "ER" ? DIFF : RESP); if (sig.length) { border = "#000000"; bw = 3; } }
      else if (st.colour === "phospho") { if (v === null) col = "#E6E6E6"; else { var up = v.filter(function (r) { return r[4] < thr && r[2] > 0; }).length, dn = v.filter(function (r) { return r[4] < thr && r[2] < 0; }).length;
          col = up && dn ? "#8C510A" : up ? "#E66100" : dn ? "#5E3C99" : "#FFFFFF"; if (up || dn) { border = "#000000"; bw = 2; } } }
      else { if (v === null || v === undefined) col = "#E6E6E6"; else { var bin = D.annFields.filter(function (f) { return "ann:" + f.col === st.colour; })[0].binary;
          col = bin ? (v > 0 ? "#1B9E77" : "#FFFFFF") : (v > 0 ? seq3(Math.log1p(v) / Math.log1p(vmax), SEQ) : "#FFFFFF"); } }
      var onF = !keep || keep[n.id], size = 6 + 16 * Math.sqrt((strength[n.id] || 0) / smax);
      // tooltip
      var t = "<b>" + esc(n.id) + "</b> (" + (n.shape === "triangle" ? "metabolite · " + esc(n.group) : "protein") + ")<br><i>" + ARMLAB[st.arm] + ", selection</i>: mean normalised response " + f3(response(n.id));
      if (sc.length) t += "<br>" + sc.slice(0, 9).map(function (r) { return (r[4] < thr ? "<b>" : "") + r[0] + " " + r[1] + " " + r[2] + ": logFC " + f3(r[3]) + ", adj. p " + f3(r[4]) + (r[4] < thr ? "</b>" : ""); }).join("<br>") + (sc.length > 9 ? "<br>… " + (sc.length - 9) + " more" : "");
      if (D.phos[n.id]) { var pcs = phosCells(n.id).filter(function (r) { return r[4] < thr; });
        t += "<br><u>MoTrPAC phosphosites</u>: " + D.phos[n.id].length + " measured, " + pcs.length + " respond in selection" + (pcs.length ? "<br>" + pcs.slice(0, 8).map(function (r) { return r[3].site + " (" + r[0] + " " + r[1] + ") logFC " + f3(r[2]) + ", adj. p " + f3(r[4]) + (r[3].kin ? " · kinase " + esc(r[3].kin) : "") + (r[3].xt ? " · O-GlcNAc site" : ""); }).join("<br>") : ""); }
      var a = D.ann[n.id]; if (a) t += "<br><u>GlyGen</u>: " + a.glygen_phosphosites + " phosphosites (" + a.glygen_kinase_sites + " with kinase) · " + (a.glycosylated ? a.glyco_sites + " glycosylation sites, " + a.glycans + " glycans" : "no glycosylation record") + (a.crosstalk_residues ? " · " + a.crosstalk_residues + " phospho = O-glyco residues" : "") + "<br>" + a.mutations + " mutations · " + a.disease + " diseases · " + a.pathways + " pathways · " + a.publications + " publications";
      var mm = D.mods.filter(function (m) { return m.members.indexOf(n.id) >= 0; })[0]; if (mm) t += "<br>module " + mm.id;
      return { id: n.id, size: size, title: t, borderWidth: onF ? bw : 0.5, font: { color: onF ? "#1A1A1A" : "rgba(0,0,0,0.06)" },
               color: { background: onF ? col : "rgba(230,230,230,0.25)", border: onF ? border : "rgba(170,170,170,0.25)", highlight: { background: col, border: "#000000" }, hover: { background: col, border: "#000000" } } };
    }));
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
               title: "<b>" + esc(k.kinase) + " → " + esc(k.substrate) + "</b> (kinase → substrate, GlyGen: " + esc(k.source) + ")<br>sites: " + esc(k.sites) +
                      "<br>MoTrPAC at these sites (" + ARMLAB[st.arm] + "): " + (rtxt.length ? "<br>" + rtxt.join("<br>") : "not measured") + (any ? "<br><b>a substrate site responds (red edge)</b>" : "") }; }));
    // 5. legend
    var L = "<b>" + ARMLAB[st.arm] + "</b> · " + on(st.omes).map(function (o) { return OMES[o]; }).join(", ") + " · " + on(st.tis).join(", ") + " · " + on(st.times).join(", ") + " &nbsp;|&nbsp; ";
    if (st.colour === "response") L += "node fill = mean normalised response " + grad(st.arm === "ER" ? DIFF : RESP) + " ±" + f3(vlim) + " · <b>black outline</b> = adj. p &lt; " + thr + " in a selected cell";
    else if (st.colour === "phospho") L += "node fill = MoTrPAC phosphosites responding in selection (adj. p &lt; " + thr + "): <span class='hk-sw2' style='background:#E66100'></span>up <span class='hk-sw2' style='background:#5E3C99'></span>down <span class='hk-sw2' style='background:#8C510A'></span>both <span class='hk-sw2' style='background:#FFFFFF'></span>measured, none <span class='hk-sw2' style='background:#E6E6E6'></span>not measured";
    else { var fl = D.annFields.filter(function (f) { return "ann:" + f.col === st.colour; })[0]; L += "node fill = " + esc(fl.label) + (fl.binary ? ": <span class='hk-sw2' style='background:#1B9E77'></span>yes <span class='hk-sw2' style='background:#FFFFFF'></span>no" : " " + grad(SEQ) + " 0 … " + vmax + " (log scale; white = 0; grey = no record / metabolite)"); }
    L += " &nbsp;|&nbsp; node size = strength over shown edges · edges: " + (st.arm === "ER" ? "colour = w_EE − w_RE " + grad(DIFF) + " ±" + f3(lim) : "colour = type, dashed = negative") + ", width = |w| (limit ±" + f3(lim) + ")" + (st.kinase ? " · <b>arrows</b> = kinase → substrate (red = a substrate site responds)" : "");
    legend.innerHTML = L;
    // 6. modules: significance for the selection, and the selected module's table
    var ms = q(".hk-mod"), cur = ms.value;
    ms.innerHTML = "<option value=''>none</option>" + D.mods.map(function (m) { var best = modCells(m).reduce(function (b, r) { return (!b || r[5] < b[5]) ? r : b; }, null);
      return "<option value='" + m.id + "'" + (m.id === cur ? " selected" : "") + ">" + m.id + " (" + m.n_prot + " prot, " + m.n_met + " met)" + (best && best[5] < thr ? " ★ FDR " + f3(best[5]) : "") + "</option>"; }).join("");
    var M = D.mods.filter(function (m) { return m.id === cur; })[0];
    if (M) { var rows = modCells(M);
      panel.innerHTML = "<b>" + M.id + "</b>: " + M.members.length + " nodes · CAMERA-PR (MoTrPAC run_cameraPR; competitive, within each ome) for " + ARMLAB[st.arm] + " in the selection:" +
        (rows.length ? "<table><tr><th>tissue</th><th>ome</th><th>time</th><th>members tested</th><th>direction</th><th>z</th><th>FDR</th></tr>" + rows.map(function (r) { return "<tr" + (r[5] < thr ? " class='hk-sig'" : "") + "><td>" + r[0] + "</td><td>" + OMES[r[1]] + "</td><td>" + r[2] + "</td><td>" + r[6] + "</td><td>" + r[4] + "</td><td>" + f3(r[3]) + "</td><td>" + f3(r[5]) + "</td></tr>"; }).join("") + "</table>" : " no tested cell (fewer than 5 members measured in the selected omes / tissues).");
    } else panel.innerHTML = "";
    bar.querySelectorAll(".hk-arm").forEach(function (b) { b.checked = b.value === st.arm; });
    if (st.collapsed) clusterAll();
  }
  function modCells(m) { var out = []; on(st.omes).forEach(function (o) { on(st.tis).forEach(function (t) { on(st.times).forEach(function (h) {
    var r = m.cam[t + "|" + o + "|" + st.arm + "|" + h]; if (r) out.push([t, o, h, r[0], r[1], r[2], r[3]]); }); }); }); return out; }
  // ---------- wiring ----------
  var chk = function (cls, obj) { bar.querySelectorAll(cls).forEach(function (c) { c.onchange = function () { obj[c.value] = c.checked ? 1 : 0; apply(); }; }); };
  chk(".hk-ome", st.omes); chk(".hk-tis", st.tis); chk(".hk-time", st.times);
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
     "364 nodes (304 proteins, 60 metabolites) · 764 edges: STRING >= 700, Rhea metabolite-protein links, shared/STRING-linked enzyme + same RefMet super class · layout as figures 15a / 15b")
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
