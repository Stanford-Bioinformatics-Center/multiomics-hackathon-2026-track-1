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
#     - three views sharing one layout: endurance (EE), resistance (RE), difference (w_EE - w_RE),
#       encoded as in figures 14a / 14b (node fill = mean normalised response, violet = down, orange = up;
#       edge colour = edge type, dashed = negative weight; difference: red = higher in EE, blue = higher in
#       RE, line style = edge type); proteins are circles, metabolites triangles (all three networks);
#     - search box, click-to-highlight neighbours, metabolite-class selector, edge-type check boxes,
#       class outlines ("bubbles", as in 15a / 15b), collapse / expand metabolite classes, hover tooltips
#       with every value (weights, responses, strengths, STRING score, Rhea reactions, link type).
#   And for each network a Cytoscape.js JSON file (.cyjs, with node positions) plus one Cytoscape style
#   file with three styles (EE, RE, difference), and plain node / edge tables.
#
# HOW TO RUN
#   After steps 10, 14 and 15:   Rscript network/17_interactive_networks.R   (about 30 seconds; needs pandoc,
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

# The legend text for each view (HTML), with CSS colour bars.
legend_html <- function(P, types) {
  bar <- function(pal) sprintf("<span class='hk-grad' style='background:linear-gradient(90deg,%s)'></span>", paste(pal, collapse = ","))
  sw <- paste(sprintf("<span class='hk-sw' style='background:%s'></span>%s", TYPE_COL[types], types), collapse = " &nbsp; ")
  arm <- function(lbl) sprintf(paste0("<b>%s</b> &nbsp;·&nbsp; node fill = mean normalised response %s −%s … +%s &nbsp;·&nbsp; node size = strength (sum |w|)",
                                      " &nbsp;·&nbsp; edge colour = type: %s &nbsp;·&nbsp; dashed = negative weight"), lbl, bar(RESP_PAL), f3(P$lim_r), f3(P$lim_r), sw)
  sty <- c("protein - protein" = "solid", "metabolite - metabolite" = "long dash", "metabolite - protein" = "dotted")[types]
  list(EE = arm("Endurance vs control"), RE = arm("Resistance vs control"),
       D = sprintf(paste0("<b>Endurance − resistance</b> &nbsp;·&nbsp; edge colour = w_EE − w_RE: higher in resistance %s higher in endurance (±%s)",
                          " &nbsp;·&nbsp; width = |difference| &nbsp;·&nbsp; line: %s &nbsp;·&nbsp; grey node size = |strength difference|"),
                   bar(DIFF_PAL), f3(P$lim_d), paste(sprintf("%s = %s", sty, types), collapse = ", ")))
}

# The controls (plain JavaScript, run once the page has drawn the network). cfg comes from R.
JS <- r"---(
function(el, x, cfg) {
  var net = document.getElementById("graph" + el.id).chart;
  var nodes = net.body.data.nodes, edges = net.body.data.edges;
  var view = "EE", focus = null, collapsed = false, showHulls = cfg.hulls.length > 0;
  var typeOn = {}; cfg.types.forEach(function (t) { typeOn[t] = true; });
  var DASH = { solid: false, neg: [5, 5], mm: [10, 5], mp: [2, 4] };
  // toolbar above the network and legend below it
  var bar = document.createElement("div"); bar.className = "hk-bar";
  var h = "<b>View</b> <button data-v='EE'>Endurance</button><button data-v='RE'>Resistance</button>" +
          "<button data-v='D'>Difference (EE − RE)</button><span class='hk-sep'></span>" +
          "<b>Find</b> <input class='hk-find' list='" + el.id + "-dl' placeholder='gene or metabolite'>" +
          "<datalist id='" + el.id + "-dl'>" + nodes.getIds().sort().map(function (i) { return "<option value=\"" + i.replace(/"/g, "&quot;") + "\">"; }).join("") + "</datalist>";
  if (cfg.classes.length) {
    h += "<span class='hk-sep'></span><b>Class</b> <select class='hk-cls'><option value=''>all</option>" +
         cfg.classes.map(function (c) { return "<option>" + c + "</option>"; }).join("") + "</select>";
  }
  if (cfg.types.length > 1) {
    h += "<span class='hk-sep'></span><b>Edges</b> " + cfg.types.map(function (t) {
      return "<label><input type='checkbox' class='hk-type' value='" + t + "' checked> " + t + "</label>"; }).join(" ");
  }
  if (cfg.hulls.length) h += "<span class='hk-sep'></span><label><input type='checkbox' class='hk-hull' checked> class outlines</label>" +
                             " <button class='hk-col'>Collapse classes</button>";
  h += "<span class='hk-sep'></span><button class='hk-reset'>Reset</button>";
  bar.innerHTML = h;
  var legend = document.createElement("div"); legend.className = "hk-legend";
  var help = document.createElement("div"); help.className = "hk-help";
  help.innerHTML = "Scroll to zoom, drag to pan (labels appear as you zoom in) · click a node to highlight it and its neighbours, click empty space to clear · hover for values" +
                   (cfg.classes.length ? " · double-click a collapsed class to open it" : "");
  el.parentNode.insertBefore(bar, el); el.parentNode.appendChild(legend); el.parentNode.appendChild(help);
  var q = function (s) { return bar.querySelector(s); };

  // collapse / expand metabolite classes (clusters placed at the centre of their members)
  function openAll() { cfg.classes.forEach(function (c) { var id = "class: " + c; if (net.isCluster(id)) net.openCluster(id); }); }
  function clusterAll() {
    cfg.classes.forEach(function (c) {
      var ids = nodes.getIds({ filter: function (n) { return n.group === c; } });
      if (ids.length < 2) return;
      var p = net.getPositions(ids), cx = 0, cy = 0;
      ids.forEach(function (i) { cx += p[i].x; cy += p[i].y; });
      net.cluster({ joinCondition: function (o) { return o.group === c; },
                    clusterNodeProperties: { id: "class: " + c, label: c + " (" + ids.length + ")", shape: "triangle", size: 20,
                      x: cx / ids.length, y: cy / ids.length, physics: false, font: { size: 18, face: "Helvetica" },
                      color: { background: "#D9D9D9", border: "#404040" }, title: "<b>" + c + "</b><br>" + ids.length + " metabolites (double-click to open)" } });
    });
  }
  // draw the current view: colours, sizes, widths, line styles, hidden edge types, then any highlight
  function apply() {
    if (collapsed) openAll();
    var keep = null;
    if (focus) { keep = {}; focus.forEach(function (i) { keep[i] = true; net.getConnectedNodes(i).forEach(function (j) { keep[j] = true; }); }); }
    nodes.update(nodes.get().map(function (n) {
      var on = !keep || keep[n.id], c = n["col" + view];
      return { id: n.id, size: n["size" + view], font: { color: on ? "#1A1A1A" : "rgba(0,0,0,0.06)" },
               color: { background: on ? c : "rgba(230,230,230,0.25)", border: on ? "#404040" : "rgba(170,170,170,0.25)",
                        highlight: { background: c, border: "#000000" }, hover: { background: c, border: "#000000" } } };
    }));
    var fset = {}; if (focus) focus.forEach(function (i) { fset[i] = true; });
    edges.update(edges.get().map(function (e) {
      var on = !focus || fset[e.from] || fset[e.to], c = e["col" + view];
      return { id: e.id, width: e["wd" + view], dashes: DASH[e["dash" + view]], hidden: !typeOn[e.etype],
               color: { color: on ? c : "rgba(200,200,200,0.12)", highlight: c, hover: c, opacity: on ? (view === "D" ? 0.95 : 0.75) : 1 } };
    }));
    legend.innerHTML = cfg.legend[view];
    bar.querySelectorAll("button[data-v]").forEach(function (b) { b.classList.toggle("hk-on", b.getAttribute("data-v") === view); });
    if (collapsed) clusterAll();
  }
  // class outlines ("bubbles"): a thick rounded stroke around each group's convex hull, drawn under the network
  function hull(pts) {
    pts.sort(function (a, b) { return a.x - b.x || a.y - b.y; });
    if (pts.length < 3) return pts;
    var cr = function (o, a, b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }, lo = [], up = [];
    pts.forEach(function (p) { while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], p) <= 0) lo.pop(); lo.push(p); });
    pts.slice().reverse().forEach(function (p) { while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], p) <= 0) up.pop(); up.push(p); });
    lo.pop(); up.pop(); return lo.concat(up);
  }
  net.on("beforeDrawing", function (ctx) {
    if (!showHulls || collapsed) return;
    cfg.hulls.forEach(function (g) {
      var ids = nodes.getIds({ filter: function (n) { return n.hull === g.key; } }); if (!ids.length) return;
      var p = net.getPositions(ids), pts = ids.map(function (i) { return { x: p[i].x, y: p[i].y }; }), hp = hull(pts);
      var path = function () { ctx.beginPath(); ctx.moveTo(hp[0].x, hp[0].y); hp.forEach(function (v) { ctx.lineTo(v.x, v.y); }); ctx.closePath(); };
      ctx.save(); ctx.lineJoin = "round"; ctx.lineCap = "round";
      path(); ctx.lineWidth = 46; ctx.strokeStyle = "rgba(90,90,90,0.55)"; ctx.stroke();
      path(); ctx.lineWidth = 43; ctx.strokeStyle = "#F3F3F3"; ctx.stroke(); ctx.fillStyle = "#F3F3F3"; ctx.fill();
      var top = Math.min.apply(null, pts.map(function (v) { return v.y; })), mx = pts.reduce(function (s, v) { return s + v.x; }, 0) / pts.length;
      ctx.font = "italic bold 15px Helvetica"; ctx.fillStyle = "#4D4D4D"; ctx.textAlign = "center"; ctx.fillText(g.label, mx, top - 28);
      ctx.restore();
    });
  });
  // wire the controls
  bar.querySelectorAll("button[data-v]").forEach(function (b) { b.onclick = function () { view = b.getAttribute("data-v"); apply(); }; });
  q(".hk-find").onchange = function () {
    var id = this.value; if (nodes.get(id) === null) return;
    if (collapsed) { collapsed = false; openAll(); var cb = q(".hk-col"); if (cb) cb.textContent = "Collapse classes"; }
    focus = [id]; apply(); net.selectNodes([id]); net.focus(id, { scale: 1.3, animation: { duration: 600 } });
  };
  if (q(".hk-cls")) q(".hk-cls").onchange = function () {
    var c = this.value; if (!c) { focus = null; apply(); net.fit({ animation: true }); return; }
    focus = nodes.getIds({ filter: function (n) { return n.group === c; } }); apply(); net.fit({ nodes: focus, animation: { duration: 600 } });
  };
  bar.querySelectorAll(".hk-type").forEach(function (cb) { cb.onchange = function () { typeOn[cb.value] = cb.checked; apply(); }; });
  if (q(".hk-hull")) q(".hk-hull").onchange = function () { showHulls = this.checked; net.redraw(); };
  if (q(".hk-col")) q(".hk-col").onclick = function () {
    if (collapsed) { collapsed = false; openAll(); this.textContent = "Collapse classes"; } else { collapsed = true; clusterAll(); this.textContent = "Expand classes"; }
    net.redraw();
  };
  q(".hk-reset").onclick = function () {
    focus = null; view = "EE"; q(".hk-find").value = ""; if (q(".hk-cls")) q(".hk-cls").value = "";
    cfg.types.forEach(function (t) { typeOn[t] = true; }); bar.querySelectorAll(".hk-type").forEach(function (cb) { cb.checked = true; });
    apply(); net.unselectAll(); net.fit({ animation: true });
  };
  net.on("click", function (p) {
    if (p.nodes.length && !net.isCluster(p.nodes[0])) { focus = [p.nodes[0]]; apply(); }
    else if (!p.nodes.length && !p.edges.length && focus) { focus = null; apply(); }
  });
  net.on("doubleClick", function (p) { if (p.nodes.length && net.isCluster(p.nodes[0])) net.openCluster(p.nodes[0]); });
  apply();
}
)---"

# Page styling for the toolbar and legend.
CSS <- tags$style(HTML("
  body { font-family: Helvetica, Arial, sans-serif; color: #1A1A1A; margin: 10px 16px; }
  .hk-bar { font-size: 13px; padding: 6px 0 8px; line-height: 2; }
  .hk-bar button { font-size: 12px; margin-right: 3px; padding: 2px 8px; border: 1px solid #999; background: #F7F7F7; border-radius: 3px; cursor: pointer; }
  .hk-bar button.hk-on { background: #1A1A1A; color: #FFF; border-color: #1A1A1A; }
  .hk-bar input.hk-find { width: 190px; font-size: 12px; } .hk-bar select { font-size: 12px; } .hk-bar label { margin-right: 6px; }
  .hk-sep { display: inline-block; width: 14px; }
  .hk-legend { font-size: 12px; padding: 6px 0 2px; } .hk-help { font-size: 11px; color: #666; }
  .hk-grad { display: inline-block; width: 90px; height: 9px; vertical-align: middle; border: 1px solid #999; margin: 0 4px; }
  .hk-sw { display: inline-block; width: 16px; height: 4px; vertical-align: middle; margin: 0 4px 0 2px; }
  div.vis-tooltip { font-family: Helvetica, Arial, sans-serif; font-size: 12px; line-height: 1.4; max-width: 420px; white-space: normal; }
"))

# Build and save one interactive page.
page <- function(P, types, classes, hulls, file, title, subtitle) {
  # vis-network node table: positions in pixels (y grows downward on screen, so it is flipped)
  vn <- P$N[, .(id = node, label = node, group = class, shape = fifelse(node_type == "metabolite", "triangle", "dot"),
                x = x * W, y = (1 - y) * H, physics = FALSE, title, hull, colEE, colRE, colD, sizeEE, sizeRE, sizeD, size = sizeEE, color = colEE)]
  # vis-network edge table
  ve <- P$E[, .(id = seq_len(.N), from = a, to = b, etype = as.character(edge_type), title, colEE, colRE, colD, wdEE, wdRE, wdD, dashEE, dashRE, dashD)]
  # settings passed to the JavaScript
  cfg <- list(types = I(types), classes = I(classes), hulls = hulls, legend = legend_html(P, types))
  # the widget: fixed positions (no physics), straight edges, labels only once zoomed in, navigation buttons
  w <- visNetwork(as.data.frame(vn), as.data.frame(ve), width = "100%", height = "780px",
                  main = list(text = title, style = "font-family:Helvetica;font-weight:bold;font-size:17px;text-align:left"),
                  submain = list(text = subtitle, style = "font-family:Helvetica;font-size:12px;color:#555;text-align:left")) |>
    visPhysics(enabled = FALSE) |> visEdges(smooth = FALSE, selectionWidth = 1.5) |>
    visNodes(borderWidth = 0.8, font = list(size = 13, face = "Helvetica"),
             scaling = list(label = list(enabled = TRUE, min = 13, max = 13, drawThreshold = 9))) |>
    visInteraction(navigationButtons = TRUE, hover = TRUE, tooltipDelay = 80, hideEdgesOnDrag = TRUE, multiselect = TRUE) |>
    onRender(JS, data = cfg)
  # the page: styles + widget, saved as one self-contained file
  w <- prependContent(w, CSS)
  saveWidget(w, file.path(HTML_DIR, file), selfcontained = TRUE, title = title)
  # the page is self-contained; remove the helper folder saveWidget leaves next to it
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
page(PJ, names(TYPE_COL), jcls, lapply(jcls, function(c) list(key = c, label = c)), "17a_joint_network.html",
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
page(PG, "protein - protein", character(0), list(), "17b_gene_network.html",
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
page(PM, "metabolite - metabolite", mcls, lapply(seq_len(nrow(mh)), function(i) list(key = mh$hull[i], label = mh$class[i])),
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
