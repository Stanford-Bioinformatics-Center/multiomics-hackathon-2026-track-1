#!/usr/bin/env Rscript
# ==============================================================================
# random_walks.R — ONE RANDOM WALK TO 3 OTHER NODES ON THE JOINT NETWORK,
#                     ENDURANCE OR RESISTANCE
# ==============================================================================
#
# PURPOSE (the question this answers)
#   Starting from any protein or metabolite, where does a short walk along the
#   joint network (step 14) lead? The edges are the same in both exercise arms;
#   only the weights differ, so the same start node can lead to different
#   places after endurance (EE) and after resistance (RE).
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. Step probabilities. From node u, the walk moves to neighbour v with
#      probability
#          P(u -> v) = sig(u,v) / sum over u's neighbours k of sig(u,k),
#          sig = sigmoid(w / s)
#      i.e. the 0-1 edge values of steps 3 / 6, rescaled so each node's
#      outgoing probabilities add up to 1. Positive weights (same-direction
#      responses) are followed more often than negative ones, but the
#      sigmoid's ceiling of 1 means no single strong edge can take (almost)
#      all the probability.
#      s = median |w| over both arms, computed SEPARATELY PER EDGE TYPE
#      (protein - protein 0.042, metabolite - metabolite 0.0070,
#      metabolite - protein 0.013), because the three types' weights are on
#      different scales (step 14, known limits); one shared s would make every
#      metabolite edge look alike (sigmoid close to 0.5). For the first two
#      types this reproduces steps 3 / 6 exactly (checked).
#   2. The function random_walk(start, arm) walks from `start` (ANY node of the
#      network; a random one if none is given) to N_STEPS = 3 other nodes and
#      returns the 4 nodes in order. `arm` only chooses whose weights set the
#      step probabilities (endurance or resistance); it is not a starting
#      point. By default it never returns to a node it has already visited: at
#      each step it chooses only among the unvisited neighbours, with their
#      step probabilities rescaled to add up to 1. A walk that runs into a dead
#      end (every neighbour already visited) is thrown away and tried again
#      (from a new random start if none was given), so the result always has
#      4 nodes, except from a start node whose part of the network has fewer
#      than 4 nodes.
#      revisit = TRUE gives a plain random walk instead (it may step back,
#      e.g. A -> B -> A -> B).
#
# HOW TO RUN
#   From the command line: ONE walk (4 nodes). Arguments, all optional: start
#   node (or "random"), arm (EE = endurance, the default, or RE = resistance),
#   seed.
#       # random start node, endurance weights
#       Rscript random_walk/random_walks.R
#       # start at PPIB, endurance weights
#       Rscript random_walk/random_walks.R PPIB
#       # start at PPIB, resistance weights
#       Rscript random_walk/random_walks.R PPIB RE
#       # random start, resistance weights, reproducible
#       Rscript random_walk/random_walks.R random RE 7
#       # quote names with spaces
#       Rscript random_walk/random_walks.R "Palmitic acid"
#   From R (defines the function, runs nothing):
#       source("random_walk/random_walks.R")
#       random_walk()                            # random start, endurance
#       random_walk("PPIB", "EE")                # e.g. PPIB, P4HB, HSP90B1, ...
#       random_walk("PPIB", "RE", seed = 1)      # reproducible
#       random_walk("ATP", "EE", revisit = TRUE) # plain walk, may step back
#
# INPUTS
#   $HACK_OUT/14_joint_edges.csv
#       node_a, node_b, edge_type, w_EE, w_RE, w_diff
#   $HACK_OUT/03_weighted_edges.csv, 06_metabolite_edges.csv
#       sig_EE / sig_RE (for the check)
#
# OUTPUTS (command-line run only; in this script's folder, random_walk/, or
# $RW_OUT; CSVs are gitignored)
#   18_transition_probs.csv   one row per edge direction: from, to, edge_type,
#                             s, w_EE, w_RE, sig_EE, sig_RE, p_EE, p_RE,
#                             p_diff (= p_EE - p_RE)
#   18_walk.csv               the walk just generated, one row per node
#                             (4 rows): step (0 = start), node, p_step (the
#                             probability of the step that reached the node,
#                             among the neighbours it could choose)
#
# KNOWN LIMITS
#   The per-type scale s is a heuristic (as in step 3), and it decides how
#   much a metabolite edge counts against a protein edge when a node has both.
#   Negative weights are treated as weak links, not as a relationship of their
#   own (use |w| instead if opposite responses should count as closeness).
#   One walk is one random draw: it shows a possible path, not the typical one.
# ==============================================================================

# Load data.table quietly.
suppressMessages(library(data.table))

# Where the pipeline's tables are (override with HACK_OUT).
OUT <- Sys.getenv(
  "HACK_OUT",
  unset = path.expand("~/Desktop/output/hackathon-2026-track1/network")
)
# Other nodes each walk moves to (so a walk has N_STEPS + 1 nodes, the start
# included).
N_STEPS <- 3
# The two arms.
ARMS <- c("EE", "RE")

# ---- step probabilities ------------------------------------------------------
# The joint network's edges (step 14), each stored once.
E <- fread(file.path(OUT, "14_joint_edges.csv"))
# Sigmoid scale per edge type: the median |w| over both arms.
E[, s := median(abs(c(w_EE, w_RE))), by = edge_type]
# Both directions of every edge (the network is undirected, the walk is not).
D <- rbind(
  E[, .(from = node_a, to = node_b, edge_type, s, w_EE, w_RE)],
  E[, .(from = node_b, to = node_a, edge_type, s, w_EE, w_RE)]
)
# The 0-1 edge values, per arm.
D[, `:=`(sig_EE = plogis(w_EE / s), sig_RE = plogis(w_RE / s))]
# Step probabilities: each node's 0-1 values divided by their sum over the
# node's neighbours.
D[, `:=`(p_EE = sig_EE / sum(sig_EE), p_RE = sig_RE / sum(sig_RE)), by = from]
# The difference between the arms (positive = the step is more likely in
# endurance).
D[, p_diff := p_EE - p_RE]
# Fixed row order, indexed by the start of the step (fast lookup of a node's
# neighbours).
setkey(D, from, to)
# Every node of the joint network: any of them can be a start node.
NODES <- sort(unique(D$from))

# Safety check: the protein - protein and metabolite - metabolite 0-1 values
# equal steps 3 and 6.
pp <- fread(file.path(OUT, "03_weighted_edges.csv"))[
  , .(from = symbol_a, to = symbol_b, sig_EE, sig_RE)
]
mm <- fread(file.path(OUT, "06_metabolite_edges.csv"))[
  , .(from = metabolite_a, to = metabolite_b, sig_EE, sig_RE)
]
chk <- merge(
  rbind(pp, mm), D,
  by = c("from", "to"), suffixes = c("_step", "")
)
stopifnot(
  nrow(chk) == nrow(pp) + nrow(mm),
  isTRUE(all.equal(chk$sig_EE, chk$sig_EE_step)),
  isTRUE(all.equal(chk$sig_RE, chk$sig_RE_step))
)
# Safety check: every node's outgoing probabilities add up to 1 in both arms.
stopifnot(
  D[, .(a = abs(sum(p_EE) - 1), b = abs(sum(p_RE) - 1)), by = from][
    , max(a, b)
  ] < 1e-12
)

# ---- the walk ----------------------------------------------------------------
# random_walk(start, arm): one walk from `start` (any node; a random one if not
# given) to n_steps other nodes, using the arm's step probabilities (the arm
# only sets the weights; it is not a place to start).
# Returns a data.table with one row per node visited, n_steps + 1 = 4 rows:
# step (0 = start), node, p_step. Which arm was used is the caller's choice, so
# it is not repeated on every row.
#   revisit = FALSE (default): only unvisited neighbours can be chosen. If a
#     walk runs into a dead end (every neighbour already visited) before
#     reaching n_steps other nodes, it is thrown away and started again, from
#     a new random node if no start was given, up to max_tries times. Only a
#     start node whose part of the network has fewer than 4 nodes can then
#     give a shorter walk (with a warning).
#   revisit = TRUE: a plain random walk; any neighbour can be chosen, so it
#     never runs into a dead end.
#   seed: set it for a reproducible walk; NULL (default) gives a new walk every
#     call.
random_walk <- function(start = NULL, arm = "EE", n_steps = N_STEPS,
                        revisit = FALSE, seed = NULL, max_tries = 1000) {
  # check the arguments
  arm <- match.arg(arm, ARMS)
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(start) && !start %in% NODES) {
    stop(
      "'", start, "' is not a node of the joint network ",
      "(see 14_joint_nodes.csv)"
    )
  }
  # try until the walk reaches n_steps other nodes (a dead end ends a try
  # early); keep the longest try
  best <- NULL
  for (try in seq_len(max_tries)) {
    # the start: the given node, or any node of the network at random (a new
    # one on every try)
    from <- if (is.null(start)) sample(NODES, 1) else start
    # the walk so far: the start node, reached with probability 1
    path <- from
    p_step <- 1
    for (k in seq_len(n_steps)) {
      # the current node's neighbours and this arm's step probabilities
      nb <- D[.(path[k]), .(to, p = get(paste0("p_", arm)))]
      # without revisits: drop visited nodes, rescale the rest to add up to 1;
      # stop if none is left
      if (!revisit) nb <- nb[!to %in% path]
      if (!nrow(nb)) break
      nb[, p := p / sum(p)]
      # choose the next node (by row number: sample() on a single value would
      # draw from 1..value instead)
      i <- sample.int(nrow(nb), 1, prob = nb$p)
      path <- c(path, nb$to[i])
      p_step <- c(p_step, nb$p[i])
    }
    if (is.null(best) || length(path) > nrow(best)) {
      best <- data.table(
        step = seq_along(path) - 1L, node = path, p_step = p_step
      )
    }
    if (length(path) == n_steps + 1) break
  }
  if (nrow(best) < n_steps + 1) {
    warning(
      "no walk of ", n_steps + 1, " different nodes from '", start, "' in ",
      max_tries, " tries (its part of the network is too small); returning ",
      nrow(best), " nodes"
    )
  }
  best
}

# ---- command-line run --------------------------------------------------------
# Only when run with Rscript (not when source()d): ONE walk, written to
# 18_walk.csv (4 rows).
#   arguments, all optional, in this order: start node (or "random"),
#   arm (EE or RE), seed
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  # start node: the first argument, or a random node if none is given (or if
  # it is "random")
  start <- if (length(args) && args[1] != "random") args[1] else NULL
  # arm: whose weights set the step probabilities (default endurance)
  arm <- if (length(args) > 1) toupper(args[2]) else "EE"
  # seed: makes the walk reproducible (default: a new walk every run)
  seed <- if (length(args) > 2) as.integer(args[3]) else NULL
  # Where the results go (override with RW_OUT): the folder the script is in.
  # *.csv is gitignored.
  here <- sub(
    "^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
  )
  RW_OUT <- Sys.getenv("RW_OUT", unset = dirname(normalizePath(here)))
  dir.create(RW_OUT, recursive = TRUE, showWarnings = FALSE)
  # the walk
  walk <- random_walk(start, arm, seed = seed)
  # Safety check: no node repeats, and every step follows an edge of the
  # network.
  x <- walk$node
  stopifnot(
    !anyDuplicated(x),
    length(x) == 1 ||
      nrow(D[.(head(x, -1), x[-1]), nomatch = NULL]) == length(x) - 1
  )
  fwrite(
    D[, .(
      from, to, edge_type, s, w_EE, w_RE, sig_EE, sig_RE, p_EE, p_RE, p_diff
    )],
    file.path(RW_OUT, "18_transition_probs.csv")
  )
  fwrite(walk, file.path(RW_OUT, "18_walk.csv"))
  cat(sprintf(
    "Random walk (%s weights), start %s:\n  %s%s\n",
    c(EE = "endurance", RE = "resistance")[match.arg(arm, ARMS)],
    if (is.null(start)) "chosen at random" else "given",
    paste(x, collapse = " -> "),
    if (length(x) < N_STEPS + 1) {
      sprintf(
        "   (only %d nodes: this part of the network is too small for %d)",
        length(x), N_STEPS + 1
      )
    } else {
      ""
    }
  ))
  cat("\nWritten: 18_transition_probs.csv, 18_walk.csv in", RW_OUT, "\n")
}
