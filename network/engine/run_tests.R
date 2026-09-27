#!/usr/bin/env Rscript
# =====================================================================================================
# engine/run_tests.R — run the exnet engine's tests from source (no installation needed)
# =====================================================================================================
# The tests (engine/tests/testthat/) check three things:
#   test-toy.R        a 5-protein + 1-metabolite example whose weights were computed by hand, the hard gate
#                     (similar but unlinked molecules get no edge), arm-specific edges and node strength;
#   test-defensive.R  bad inputs fail early with a classed error (unknown arm, gaps in a column, infinite
#                     values, self-loops, duplicate edges, mismatched arms or dimensions, missing molecules);
#   test-pipeline.R   the engine reproduces every edge weight of the pipeline's joint network (step 14) to
#                     within 1e-12 (skipped if the pipeline has not been run).
# HOW TO RUN:  Rscript network/engine/run_tests.R      (run_all.sh step 14t; stops with an error on any failure)
# =====================================================================================================
# Folder of this script (the package root), wherever the repository is cloned.
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
# Load the package from source and run every test; stop_on_failure makes a failing test fail the pipeline.
testthat::test_local(here, reporter = "summary", stop_on_failure = TRUE)
