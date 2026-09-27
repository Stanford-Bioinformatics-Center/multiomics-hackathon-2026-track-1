# =====================================================================================================
# toy.R — a five-protein, one-metabolite example small enough to check by hand (tests, method figure)
# =====================================================================================================
# Four response dimensions (muscle RNA, muscle protein, blood RNA, blood protein, all at one time point) stand
# in for the real 16; the metabolite has two (muscle, blood) and is doubled to four for its protein edge.
# Physical edges: A-B, B-C, A-C, C-D (STRING) and M-D (Rhea). E responds almost like A but has NO physical
# link, so it never gets an edge: the hard gate at work. Hand-computed weights are in ?toy_example.
# =====================================================================================================

#' The toy example.
#' @return list(gene_EE, gene_RE, metab_EE, metab_RE (Embedding), hard_pp, hard_mp (PhysicalEdges),
#'   expected (data.frame of hand-computed weights)).
toy_example <- function() {
  dims <- c("muscle_rna_4h", "muscle_prot_4h", "blood_rna_4h", "blood_prot_4h")
  g <- function(...) { m <- rbind(...); colnames(m) <- dims; m }
  # Endurance responses (normalised log fold changes) of the five proteins.
  gEE <- g(A = c(0.6, 0.4, 0.1, 0.0), B = c(0.5, 0.5, 0.0, 0.1), C = c(-0.2, -0.1, 0.3, 0.2), D = c(0.0, 0.1, -0.4, -0.3), E = c(0.6, 0.5, 0.1, 0.1))
  # Resistance responses.
  gRE <- g(A = c(0.1, 0.0, 0.5, 0.4), B = c(-0.2, -0.1, 0.4, 0.5), C = c(0.5, 0.6, 0.1, 0.0), D = c(0.4, 0.3, 0.0, -0.1), E = c(0.1, 0.1, 0.5, 0.4))
  # The metabolite: one value per tissue (muscle, blood).
  mEE <- matrix(c(0.4, -0.6), 1, dimnames = list("M", c("muscle_metab_4h", "blood_metab_4h")))
  mRE <- matrix(c(-0.1, 0.2), 1, dimnames = list("M", c("muscle_metab_4h", "blood_metab_4h")))
  list(gene_EE = embedding(gEE, "EE"), gene_RE = embedding(gRE, "RE"),
       metab_EE = embedding(mEE, "EE", "metabolite"), metab_RE = embedding(mRE, "RE", "metabolite"),
       hard_pp = physical_edges(c("A", "B", "A", "C"), c("B", "C", "C", "D"), "STRING", "combined score >= 700"),
       hard_mp = physical_edges("M", "D", "Rhea", "D catalyses a reaction with M"),
       # Worked by hand, e.g. A-B endurance: 0.6*0.5 + 0.4*0.5 + 0.1*0 + 0*0.1 = 0.50.
       expected = data.frame(a = c("A", "B", "A", "C", "M"), b = c("B", "C", "C", "D", "D"),
                             w_EE = c(0.50, -0.13, -0.13, -0.19, 0.46), w_RE = c(0.38, -0.12, 0.10, 0.38, -0.09), stringsAsFactors = FALSE))
}
