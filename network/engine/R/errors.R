# =====================================================================================================
# errors.R — classed errors for the exercise network engine (exnet)
# =====================================================================================================
# Every problem the engine detects is raised as a condition of class c(<specific>, "exnet_error", "error"),
# so a caller (or a test) can catch exactly the failure it expects, and the message always says WHAT is
# wrong and WHERE, instead of failing later with an obscure error.
# =====================================================================================================

#' Raise a classed exnet error.
#' @param message `character(1)` human-readable description of the problem.
#' @param class `character(1)` specific error class, e.g. "exnet_input_error".
#' @return never returns; signals the error.
exnet_abort <- function(message, class = "exnet_input_error") {
  # Build a condition object that carries our class names plus the standard R error classes.
  cond <- structure(class = c(class, "exnet_error", "error", "condition"), list(message = message, call = sys.call(-1)))
  # Signal it (stops execution unless the caller handles this class).
  stop(cond)
}

#' Run a checkmate check and convert a failure into a classed exnet error.
#' @param res result of a checkmate `check_*()` call: TRUE or a character message.
#' @param what `character(1)` name of the argument being checked (for the message).
#' @return invisible TRUE when the check passed.
exnet_check <- function(res, what) {
  # checkmate check_* functions return TRUE, or a string describing the failure.
  if (!isTRUE(res)) exnet_abort(sprintf("invalid '%s': %s", what, res))
  invisible(TRUE)
}
