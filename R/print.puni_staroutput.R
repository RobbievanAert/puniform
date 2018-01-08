#' print.puni_staroutput
#'
#' Print method for objecs of class \code{puni_staroutput}
#'
#' @param x an object of class \code{puni_staroutput}
#' @param ... Additional arguments that can be passed to the function
#'
#' @return The \code{print.puni_staroutput} function does not return an object.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @export

print.puni_staroutput <- function(x, ...) {
  cat("\n")
  cat("Method: ", x$method, " (k = ", x$k, "; ksig = ", x$ksig, ")", sep = "")
  cat("\n")
  cat("\n")
  cat("Estimating effect size p-uniform")
  cat("\n")
  cat("\n")
  x$pval.0 <- ifelse(x$pval.0 < 0.001, "  <.001", round(x$pval.0, 4))
  print(format(data.frame(est = round(x$est, 4), ci.lb = round(x$ci.lb, 4),
                          ci.ub = round(x$ci.ub, 4), L.0 = round(x$L.0, 4),
                          pval = x$pval.0, row.names = ""), width = 9))
  if (x$method == "P" | x$method == "LNP")
  {
    cat("\n")
    cat("Note:")
    cat("\n")
    cat("- Test of no effect is not available for method", x$method)
    cat("\n")
  }
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Estimating between-study variance p-uniform")
  cat("\n")
  cat("\n")
  x$pval.het <- ifelse(x$pval.het < 0.001, "  <.001", round(x$pval.het, 4))
  if (x$method == "ML")
  {
    print(format(data.frame(tau2 = round(x$tau2, 4), tau2.lb = round(x$tau2.lb, 4),
                            tau2.ub = round(x$tau2.ub, 4), L.het = round(x$L.het, 4),
                            pval = x$pval.het, row.names = ""), width = 9))
  } else if (x$method == "P" | x$method == "LNP")
  {
    x$pval.boot <- ifelse(x$pval.boot < 0.001, "  <.001", round(x$pval.boot, 4))
    print(format(data.frame(tau2 = round(x$tau2, 4), tau2.lb = round(x$tau2.lb, 4),
                            tau2.ub = round(x$tau2.ub, 4), L.het = round(x$L.het, 4),
                            pval = x$pval.het, pval.boot = x$pval.boot, 
                            row.names = ""), width = 9))
  }
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Publication bias test p-uniform")
  cat("\n")
  cat("\n")
  x$pval.pb <- ifelse(x$pval.pb < 0.001, "  <.001", round(x$pval.pb, 4))
  print(format(data.frame(L.pb = round(x$L.pb, 4), pval = x$pval.pb,
                          row.names = ""), width = 9))
  if (x$method == "P" | x$method == "LNP")
  {
    cat("\n")
    cat("Note:")
    cat("\n")
    cat("- Publication bias test is not yet implemented for method", x$method)
    cat("\n")
  }
}