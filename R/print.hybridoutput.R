#' print.hybridoutput
#'
#' Print method for objecs of class \code{hybridoutput}
#'
#' @param x an object of class \code{hybridoutput}
#' @param ... Additional arguments that can be passed to the function
#'
#' @return The \code{print.hybridoutput} function does not return an object.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @export

print.hybridoutput <- function(x, ...) {
  cat("\n")
  cat("Results Hybrid method")
  cat("\n")
  cat("\n")
  x$pval.hy <- ifelse(x$pval.hy < 0.001, "  <.001", round(x$pval.hy, 4))
  print(format(data.frame(est = round(x$est.hy, 4), x = round(x$x.hy, 4), pval = x$pval.hy, 
                          ci.lb = round(x$ci.lb.hy, 4), ci.ub = round(x$ci.ub.hy, 4),
                          row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Results Hybrid^0 method")
  cat("\n")
  cat("\n")
  x$pval.hy0 <- ifelse(x$pval.hy0 < 0.001, "  <.001", round(x$pval.hy0, 4))
  print(format(data.frame(est = round(x$est.hy0, 4), x = round(x$x.hy0, 4), pval = x$pval.hy0, 
                          ci.lb = round(x$ci.lb.hy0, 4), ci.ub = round(x$ci.ub.hy0, 4), 
                          row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Results Hybrid^R method")
  cat("\n")
  cat("\n")
  x$pval.hyr <- ifelse(x$pval.hyr < 0.001, "  <.001", round(x$pval.hyr, 4))
  x$pval.o <- ifelse(x$pval.o < 0.001, "<.001", round(x$pval.o, 4))
  if (x$measure == "COR") {
    print(format(data.frame(est = round(x$est.hyr, 4), zval = round(x$stat.hyr, 4), 
                            pval = x$pval.hyr, ci.lb = round(x$ci.lb.hyr, 4), 
                            ci.ub = round(x$ci.ub.hyr, 4), row.names = ""), width = 9))
  } else {
    print(format(data.frame(est = round(x$est.hyr, 4), tval = round(x$stat.hyr, 4), 
                            pval = x$pval.hyr, ci.lb = round(x$ci.lb.hyr, 4), 
                            ci.ub = round(x$ci.ub.hyr, 4), row.names = ""), width = 9))
  }
  cat("\n")
  cat("- Two-tailed p-value original study:", x$pval.o)
  cat("\n")
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Fixed-effect meta-analysis")
  cat("\n")
  cat("\n")
  x$pval.fe <- ifelse(x$pval.fe < 0.001, "  <.001", round(x$pval.fe, 4))
  print(format(data.frame(est = round(x$est.fe, 4), se = round(x$se.fe, 4), 
                          zval = round(x$zval.fe, 4), pval = x$pval.fe, ci.lb = round(x$ci.lb.fe, 4),
                          ci.ub = round(x$ci.ub.fe, 4), row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Replication")
  cat("\n")
  cat("\n")
  x$pval.repl <- ifelse(x$pval.repl < 0.001, "  <.001", round(x$pval.repl, 4))
  if (x$measure == "COR") {
    print(format(data.frame(est = round(x$est.repl, 4), se = round(x$se.repl, 4),
                            zval = round(x$stat.repl, 4), pval = x$pval.repl, 
                            ci.lb = round(x$ci.lb.repl, 4), ci.ub = round(x$ci.ub.repl, 4), 
                            row.names = ""), width = 9))
  } else {
    print(format(data.frame(est = round(x$est.repl, 4), se = round(x$se.repl, 4), 
                            tval = round(x$stat.repl, 4), pval = x$pval.repl, 
                            ci.lb = round(x$ci.lb.repl, 4), ci.ub = round(x$ci.ub.repl, 4), 
                            row.names = ""), width = 9))
  }
  cat("\n")
} 