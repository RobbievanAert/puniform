#' print.hybridoutput
#'
#' @export

print.hybridoutput <- function(w) {
  cat("\n")
  cat("Results Hybrid method")
  cat("\n")
  cat("\n")
  w$pval.hy <- ifelse(w$pval.hy < .001, "  <.001", round(w$pval.hy, 4))
  print(format(data.frame(est = round(w$est.hy, 4), x = round(w$x.hy, 4), pval = w$pval.hy, ci.lb = round(w$ci.lb.hy, 4), ci.ub = round(w$ci.ub.hy, 4), row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Results Hybrid^0 method")
  cat("\n")
  cat("\n")
  w$pval.hy0 <- ifelse(w$pval.hy0 < .001, "  <.001", round(w$pval.hy0, 4))
  print(format(data.frame(est = round(w$est.hy0, 4), x = round(w$x.hy0, 4), pval = w$pval.hy0, ci.lb = round(w$ci.lb.hy0, 4), ci.ub = round(w$ci.ub.hy0, 4), row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Results Hybrid^R method")
  cat("\n")
  cat("\n")
  w$pval.hyr <- ifelse(w$pval.hyr < .001, "  <.001", round(w$pval.hyr, 4))
  w$pval.o <- ifelse(w$pval.o < .001, "<.001", round(w$pval.o, 4))
  if (w$measure == "COR") {
    print(format(data.frame(est = round(w$est.hyr, 4), zval = round(w$stat.hyr, 4), pval = w$pval.hyr, ci.lb = round(w$ci.lb.hyr, 4), ci.ub = round(w$ci.ub.hyr, 4), row.names = ""), width = 9))
  } else { print(format(data.frame(est = round(w$est.hyr, 4), tval = round(w$stat.hyr, 4), pval = w$pval.hyr, ci.lb = round(w$ci.lb.hyr, 4), ci.ub = round(w$ci.ub.hyr, 4), row.names = ""), width = 9)) }
  cat("\n")
  cat("- Two-tailed p-value original study:", w$pval.o)
  cat("\n")
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Fixed-effect meta-analysis")
  cat("\n")
  cat("\n")
  w$pval.fe <- ifelse(w$pval.fe < .001, "  <.001", round(w$pval.fe, 4))
  print(format(data.frame(est = round(w$est.fe, 4), se = round(w$se.fe, 4), zval = round(w$zval.fe, 4),
                          pval = w$pval.fe, ci.lb = round(w$ci.lb.fe, 4), ci.ub = round(w$ci.ub.fe, 4),
                          row.names = ""), width = 9))
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  cat("Replication")
  cat("\n")
  cat("\n")
  w$pval.repl <- ifelse(w$pval.repl < .001, "  <.001", round(w$pval.repl, 4))
  if (w$measure == "COR") {
    print(format(data.frame(est = round(w$est.repl, 4), se = round(w$se.repl, 4), zval = round(w$stat.repl, 4),
                            pval = w$pval.repl, ci.lb = round(w$ci.lb.repl, 4), ci.ub = round(w$ci.ub.repl, 4),
                            row.names = ""), width = 9))
  } else { print(format(data.frame(est = round(w$est.repl, 4), se = round(w$se.repl, 4), tval = round(w$stat.repl, 4),
                                   pval = w$pval.repl, ci.lb = round(w$ci.lb.repl, 4), ci.ub = round(w$ci.ub.repl, 4),
                                   row.names = ""), width = 9)) }
  cat("\n")
}
