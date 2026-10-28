#' @export

print.puniformoutput <- function(x, ...) {
  cat("\n")
  cat("Method:", x$method)
  cat("\n")
  cat("\n")
  cat("Effect size estimation p-uniform")
  cat("\n")
  cat("\n")
  if (x$method == "LNP" | x$method == "LN1MINP" | x$method == "P" | x$method == "ML") {
    x$pval.0 <- ifelse(x$pval.0 < 0.001, "  <.001", round(x$pval.0,
                                                          4))
    print(format(data.frame(est = round(x$est, 4), ci.lb = round(x$ci.lb, 4),
                            ci.ub = round(x$ci.ub, 4), L.0 = round(x$L.0, 4),
                            pval = x$pval.0, ksig = x$ksig, row.names = ""), width = 9))
    cat("\n")
  }
  if (x$method == "KS" | x$method == "AD") {
    print(format(data.frame(est = round(x$est, 4), ci.lb = round(x$ci.lb, 4),
                            ci.ub = round(x$ci.ub, 4), L.0 = NA, pval = NA, ksig = x$ksig,
                            row.names = ""), width = 9))
    cat("\n")
    cat("Notes:")
    cat("\n")
    cat("- Confidence interval does not exist for method", x$method)
    cat("\n")
    cat("- Test of an effect does not exist for method", x$method)
    if (x$approx.est == 1) {
      cat("\n")
      cat("- Approximation used for estimating effect size")
    }
    cat("\n")
    cat("\n")
  } else if (x$approx.est == 1 & x$approx.ci.lb == 1 & is.na(x$est) & is.na(x$ci.lb)) {
    cat("Notes:")
    cat("\n")
    cat("- Approximation used for estimating effect size and ci.lb, est and ci.lb <",
        x$ext.lb)
    if (x$method == "P") {
      cat("\n")
      cat("- p-value approximated with normal distribution")
    }
    cat("\n")
    cat("\n")
  } else if (x$approx.ci.lb == 1 & is.na(x$ci.lb)) {
    cat("Notes:")
    cat("\n")
    cat("- Approximation used for estimating ci.lb, ci.lb <", x$ext.lb)
    if (x$method == "P") {
      cat("\n")
      cat("- p-value approximated with normal distribution")
    }
    cat("\n")
    cat("\n")
  } else if (x$approx.est == 1 & x$approx.ci.lb == 1 & is.na(x$est) == FALSE & is.na(x$ci.lb) == FALSE) {
    cat("Notes:")
    cat("\n")
    cat("- Approximation used for estimating effect size and ci.lb")
    if (x$method == "P") {
      cat("\n")
      cat("- p-value approximated with normal distribution")
    }
    cat("\n")
    cat("\n")
  } else if (x$approx.ci.lb == 1 & is.na(x$ci.lb) == FALSE) {
    cat("Notes:")
    cat("\n")
    cat("- Approximation used for estimating ci.lb")
    if (x$method == "P") {
      cat("\n")
      cat("- p-value approximated with normal distribution")
    }
    cat("\n")
    cat("\n")
  }
  cat("===")
  cat("\n")
  cat("\n")
  cat("Publication bias test p-uniform")
  cat("\n")
  cat("\n")
  if (x$method == "LNP" | x$method == "LN1MINP" | x$method == "P" | x$method == "ML") {
    x$pval.pb <- ifelse(x$pval.pb < 0.001, "  <.001", round(x$pval.pb,
                                                            4))
    print(format(data.frame(L.pb = round(x$L.pb, 4), pval = x$pval.pb,
                            row.names = ""), width = 9))
  } else if (x$method == "P") {
    cat("\n")
    cat("p-value approximated with normal distribution")
    cat("\n")
  } else if (x$method == "KS" | x$method == "AD") {
    cat("Publication bias test does not exist for method", x$method)
    cat("\n")
  }
  cat("\n")
  if (x$approx.pb == 1) {
    cat("P(Z>=z) and P(Z>=zcv) were approximated")
  }
  cat("===")
  cat("\n")
  cat("\n")
  cat("Fixed-effect meta-analysis")
  cat("\n")
  cat("\n")
  x$Qpval <- ifelse(x$Qpval < 0.001, "  <.001", round(x$Qpval, 4))
  x$pval.fe <- ifelse(x$pval.fe < 0.001, "  <.001", round(x$pval.fe,
                                                          4))
  print(format(data.frame(est.fe = round(x$est.fe, 4), se.fe = round(x$se.fe, 4),
                          zval.fe = round(x$zval.fe, 4), pval.fe = x$pval.fe,
                          ci.lb.fe = round(x$ci.lb.fe, 4), ci.ub.fe = round(x$ci.ub.fe, 4),
                          Qstat = round(x$Qstat, 4), Qpval = x$Qpval, row.names = ""), width = 9))
  cat("\n")
}
