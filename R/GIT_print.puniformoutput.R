

print.puniformoutput <- function(x) {
  cat("\n")
  cat("Method:", x$method)
  cat("\n")
  cat("\n")
  cat("Effect size estimation p-uniform")
  cat("\n")
  cat("\n")
  print(format(data.frame(estimate = x$est, ci.lb = x$ci.lb, ci.ub = x$ci.ub, ksig = x$ksig, row.names = ""), width = 9))
  cat("\n")
  if(x$approx.est == 1 & x$approx.ci.lb == 1) {
    cat("Approximation used for estimating effect size and ci.lb")
    cat("\n")
    cat("\n")
  } else if(x$approx.ci.lb == 1) {
    cat("Approximation used for estimating ci.lb")
    cat("\n")
    cat("\n")
  }
  cat("===")
  cat("\n")
  cat("\n")
  cat("Test of an effect p-uniform")
  cat("\n")
  cat("\n")
  if(x$method == "LNP" | x$method == "LN1MINP") { print(format(data.frame(L.0 = x$L.0, pval = x$pval.0, row.names = ""), width = 9)) }
  if(x$method == "P") {
    print(format(data.frame(z.0 = x$L.0, pval = x$pval.0, row.names = ""), width = 9))
    cat("\n")
    cat("p-value approximated with normal distribution")
    cat("\n")
  }
  if(x$method == "KS" | x$method == "AD") {
    cat("Test of an effect does not exist for methods KS and AD")
    cat("\n")
  }
  cat("\n")
  if(x$approx.0.imp == 1) {
    cat("Imputation of minimum for calculating transformed p-value")
    cat("\n")
    cat("\n")
  }
  cat("===")
  cat("\n")
  cat("\n")
  cat("Publication bias test p-uniform")
  cat("\n")
  cat("\n")
  if(x$method == "LNP" | x$method == "LN1MINP") { print(format(data.frame(L.pb = x$L.pb, pval = x$pval.pb, row.names = ""), width = 9)) }
  if(x$method == "P") {
    print(format(data.frame(z.pb = x$L.pb, pval = x$pval.pb, row.names = ""), width = 9))
    cat("\n")
    cat("p-value approximated with normal distribution")
    cat("\n")
  }
  if(x$method == "KS" | x$method == "AD") {
    cat("Publication bias test does not exist for methods KS and AD")
    cat("\n")
  }
  cat("\n")
  if(x$approx.pb == 1) {
    cat("P(Z>=z) and P(Z>=zcv) were approximated")
  }
  cat("===")
  cat("\n")
  cat("\n")
  cat("Fixed-effect meta-analysis")
  cat("\n")
  cat("\n")
  print(format(data.frame(estimate = x$est.fe, se = x$se.fe, ci.lb = x$ci.lb.fe, ci.ub = x$ci.ub.fe, Qstat. = x$Qstat., Qpval = x$Qpval, row.names = ""), width = 9))
  cat("\n")
}
