#' print.puniformoutput
#'
#' Print method for objecs of class \code{puniformoutput}
#'
#' @param x an object of class \code{puniformoutput}
#'
#' @return The print function does not return an object.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @export

print.puniformoutput <- function(x) {
  cat("\n")
  cat("Method:", x$method)
  cat("\n")
  cat("\n")
  cat("Effect size estimation p-uniform")
  cat("\n")
  cat("\n")
  print(format(data.frame(est = round(x$est, 4), ci.lb = round(x$ci.lb, 4), ci.ub = round(x$ci.ub, 4), ksig = x$ksig, row.names = ""), width = 9))
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
  if(x$method == "LNP" | x$method == "LN1MINP") { print(format(data.frame(L.0 = round(x$L.0, 4), pval = round(x$pval.0, 4), row.names = ""), width = 9)) }
  if(x$method == "P") {
    print(format(data.frame(z.0 = round(x$L.0, 4), pval = round(x$pval.0, 4), row.names = ""), width = 9))
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
  if(x$method == "LNP" | x$method == "LN1MINP") { print(format(data.frame(L.pb = round(x$L.pb, 4), pval = round(x$pval.pb, 4), row.names = ""), width = 9)) }
  if(x$method == "P") {
    print(format(data.frame(z.pb = round(x$L.pb, 4), pval = round(x$pval.pb, 4), row.names = ""), width = 9))
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
  print(format(data.frame(est.fe = round(x$est.fe, 4), se.fe = round(x$se.fe, 4), zval.fe = round(x$zval.fe, 4),
                          pval.fe = round(x$pval.fe, 4), ci.lb.fe = round(x$ci.lb.fe, 4), ci.ub.fe = round(x$ci.ub.fe, 4),
                          Qstat = round(x$Qstat, 4), Qpval = round(x$Qpval, 4), row.names = ""), width = 9))
  cat("\n")
}
