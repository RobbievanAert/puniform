#' fe_ma
#'
#' Function that conducts fixed-effect meta-analysis
#'
#' @param yi A vector of standardized effect sizes
#' @param vi A vector of sampling variances belonging to the standardized effect sizes (\code{yi})
#'
#' @details The \code{fe.ma} function can be used for conducting fixed-effect meta-analysis
#' on a set of studies
#'
#' @return
#' \item{est.fe}{effect size estimate of fixed-effect meta-analysis}
#' \item{se.fe}{standard error of estimate of fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound 95\% confidence interval}
#' \item{ci.ub.fe}{upper bound 95\% confidence interval}
#' \item{zval.fe}{z-value of test of no effect}
#' \item{pval.fe}{two-tailed p-value of test of no effect}
#' \item{pval.fe.one}{one-tailed p-value of test of no effect}
#' \item{Qstat}{test statistic of the Q-test}
#' \item{Qpval}{p-value of the Q-test}
#'
#' @export

### Function for fixed-effect meta-analysis
fe_ma <- function(yi, vi) {
  wi <- 1/vi # Weight per study
  est.fe <- sum(yi*wi)/sum(wi) # FE meta-analytic estimate
  se.fe <- sqrt(1/sum(wi)) # Standard error of meta-analytic estimate
  ci.lb.fe <- est.fe-qnorm(0.975)*se.fe # Lower bound CI meta-analytical estimate
  ci.ub.fe <- est.fe+qnorm(0.975)*se.fe # Upper bound CI meta-analytical estimate
  zval.fe <- est.fe/se.fe # Z-value for test of no effect
  pval.fe.one <- pnorm(zval.fe, lower.tail = FALSE) # Compute one-sided p-value
  pval.fe <- ifelse(pval.fe.one > 0.5, (1-pval.fe.one)*2, pval.fe.one*2) # Compute two-tailed p-value
  Qstat <- sum(wi * (yi - est.fe)^2) # Q-statistic
  Qpval <- pchisq(Qstat, df = length(yi) - 1, lower.tail = FALSE) # p-value of Q-statistic

  return(data.frame(est.fe = est.fe, se.fe = se.fe, ci.lb.fe = ci.lb.fe,
                    ci.ub.fe = ci.ub.fe, zval.fe = zval.fe, pval.fe = pval.fe,
                    pval.fe.one = pval.fe.one, Qstat = Qstat, Qpval = Qpval))
}
