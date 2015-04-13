#' p-uniform
#'
#' Function to apply p-uniform method for one-sample means, two-sample means, and correlations. \cr
#' \cr
#' Please note that the method is still in development and that this is a beta version. If you suspect a bug, please send me an email (\email{R.C.M.vanAert@@tilburguniversity.edu}).
#'
#' @param mi A vector of group means for one-sample means
#' @param ri A vector of raw correlations
#' @param ni A vector of sample sizes for one-sample means and correlations
#' @param sdi A vector of standard deviations for one-sample means
#' @param m1i A vector of means in group 1 for two-sample means
#' @param m2i A vector of means in group 2 for two-sample means
#' @param n1i A vector of sample sizes in group 1 for two-sample means
#' @param n2i A vector of sample sizes in group 2 for two-sample means
#' @param sd1i A vector of standard deviations in group 1 for two-sample means
#' @param sd2i A vector of standard deviations in group 2 for two-sample means
#' @param alpha A integer specifying the alpha level as used in primary studies (default is 0.05)
#' @param side A character indicating the direction of the tested hypothesis in the primary studies (either "right" or "left")
#' @param method A character indicating the method to be used ("P" (default), "LNP", "LN1MINP", "KS", or "AD")
#' @param plot A logical indicating whether a plot showing the relation between observed and expected p-values has to be rendered (default is TRUE)
#'
#' @details Only one effect size measure can be specified at a time. A combination of effect size measures is not possible.
#'
#' @return
#' \item{est}{p-uniform's effect size estimate}
#' \item{ci.lb}{lower bound of p-uniform's confidence interval}
#' \item{ci.ub}{upper bound of p-uniform's confidence interval}
#' \item{ksig}{number of significant studies}
#' \item{L.0}{test statistic of p-uniform's test of null-hypothesis of no effect}
#' \item{pval.0}{one-tailed p-value of p-uniform's test of null-hypothesis of no effect}
#' \item{L.pb}{test statistic of p-uniform's publication bias test}
#' \item{pval.pb}{one-tailed p-value of p-uniform's publication bias test}
#' \item{est.fe}{effect size estimate based on traditional fixed-effect meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional fixed-effect meta-analysis}
#' \item{ci.ub.fe}{ci.ub.fe upper bound of confidence interval based on traditional fixed-effect meta-analysis}
#' \item{Qstat.}{test statistic of the Q-test for testing the null-hypothesis of homogeneity}
#' \item{Qpval}{one-tailed p-value of the Q-test}
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @note The p-uniform method is described in: \cr
#' \cr
#' van Assen, M. A. L. M., van Aert, R. C. M., & Wicherts, J. M. (2014).
#' Meta-analysis using effect size distributions of only statistically significant studies. Psychological Methods,
#' Advance online publication. doi: http://dx.doi.org/10.1037/met0000025
#'
#' @examples ### Load data from meta-analysis by McCall and Carriger (1993)
#' data(data.mccall93)
#'
#' ### Apply p-uniform method to get the same results as in van Assen et al. (2014)
#' puniform(ri = data.mccall93$ri, ni = data.mccall93$ni, alpha = .025, side = "right", method = "LNP", plot = TRUE)
#'
#' ### Note that the results of the publication bias test of p-uniform are not exactly equal to the results as stated in van Assen et al. (2014).
#' ### This is caused by a minor mistake in the analyses in van Assen et al. (2014).
#'
#' @export

puniform <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, alpha = .05, side, method = "P", plot = TRUE) {

  if(!missing("mi") & !missing("ni") & !missing("sdi")) {
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, side = side, measure = measure)
  } else if(!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") & !missing("sd1i") & !missing("sd2i")) {
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i, side = side, measure = measure)
  } else if(!missing("ri") & !missing("ni")) {
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, side = side, measure = measure)
  }

  res1 <- pubbias(es = es, alpha = alpha, method = method)

  if(res1$ksig == 0) {
    stop("No significant studies on the specified side")
  }

  res2 <- testeffect(zval = res1$data$zval, ksig = res1$ksig, alpha = alpha, method = method)

  res3 <- esest(yi = res1$data$yi, vi = res1$data$vi, zval = res1$data$zval, ksig = res1$ksig, alpha = alpha, method = method)

  if(plot == TRUE) { plottrans(tr.q = res3$tr.q, ksig = res1$ksig) }

  res5 <- transform(res1 = res1, res3 = res3, side = side, measure = measure)

  x <- list(method = method, est = res5$est, ci.lb = res5$ci.lb, ci.ub = res5$ci.ub, ksig = res1$ksig, approx.est = res3$approx.est,
            approx.ci.lb = res3$approx.ci.lb, L.0 = res2$L.0, pval.0 = res2$pval.0, approx.0.imp = res2$approx.0.imp, L.pb = res1$L.pb, pval.pb = res1$pval.pb,
            approx.pb = res1$approx.pb, est.fe = res5$est.fe, se.fe = res5$se.fe, zval.fe = res5$zval.fe, pval.fe = res1$pval.fe, ci.lb.fe = res5$ci.lb.fe,
            ci.ub.fe = res5$ci.ub.fe, Qstat = res1$Qstat, Qpval = res1$Qpval)

  class(x) <- "puniformoutput"
  return(x)

}
