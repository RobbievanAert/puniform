#' p-uniform
#'
#' Function to apply p-uniform method for one-sample mean, two-independent means, and one raw correlation coefficient as described in van Assen, van Aert, and Wicherts (2015) and van Aert, Wicherts, and van Assen (2015). \cr
#' \cr
#' Please note that the method is still in development and that this is a beta version. If you suspect a bug, please send me an email (\email{R.C.M.vanAert@@tilburguniversity.edu}).
#'
#' @param mi A vector of group means for one-sample means
#' @param ri A vector of raw correlations
#' @param ni A vector of sample sizes for one-sample means and correlations
#' @param sdi A vector of standard deviations for one-sample means
#' @param m1i A vector of means in group 1 for two-independent means
#' @param m2i A vector of means in group 2 for two-independent means
#' @param n1i A vector of sample sizes in group 1 for two-independent means
#' @param n2i A vector of sample sizes in group 2 for two-independent means
#' @param sd1i A vector of standard deviations in group 1 for two-independent means
#' @param sd2i A vector of standard deviations in group 2 for two-independent means
#' @param tobs A vector of t-values
#' @param alpha A integer specifying the alpha level as used in primary studies (default is 0.05).
#' @param side A character indicating whether the effect sizes in the primary studies are in the right-tail of the distribution (i.e., positive)  or in the left-tail of the distribution (i.e., negative) (either \code{"right"} or \code{"left"})
#' @param method A character indicating the method to be used (\code{"P"} (default), \code{"LNP"}, \code{"LN1MINP"}, \code{"KS"}, or \code{"AD"})
#' @param plot A logical indicating whether a plot showing the relation between observed and expected p-values has to be rendered (default is \code{TRUE})
#'
#' @details Three different effect sizes can be used as input for the \code{puniform} function: one-sample means, two-independent means, and raw correlation coefficients.
#' Analyzing one-sample means and two-independent means can be done by either providing the function group means (\code{mi} or \code{m1i} and \code{m2i}), standard deviations
#' (\code{sdi} or \code{sd1i} and \code{sd2i}), and sample sizes (\code{ni} or \code{n1i} and \code{n2i}) or t-values (\code{tobs}) and sample sizes (\code{ni} or \code{n1i}
#' and \code{n2i}). Both options should be accompanied with input for the arguments \code{side} and \code{method}. See the Example section for examples. Raw correlation
#' coefficients can be analyzed by supplying \code{ri} and \code{ni} to the \code{puniform} function next to input for the arguments \code{side} and \code{method}.
#'
#' P-uniform assumes that two-tailed hypothesis tests were conducted in the primary studies. In case one-tailed hypothesis tests were conducted in the primary studies, the
#' alpha level has to be multiplied by two. For example, if one-tailed hypothesis tests were conducted with an alpha level of .05, an alpha of 0.1 has to be
#' submitted to p-uniform.
#'
#' Note that only one effect size measure can be specified at a time. A combination of effect size measures usually causes true heterogeneity among effect sizes and including
#' different effect size measures is therefore not recommended.
#'
#' Five different estimators can be used when applying p-uniform. The \code{P} method is based on the distribution of the sum of independent uniformly distributed random
#' variables (Irwin-Hall distribution) and is the recommended estimator (van Aert et al., 2015). The \code{LNP} estimator refers to Fisher’s method (1950, Chapter 4)
#' for combining p-values and the \code{LN1MINP} estimator first computes 1 – p-value in each study before applying Fisher’s method on these transformed p-values
#' (van Assen et al., 2015). \code{KS} and \code{AD} respectively use the Kolmogorov-Smirnov test (Massey, 1951) and the Anderson-Darling test (Anderson & Darling, 1954)
#' for testing whether the (conditional) p-values follow a uniform distribution.
#'
#' @return
#' \item{est}{p-uniform's effect size estimate}
#' \item{ci.lb}{lower bound of p-uniform's confidence interval}
#' \item{ci.ub}{upper bound of p-uniform's confidence interval}
#' \item{ksig}{number of significant studies}
#' \item{L.0}{test statistic of p-uniform's test of null-hypothesis of no effect (for method \code{"P"} a z-value)}
#' \item{pval.0}{one-tailed p-value of p-uniform's test of null-hypothesis of no effect}
#' \item{L.pb}{test statistic of p-uniform's publication bias test}
#' \item{pval.pb}{one-tailed p-value of p-uniform's publication bias test}
#' \item{est.fe}{effect size estimate based on traditional fixed-effect meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional fixed-effect meta-analysis}
#' \item{zval.fe}{test statistic of the null-hypothesis of no effect based on traditional fixed-effect meta-analysis}
#' \item{pval.fe}{one-tailed p-value of the null-hypothesis of no effect based on traditional fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional fixed-effect meta-analysis}
#' \item{ci.ub.fe}{ci.ub.fe upper bound of confidence interval based on traditional fixed-effect meta-analysis}
#' \item{Qstat}{test statistic of the Q-test for testing the null-hypothesis of homogeneity}
#' \item{Qpval}{one-tailed p-value of the Q-test}
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references Anderson, T. W., & Darling, D. A. (1954). A test of goodness of fit. Journal of the American Statistical Association, 49(268), 765-769.
#' @references Fisher, R. A. (1950). Statistical methods for research workers (11th ed.). London: Oliver & Boyd.
#' @references Massey, F. J. (1951). The Kolmogorov-Smirnov test for goodness of fit. Journal of the American Statistical Association, 46(253), 68-78.
#' @references Van Aert, R. C. M., Wicherts, J. M., & Van Assen, M. A. L. M. (2015). Conducting meta-analyses on p-values: Reservations and recommendations for applying p-uniform and p-curve. Manuscript submitted for publication.
#' @references Van Assen, M. A. L. M., Van Aert, R. C. M., & Wicherts, J. M. (2015). Meta-analysis using effect size distributions of only statistically significant studies. Psychological Methods, 20(3), 293-309. doi: http://dx.doi.org/10.1037/met0000025
#'
#' @examples ### Load data from meta-analysis by McCall and Carriger (1993)
#' data(data.mccall93)
#'
#' ### Apply p-uniform method to get the same results as in van Assen et al. (2015)
#' puniform(ri = data.mccall93$ri, ni = data.mccall93$ni, alpha = .025, side = "right", method = "LNP", plot = TRUE)
#'
#' ### Note that the results of the publication bias test of p-uniform are not exactly equal to the results as stated in van Assen et al. (2015).
#' ### This is caused by a minor mistake in the analyses in van Assen et al. (2015).
#'
#' ### Generate some example data
#' set.seed(123)
#' ni <- 100
#' sdi <- 1
#' mi <- rnorm(8, mean = 0.2, sd = sdi/sqrt(ni))
#' tobs <- mi/(sdi/sqrt(ni))
#'
#' ### Apply p-uniform method based on sample means
#' puniform(mi = mi, ni = ni, sdi = sdi, alpha = 0.05, side = "right", method = "P", plot = FALSE)
#'
#' ### Apply p-uniform method based on t-values
#' puniform(ni = ni, tobs = tobs, alpha = 0.05, side = "right", method = "P", plot = FALSE)
#'
#' @export

puniform <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, alpha = .05, side, method, plot = TRUE) {

  if (!missing("mi") & !missing("ni") & !missing("sdi")) {
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("ni") & !missing("tobs")) {
    measure <- "MT"
    es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") & !missing("sd1i") & !missing("sd2i")) {
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) {
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("ri") & !missing("ni")) {
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side, measure = measure)
  }

  res1 <- pubbias(es = es, alpha = alpha/2, method = method)

  if(res1$ksig == 0) {
    stop("No significant studies on the specified side")
  }

  res2 <- testeffect(zval = res1$data$zval, zcv = res1$data$zcv, ksig = res1$ksig, method = method)

  res3 <- esest(yi = res1$data$yi, vi = res1$data$vi, zval = res1$data$zval, zcv = res1$data$zcv, ksig = res1$ksig, method = method)

  if(plot == TRUE) { plottrans(tr.q = res3$tr.q, ksig = res1$ksig) }

  res5 <- transform(res1 = res1, res3 = res3, side = side, measure = measure)

  x <- list(method = method, est = res5$est, ci.lb = res5$ci.lb, ci.ub = res5$ci.ub, ksig = res1$ksig, approx.est = res3$approx.est,
            approx.ci.lb = res3$approx.ci.lb, ext.lb = res3$ext.lb, L.0 = res2$L.0, pval.0 = res2$pval.0, approx.0.imp = res2$approx.0.imp,
            L.pb = res1$L.pb, pval.pb = res1$pval.pb, approx.pb = res1$approx.pb, est.fe = res5$est.fe, se.fe = res5$se.fe, zval.fe = res5$zval.fe,
            pval.fe = res1$pval.fe, ci.lb.fe = res5$ci.lb.fe, ci.ub.fe = res5$ci.ub.fe, Qstat = res1$Qstat, Qpval = res1$Qpval)

  class(x) <- "puniformoutput"
  return(x)

}
