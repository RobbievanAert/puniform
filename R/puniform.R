#' p-uniform
#'
#' Function to apply p-uniform method for one-sample mean, two-independent means,
#' and one raw correlation coefficient as described in Van Assen, Van Aert, and
#' Wicherts (2015) and Van Aert, Wicherts, and Van Assen (2015).
#' \cr
#' \cr
#' Please note that this package is still under development and that this is a beta
#' version. If you suspect a bug, please send me an email (\email{R.C.M.vanAert@@tilburguniversity.edu}).
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
#' @param yi A vector of standardized effect sizes (see Details)
#' @param vi A vector of sampling variances belonging to the standardized effect
#' sizes (\code{yi})
#' @param alpha A integer specifying the alpha level as used in primary studies
#' (default is 0.05).
#' @param side A character indicating whether the effect sizes in the primary studies
#' are in the right-tail of the distribution (i.e., positive)  or in the left-tail
#' of the distribution (i.e., negative) (either \code{"right"} or \code{"left"})
#' @param method A character indicating the method to be used (\code{"P"} (default),
#' \code{"LNP"}, \code{"LN1MINP"}, \code{"KS"}, or \code{"AD"})
#' @param plot A logical indicating whether a plot showing the relation between
#' observed and expected p-values has to be rendered (default is \code{TRUE})
#'
#' @details Three different effect size measures can be used as input for the \code{puniform}
#' function: one-sample means, two-independent means, and raw correlation coefficients.
#' Analyzing one-sample means and two-independent means can be done by either providing
#' the function group means (\code{mi} or \code{m1i} and \code{m2i}), standard deviations
#' (\code{sdi} or \code{sd1i} and \code{sd2i}), and sample sizes (\code{ni} or
#' \code{n1i} and \code{n2i}) or t-values (\code{tobs}) and sample sizes (\code{ni}
#' or \code{n1i} and \code{n2i}). Both options should be accompanied with input
#' for the arguments \code{side}, \code{method}, and \code{alpha}. See the Example section for
#' examples. Raw correlation coefficients can be analyzed by supplying \code{ri}
#' and \code{ni} to the \code{puniform} function next to input for the arguments
#' \code{side}, \code{method}, and \code{alpha}.
#'
#' It is also possible to specify the standardized effect sizes and its sampling
#' variances directly via the \code{yi} and \code{vi} arguments. However, extensive
#' knowledge about computing standardized effect sizes and its sampling variances
#' is required and specifying standardized effect sizes and sampling variances is
#' not recommended to be used if the p-values in the primary studies are not computed
#' with a z-test. In case the p-values in the primary studies were computed with,
#' for instance, a t-test, the p-values of a z-test and t-test do not exactly
#' coincide and studies may be incorrectly included in the analyses. Furthermore,
#' critical values in the primary studies cannot be transformed to critical z-values
#' if \code{yi} and \code{vi} are used as input. This yields less accurate results.
#'
#' The \code{puniform} function assumes that two-tailed hypothesis tests were conducted
#' in the primary studies. In case one-tailed hypothesis tests were conducted in the primary studies,
#' the alpha level has to be multiplied by two. For example, if one-tailed hypothesis
#' tests were conducted with an alpha level of .05, an alpha of 0.1 has to be
#' submitted to p-uniform.
#'
#' Note that only one effect size measure can be specified at a time. A combination
#' of effect size measures usually causes true heterogeneity among effect sizes and
#' including different effect size measures is therefore not recommended.
#'
#' Five different estimators can be used when applying p-uniform. The \code{P} method
#' is based on the distribution of the sum of independent uniformly distributed random
#' variables (Irwin-Hall distribution) and is the recommended estimator (Van Aert et al., 2015).
#' The \code{LNP} estimator refers to Fisher’s method (1950, Chapter 4) for combining
#' p-values and the \code{LN1MINP} estimator first computes 1 – p-value in each
#' study before applying Fisher’s method on these transformed p-values
#' (Van Assen et al., 2015). \code{KS} and \code{AD} respectively use the Kolmogorov-Smirnov
#' test (Massey, 1951) and the Anderson-Darling test (Anderson & Darling, 1954)
#' for testing whether the (conditional) p-values follow a uniform distribution.
#'
#' @return
#' \item{est}{p-uniform's effect size estimate}
#' \item{ci.lb}{lower bound of p-uniform's confidence interval}
#' \item{ci.ub}{upper bound of p-uniform's confidence interval}
#' \item{ksig}{number of significant studies}
#' \item{L.0}{test statistic of p-uniform's test of null-hypothesis of no effect
#' (for method \code{"P"} a z-value)}
#' \item{pval.0}{one-tailed p-value of p-uniform's test of null-hypothesis of no effect}
#' \item{L.pb}{test statistic of p-uniform's publication bias test}
#' \item{pval.pb}{one-tailed p-value of p-uniform's publication bias test}
#' \item{est.fe}{effect size estimate based on traditional fixed-effect meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional
#' fixed-effect meta-analysis}
#' \item{zval.fe}{test statistic of the null-hypothesis of no effect based on traditional
#' fixed-effect meta-analysis}
#' \item{pval.fe}{one-tailed p-value of the null-hypothesis of no effect based on
#' traditional fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional fixed-effect
#' meta-analysis}
#' \item{ci.ub.fe}{ci.ub.fe upper bound of confidence interval based on traditional
#' fixed-effect meta-analysis}
#' \item{Qstat}{test statistic of the Q-test for testing the null-hypothesis of homogeneity}
#' \item{Qpval}{one-tailed p-value of the Q-test}
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references Anderson, T. W., & Darling, D. A. (1954). A test of goodness of fit.
#' Journal of the American Statistical Association, 49(268), 765-769.
#' @references Fisher, R. A. (1950). Statistical methods for research workers (11th ed.).
#' London: Oliver & Boyd.
#' @references Massey, F. J. (1951). The Kolmogorov-Smirnov test for goodness of fit.
#' Journal of the American Statistical Association, 46(253), 68-78.
#' @references Van Aert, R. C. M., Wicherts, J. M., & Van Assen, M. A. L. M. (in press).
#' Conducting meta-analyses on p-values: Reservations and recommendations for applying p-uniform and p-curve. Perspectives on Psychological Science.
#' @references Van Assen, M. A. L. M., Van Aert, R. C. M., & Wicherts, J. M. (2015).
#' Meta-analysis using effect size distributions of only statistically significant studies.
#' Psychological Methods, 20(3), 293-309. doi: http://dx.doi.org/10.1037/met0000025
#'
#' @examples ### Load data from meta-analysis by McCall and Carriger (1993)
#' data(data.mccall93)
#'
#' ### Apply p-uniform method to get the same results as in Van Assen et al. (2015)
#' puniform(ri = data.mccall93$ri, ni = data.mccall93$ni, alpha = .05, side = "right", method = "LNP", plot = TRUE)
#'
#' ### Note that the results of p-uniform's publication bias test are not exactly equal
#' ### to the results as stated in Van Assen et al. (2015).
#' ### This is caused by a small mistake in the analyses of Van Assen et al. (2015).
#'
#' ### Generate some example data for one-sample means design
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

puniform <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, yi, vi,
                     alpha = 0.05, side, method, plot = FALSE) {

  ##### COMPUTE EFFECT SIZE, VARIANCE, AND Z-VALUES PER STUDY #####
  if (!missing("mi") & !missing("ni") & !missing("sdi")) { # Mean unknown sigma
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side,
                    measure = measure)
  } else if (!missing("ni") & !missing("tobs")) {
    measure <- "MT"
    es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
             !missing("sd1i") & !missing("sd2i")) { # Mean difference unknown sigma
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                    sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) { # Mean difference unknown sigma with observed t-value
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2, side = side,
                    measure = measure)
  } else if (!missing("ri") & !missing("ni")) { # Correlation
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("yi") & !missing("vi")) { # User-specified standardized effect sizes
    measure <- "SPE"
    es <- escompute(yi = yi, vi = vi, alpha = alpha/2, side = side, measure = measure)
  }

  ##### FIXED-EFFECT META-ANALYSIS #####
  res.fe <- fe.ma(yi = es$yi, vi = es$vi)

  ##### PUBLICATION BIAS TEST #####
  res1 <- pubbias(es = es, alpha = alpha/2, method = method, est.fe = res.fe$est.fe)

  if (res1$ksig == 0) { # If there are no significant studies return an error message
    stop("No significant studies on the specified side")
  }

  ##### TEST OF AN EFFECT #####
  res2 <- testeffect(zval = res1$data$zval, zcv = res1$data$zcv, ksig = res1$ksig,
                     method = method)

  ##### EFFECT SIZE ESTIMATION #####
  res3 <- esest(yi = res1$data$yi, vi = res1$data$vi, zval = res1$data$zval, zcv = res1$data$zcv,
                ksig = res1$ksig, method = method)

  ##### PLOT ILLUSTRATING RELATIONSHIP BETWEEN OBSERVED AND EXPECTED P-VALUES #####
  if (plot == TRUE) {
    plottrans(tr.q = res3$tr.q, ksig = res1$ksig)
  }

  ##### MIRROR OR TRANSFORM RESULTS #####
  res5 <- transform(res.fe = res.fe, res1 = res1, res3 = res3, side = side, measure = measure)

  ##### CREATE OUTPUT #####
  x <- list(method = method, est = res5$est, ci.lb = res5$ci.lb, ci.ub = res5$ci.ub,
            ksig = res1$ksig, approx.est = res3$approx.est, approx.ci.lb = res3$approx.ci.lb,
            ext.lb = res3$ext.lb, L.0 = res2$L.0, pval.0 = res2$pval.0,
            approx.0.imp = res2$approx.0.imp, L.pb = res1$L.pb, pval.pb = res1$pval.pb,
            approx.pb = res1$approx.pb, est.fe = res5$est.fe, se.fe = res5$se.fe,
            zval.fe = res5$zval.fe, pval.fe = res.fe$pval.fe.one, ci.lb.fe = res5$ci.lb.fe,
            ci.ub.fe = res5$ci.ub.fe, Qstat = res.fe$Qstat, Qpval = res.fe$Qpval)

  class(x) <- "puniformoutput"
  return(x)

}