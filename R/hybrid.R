#' hybrid
#'
#' Function to statistically combine an original study and replication by means
#' of the hybrid methods and fixed-effect meta-analysis as described in van Aert
#' and van Assen (2017).
#' \cr
#' \cr
#' Please note that this package is still under development and that this is a
#' beta version. If you suspect a bug, please send me an email
#' (\email{R.C.M.vanAert@@tilburguniversity.edu}).
#'
#' @param mi A vector of group means for one-sample means
#' @param ri A vector of raw correlations
#' @param ni A vector of sample sizes for one-sample means and correlations
#' @param sdi A vector of standard deviations for one-sample means
#' @param m1i A vector of means in group 1 for two-independent means
#' @param m2i A vector of means in group 2 for two-independent means
#' @param n1i A vector of sample sizes in group 1 for two-independent means
#' @param n2i A vector of sample sizes in group 2 for two-independent means
#' @param sd1i A vector of standard deviations in group 1 for two-independent
#' means
#' @param sd2i A vector of standard deviations in group 2 for two-independent
#' means
#' @param tobs A vector of t-values
#' @param alpha A integer specifying the alpha level as used in the original
#' study (default is 0.05).
#' @param side A character indicating whether the observed effect size of the
#' original study is in the right-tail of the distribution (i.e., positive) or
#' in the left-tail of the distribution (i.e., negative) (either \code{"right"}
#' or \code{"left"})
#'
#' @details Three different effect sizes can be used as input for the
#' \code{hybrid} function: one-sample means, two-independent means, and raw
#' correlation coefficients. Analyzing one-sample means and two-independent
#' means can be done by either providing the function group means (\code{mi} or
#' \code{m1i} and \code{m2i}), standard deviations (\code{sdi} or \code{sd1i}
#' and \code{sd2i}), and sample sizes (\code{ni} or \code{n1i} and \code{n2i})
#' or t-values (\code{tobs}) and sample sizes (\code{ni} or \code{n1i}
#' and \code{n2i}). Both options should be accompanied with input for the
#' arguments \code{side} and \code{method}. See the Example section for an
#' example. Raw correlation coefficients can be analyzed by supplying \code{ri}
#'  and \code{ni} to the \code{puniform} function next to input for the
#'arguments \code{side} and \code{method}. The vectors containing data of the
#'original study and replication should always be of length two; the first
#'element should contain information about the original study and the second
#'element should contain information about the replication.
#'
#' The hybrid methods assume that the original study is statistically
#' significant and a two-tailed hypothesis test was conducted in the original
#' study. In case a one-tailed hypothesis tests was conducted in the original
#' study, the alpha level has to be multiplied by two. For example, if a
#' one-tailed hypothesis test was conducted with an alpha level of .05, an alpha
#'  of 0.1 has to be entered into the \code{hybrid} function.
#'
#' @return
#' \item{est.hy}{effect size estimate of hybrid method}
#' \item{ci.lb.hy}{lower bound of hybrid method's confidence interval}
#' \item{ci.ub.hy}{upper bound of hybrid method's confidence interval}
#' \item{x.hy}{test statistic of hybrid method's test of null-hypothesis of no
#' effect}
#' \item{pval.hy}{two-tailed p-value of hybrid method's test of null-hypothesis
#'  of no effect}
#' \item{measure}{effect size measure}
#' \item{est.hyr}{effect size estimate of hybridR method}
#' \item{ci.lb.hyr}{lower bound of hybridR method's confidence interval}
#' \item{ci.ub.hyr}{upper bound of hybridR method's confidence interval}
#' \item{stat.hyr}{test statistic of hybridR method's test of null-hypothesis of
#'  no effect}
#' \item{pval.hyr}{two-tailed p-value of hybridR method's test of
#' null-hypothesis of no effect}
#' \item{pval.o}{two-tailed p-value of original study}
#' \item{est.hy0}{effect size estimate of hybrid0 method}
#' \item{ci.lb.hy0}{lower bound of hybrid0 method's confidence interval}
#' \item{ci.ub.hy0}{upper bound of hybrid0 method's confidence interval}
#' \item{x.hy0}{test statistic of hybrid0 method's test of null-hypothesis of no
#' effect}
#' \item{pval.hy0}{two-tailed p-value of hybrid0 method's test of
#' null-hypothesis of no effect}
#' \item{est.fe}{effect size estimate based on traditional fixed-effect
#' meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional
#' fixed-effect meta-analysis}
#' \item{zval.fe}{test statistic of the null-hypothesis of no effect based on
#' traditional fixed-effect meta-analysis}
#' \item{pval.fe}{two-tailed p-value of the null-hypothesis of no effect based
#' on traditional fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional
#' fixed-effect meta-analysis}
#' \item{ci.ub.fe}{upper bound of confidence interval based on
#' traditional fixed-effect meta-analysis}
#' \item{est.repl}{effect size estimate of replication}
#' \item{se.repl}{standard error of replication's effect size estimate}
#' \item{ci.lb.repl}{lower bound of replication's confidence interval}
#' \item{ci.ub.repl}{upper bound of replication's confidence interval}
#' \item{stat.repl}{test statistic of replication for testing null-hypothesis of
#'  no effect}
#' \item{pval.repl}{two-tailed p-value of replication for testing
#' null-hypothesis of no effect}
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R. C. M., & van Assen, M. A. L. M. (2017). Examining
#' reproducibility in psychology: A hybrid method for statistically combining a
#' biased original study and replication. Behavior Research Methods. doi:10.3758/s13428-017-0967-6
#'
#' @examples
#' ### Apply hybrid function to example on page 5 of van Aert and van Assen (2017).
#'
#' hybrid(tobs = c(2.211,1.04), n1i = c(40,80), n2i = c(40,80), alpha = .05, side = "right")
#'
#' @export

hybrid <- function(m1i, m2i, mi, ri, sd1i, sd2i, sdi, n1i, n2i, ni, tobs,
                   alpha, side) {

  if (!missing("mi") & !missing("ni") & !missing("sdi")) {
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side,
                    measure = measure)
    res.repl <- repl(es = es, mi = mi, sdi = sdi, ni = ni, measure = measure,
                     side = side)
  } else if (!missing("ni") & !missing("tobs")) {
    measure <- "MT"
    es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side,
                    measure = measure)
    res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
  } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
             !missing("sd1i") & !missing("sd2i")) {
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                    sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i,
                     sd1i = sd1i, sd2i = sd2i, measure = measure, side = side)
  } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) {
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2,
                    side = side, measure = measure)
    res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
  } else if (!missing("ri") & !missing("ni")) {
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side,
                    measure = measure)
    res.repl <- repl(es = es, measure = measure, side = side)
  }

  if (es$pval[1] > alpha/2) {
    stop("Original study is not statistically significant")
  }

  res1 <- hy(es = es, measure = measure, side = side, alpha = alpha/2)

  if (es$pval[1] < alpha/4) {
    res2 <- data.frame(est.hyr = res1$est.hy, ci.lb.hyr = res1$ci.lb.hy,
                       ci.ub.hyr = res1$ci.ub.hy, stat.hyr = res1$x.hy, pval.hyr = res1$pval.hy,
                       pval.o = res.repl$pval.o)
  } else {
    res2 <- data.frame(est.hyr = res.repl$est.repl, ci.lb.hyr = res.repl$ci.lb.repl,
                       ci.ub.hyr = res.repl$ci.ub.repl, stat.hyr = res.repl$stat.repl,
                       pval.hyr = res.repl$pval.repl, pval.o = res.repl$pval.o)
  }

  res3 <- hy0(es = es, res1 = res1, alpha = alpha/2)

  res4 <- fe_ma(yi = es$yi, vi = es$vi)

  ### Transform results of fixed-effect meta-analysis
  if (measure == "COR") { # Back transform Fisher z to correlation
    est.fe <- (exp(2*res4$est.fe) - 1)/(exp(2*res4$est.fe) + 1)
    se.fe <- (exp(2*res4$se.fe) - 1)/(exp(2*res4$se.fe) + 1)
    ci.lb.fe <- (exp(2*res4$ci.lb.fe) - 1)/(exp(2*res4$ci.lb.fe) + 1)
    ci.ub.fe <- (exp(2*res4$ci.ub.fe) - 1)/(exp(2*res4$ci.ub.fe) + 1)
    zval.fe <- res4$zval.fe

    if (side == "left") { # Re-mirror estimates
      est.fe <- est.fe*-1
      tmp <- ci.ub.fe
      ci.ub.fe <- ci.lb.fe*-1
      ci.lb.fe <- tmp*-1
      zval.fe <- zval.fe*-1
    }

  } else if (side == "left" & measure != "COR") { # Re-mirror estimates
    est.fe <- res4$est.fe*-1
    tmp <- res4$ci.ub.fe
    ci.ub.fe <- res4$ci.lb.fe*-1
    ci.lb.fe <- tmp*-1
    zval.fe <- res4$zval.fe*-1
  } else {
    est.fe <- res4$est.fe
    ci.ub.fe <- res4$ci.ub.fe
    ci.lb.fe <- res4$ci.lb.fe
    zval.fe <- res4$zval.fe
  }

  x <- list(est.hy = res1$est, ci.lb.hy = res1$ci.lb, ci.ub.hy = res1$ci.ub,
            x.hy = res1$x, pval.hy = res1$pval, measure = measure, est.hyr = res2$est.hyr,
            ci.lb.hyr = res2$ci.lb.hyr, ci.ub.hyr = res2$ci.ub.hyr, stat.hyr = res2$stat.hyr,
            pval.hyr = res2$pval.hyr, pval.o = res2$pval.o, est.hy0 = res3$est.hy0,
            ci.lb.hy0 = res3$ci.lb.hy0, ci.ub.hy0 = res3$ci.ub.hy0, x.hy0 = res3$x.hy0,
            pval.hy0 = res3$pval.hy0, est.fe = est.fe, se.fe = res4$se.fe,
            ci.lb.fe = ci.lb.fe, ci.ub.fe = ci.ub.fe, zval.fe = zval.fe,
            pval.fe = res4$pval.fe, est.repl = res.repl$est.repl, se.repl = res.repl$se.repl,
            ci.lb.repl = res.repl$ci.lb.repl, ci.ub.repl = res.repl$ci.ub.repl,
            stat.repl = res.repl$stat.repl, pval.repl = res.repl$pval.repl)
  class(x) <- "hybridoutput"
  return(x)
}
