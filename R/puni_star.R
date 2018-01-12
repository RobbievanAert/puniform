#' p-uniform*
#'
#' Function to apply the p-uniform* method for one-sample mean, two-independent means,
#' and one raw correlation coefficient as described in van Aert and van Assen (2018).
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
#' sizes (see Details)
#' @param alpha A integer specifying the alpha level as used in primary studies
#' (default is 0.05).
#' @param side A character indicating whether the effect sizes in the primary studies
#' are in the right-tail of the distribution (i.e., positive)  or in the left-tail
#' of the distribution (i.e., negative) (either \code{"right"} or \code{"left"})
#' @param method A character indicating the method to be used \code{"ML"} (default),
#' \code{"P"}, or \code{"LNP"}
#' @param boot A logical indicating whether the p-value of testing whether the 
#' between-study variance is zero for methods \code{P} and \code{LNP} should be 
#' obtained by means of a parametric bootstrap. The default value is FALSE.
#' @param control An optional list of elements that give the user more control 
#' over the optimization and root-finding algorithms (see Note) 
#'
#' @details Three different effect size measures can be used as input for the \code{puni_star}
#' function: one-sample means, two-independent means, and raw correlation coefficients.
#' Analyzing one-sample means and two-independent means can be done by either providing
#' the function group means (\code{mi} or \code{m1i} and \code{m2i}), standard deviations
#' (\code{sdi} or \code{sd1i} and \code{sd2i}), and sample sizes (\code{ni} or
#' \code{n1i} and \code{n2i}) or t-values (\code{tobs}) and sample sizes (\code{ni}
#' or \code{n1i} and \code{n2i}). Both options should be accompanied with input
#' for the arguments \code{side}, \code{method}, and \code{alpha}. See the Example section for
#' examples. Raw correlation coefficients can be analyzed by supplying \code{ri}
#' and \code{ni} to the \code{puni_star} function next to input for the arguments
#' \code{side}, \code{method}, and \code{alpha}.
#'
#' It is also possible to specify the standardized effect sizes and its sampling
#' variances directly via the \code{yi} and \code{vi} arguments. However, extensive
#' knowledge about computing standardized effect sizes and its sampling variances
#' is required and specifying standardized effect sizes and sampling variances is
#' not recommended to be used if the p-values in the primary studies are not computed
#' with a z-test. In case the p-values in the primary studies were computed with,
#' for instance, a t-test, the p-values of a z-test and t-test do not exactly
#' coincide and studies may be incorrectly included as a statistically significant or 
#' nonsignificant effect size. Furthermore, critical values in the primary studies 
#' are not transformed to critical z-values if \code{yi} and \code{vi} are used 
#' as input. This yields less accurate results.
#'
#' The \code{puni_star} function assumes that two-tailed hypothesis tests were conducted
#' in the primary studies. In case one-tailed hypothesis tests were conducted in 
#' the primary studies, the submitted \code{alpha} argument to the \code{puni_star} 
#' function has to be multiplied by two. For example, if one-tailed hypothesis tests were 
#' conducted with an alpha level of .05, an alpha of 0.1 has to be submitted to 
#' the \code{puni_star} function.
#'
#' Note that only one effect size measure can be specified at a time. A combination
#' of effect size measures usually causes true heterogeneity among effect sizes and
#' including different effect size measures is therefore not recommended.
#'
#' \bold{Selecting a method}
#'
#' Three different methods are currently implemented in the \code{puni_star} function. 
#' The \code{ML} method refers to maximum likelihood estimation of the effect size 
#' and the between-study variance. Profile likelihood confidence intervals around 
#' the estimates are computed by means of inverting the likelihood-ratio test. 
#' Likelihood-ratio tests are used for the publication bias test and testing the 
#' null hypotheses of no effect and no between-study variance. The \code{ML} method 
#' is the recommended method for applying p-uniform*. 
#' 
#' The two other methods (\code{P} and \code{LNP}) are moment based estimators. 
#' The method \code{P} is based on the distribution of the sum of independent 
#' uniformly distributed random variables (Irwin-Hall distribution) and the 
#' \code{LNP} method refers to Fisher's method (1950, Chapter 4). For these methods, 
#' a p-value for testing the null hypothesis of no between-study variance can also be 
#' obtained by means of a parametric bootstrap. This is necessary since the data 
#' is otherwise first used for estimating the effect size in the procedure for testing 
#' the null hypothesis of no between-study variance and then also used for computing 
#' a p-value. The test of no effect is not available for the methods \code{P} and \code{LNP} 
#' and the publication bias test for these methods is not yet implemented.    
#'
#' @return
#' \item{est}{p-uniform*'s effect size estimate}
#' \item{ci.lb}{lower bound of p-uniform*'s 95\% confidence interval of the effect size}
#' \item{ci.ub}{upper bound of p-uniform*'s 95\% confidence interval of the effect size}
#' \item{L.0}{test statistic of p-uniform*'s test of the null hypothesis of no effect}
#' \item{pval.0}{one-tailed p-value of p-uniform*'s test of null hypothesis of no effect}
#' \item{tau2}{p-uniform*'s estimate of the between-study variance}
#' \item{tau2.lb}{lower bound of p-uniform*'s 95\% confidence interval of the 
#' between-study variance}
#' \item{tau2.ub}{upper bound of p-uniform*'s 95\% confidence interval of the 
#' between-study variance}
#' \item{L.het}{test statistic of p-uniform*'s test of the null hypothesis of no 
#' between-study variance}
#' \item{pval.het}{one-tailed p-value of p-uniform*'s test of null hypothesis of 
#' no between-study variance}
#' \item{pval.boot}{one-tailed p-value of p-uniform*'s test of null hypothesis 
#' of no between-study variance obtained with a parametric bootstrap}
#' \item{L.pb}{test statistic of p-uniform*'s publication bias test}
#' \item{pval.pb}{one-tailed p-value of p-uniform*'s publication bias test}
#' \item{...}{a number of additional elements}
#'
#' @note The \code{control} argument in the \code{puni_star} function is an optional 
#' argument that gives the user more control over the optimzation and root-finding 
#' algorithms. This can be especially useful if estimation of the method does not 
#' converge and NAs are returned by the function. The \code{control} argument should 
#' be specified as a list containing one or more elements. For example, 
#' \code{control = list(verbose = TRUE)} Default values are used if an element is 
#' not specified. The following elements can be specified by the user:
#' 
#' \itemize{
#' \item{\code{stval.tau:}}{ An integer that is the starting value of tau for estimation.
#' This starting value is used for the methods \code{ML}, \code{P}, and \code{LNP} and 
#' its default value is 0.}
#' \item{\code{int:}}{ A vector of length two that indicates the lower and upper 
#' bound of the interval that is used for estimating the effect size. The effect 
#' size estimate should be included in this interval. This interval is used for the 
#' methods \code{ML}, \code{P}, and \code{LNP} and its default values are (-2, 2).}
#' \item{\code{tau.int:}}{ A vector of length two that indicates the lower and upper 
#' bound of the interval that is used for estimating the between-study variance. 
#' The estimate of the between-study variance should be included in this interval. 
#' This #' interval is used for the methods \code{ML}, \code{P}, and \code{LNP} 
#' and its default values are (0, 2).}
#' \item{\code{est.ci:}}{ A vector of length two indicating the values that are 
#' added to the estimate of the effect size for computing the 95\% confidence 
#' intervals. This vector is used for the methods \code{ML}, \code{P}, and \code{LNP} 
#' and its default values are (4, 3). To give an example, estimates for the lower 
#' and upper bound around the effect size estimate are searched on the interval 
#' (est-4, est) and (est, est+3), respectively.}
#' \item{\code{tau.ci:}}{ A vector of length two indicating the values that are 
#' added to the estimate of the between-study variance for computing the 95\% confidence 
#' intervals. This vector is used for the methods \code{ML}, \code{P}, and \code{LNP} 
#' and its default values are (1, 3).}
#' \item{\code{tol:}}{ A number indicating the desired accuracy of the estimates. 
#' This number is used for the methods \code{ML}, \code{P}, and \code{LNP} and its 
#' default value is 0.001.} 
#' \item{\code{max.iter:}}{ An integer indicating the maximum number of iterations 
#' that is used for estimating the effect size and between-study variance. This 
#' number is used for the methods \code{ML}, \code{P}, and \code{LNP} and its default 
#' value is 1000.}
#' \item{\code{verbose:}}{ A logical indicating whether information should be printed 
#' about the algorithm for estimating the effect size and between-study variance. 
#' This logical is used for the methods \code{ML}, \code{P}, and \code{LNP} and 
#' its default value is FALSE.}
#' \item{\code{reps:}}{ An integer indicating the number of bootstrap replications
#' for computing the bootstrapped p-value for the test of no between-study variance.
#' This integer is used for the methods \code{P} and \code{LNP} and its default value 
#' is 1000.}
#' } 
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references Fisher, R.A. (1950). Statistical methods for research workers (11th ed.).
#' London: Oliver & Boyd.
#' @references van Aert, R.C.M., & van Assen, M.A.L.M. (2016). Manuscript in preparation.
#'
#' @examples ### Generate data for one-sample mean with mu = 0.2 and tau^2 = 0.01
#' set.seed(123)
#' ni <- rep(50, 25)
#' sdi <- rep(1, 25)
#' ui <- rnorm(25, mean = 0.2, sd = 0.1)
#' mi <- rnorm(25, mean = ui, sd = sdi/sqrt(ni))
#' tobs <- mi/(sdi/sqrt(ni))
#'
#' ### Apply p-uniform* method using sample means
#' puni_star(mi = mi, ni = ni, sdi = sdi, alpha = 0.05, side = "right", method = "ML")
#'
#' ### Apply p-uniform* method using t-values
#' puni_star(tobs = tobs, ni = ni, alpha = 0.05, side = "right", method = "ML")
#'
#' @export

puni_star <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, yi, vi, 
                      alpha = 0.05, side, method = "ML", boot = FALSE, control)
{
  
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
  
  ### Default values for optimizing (ML) and root-finding procedures (P and LNP)
  con <- list(stval.tau = 0,     # Starting value of tau for estimation (ML, P, LNP)
              int = c(-2, 2),    # Interval that is used for estimating ES (ML, P, LNP)
              tau.int = c(0, 2), # Interval that is used for estimating tau (ML, P, LNP)
              ### Values that are added to the estimates of the ES and tau for estimating 
              # CIs. For example, for CIs around ES estimate lb is searched for on the 
              # interval c(est-4, est) and ub c(est, est+3)
              est.ci = c(4, 3),
              tau.ci = c(1, 3),
              tol = 0.001,       # Desired accuracy for the optimizing (ML) and root-finding procedures (P, LNP)
              max.iter = 1000,   # Maximum number of iterations for the optimizing (ML) and root-finding procedures (P, LNP)
              verbose = FALSE,   # If verbose = TRUE output is printed about estimation procedures for ES and tau (ML, P, LNP)
              reps = 1000)       # Number of bootstrap replications for computing bootstrapped p-value test of heterogeneity (P, LNP)
  
  ### Check if user has specified values in control and if yes replace values in con
  if (missing(control) == FALSE)
  {
    con.pos <- pmatch(names(control), names(con))
    con[con.pos] <- control[1:length(con.pos)]
  }
  
  ##### EFFECT SIZE ESTIMATION #####
  res.es <- esest_nsig(yi = es$yi, vi = es$vi, ycv = es$zcv*sqrt(es$vi), 
                       method = method, con = con)
  
  ##### TEST OF AN EFFECT #####
  res.null <- testeffect_nsig(yi = es$yi, vi = es$vi, est = res.es$est,
                              tau.est = res.es$tau.est, ycv = es$zcv*sqrt(es$vi),
                              method = method)
  
  ##### TEST OF NO BETWEEN-STUDY VARIANCE #####
  res.hetero <- testhetero(yi = es$yi, vi = es$vi, est = res.es$est,
                           tau.est = res.es$tau.est, ycv = es$zcv*sqrt(es$vi),
                           method = method, boot = boot, con = con)
  
  ##### PUBLICATION BIAS TEST #####
  res.pub <- pubbias_nsig(yi = es$yi, vi = es$vi, ycv = es$zcv*sqrt(es$vi),
                          est = res.es$est, tau.est = res.es$tau.est, method = method)
  
  ##### MIRROR OR TRANSFORM RESULTS #####
  res.trans <- transform_nsig(res.es = res.es, side = side)
  
  ##### CREATE OUTPUT #####
  x <- list(method = method, k = length(es$yi), ksig = sum(es$pval < alpha/2), 
            est = res.trans$est, ci.lb = res.trans$lb, ci.ub = res.trans$ub, 
            L.0 = res.null$L.0, pval.0 = res.null$pval.0, tau2 = res.es$tau.est^2, 
            tau2.lb = res.es$tau.lb^2, tau2.ub = res.es$tau.ub^2, 
            L.het = res.hetero$L.het, pval.het = res.hetero$pval.het, 
            pval.boot = res.hetero$pval.boot, L.pb = res.pub$L.pb, 
            pval.pb = res.pub$pval.pb)
  
  class(x) <- "puni_staroutput"
  
  return(x)
  
}