#' p-uniform*
#'
#' Function to apply the p-uniform* method for one-sample mean, two-independent means,
#' and one raw correlation coefficient as described in van Aert and van Assen (2023).
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
#' @param mods A one-sided formula to specify the moderators to include. For 
#' example \code{x1} can be included as moderator by specifying \code{mods = ~ x1}
#' @param alpha A numerical value specifying the alpha level as used in primary studies
#' (default is 0.05 but see Details).
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
#' examples. Raw correlation coefficients can be analyzed by supplying the raw
#' correlation coefficients \code{ri} and sample sizes and \code{ni} to the 
#' \code{puni_star} function next to input for the arguments \code{side}, 
#' \code{method}, and \code{alpha}. Note that the method internally transforms the 
#' raw correlation coefficients to Fisher's z correlation coefficients. The output
#' of the function also shows the results for the Fisher's z correlation coefficient. 
#' Hence, the results need to be transformed to raw correlation coefficients if 
#' this is preferred by the user.
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
#' \bold{Selecting an estimator}
#'
#' Three different estimators are currently implemented in the \code{puni_star} function. 
#' The \code{ML} estimator refers to maximum likelihood estimation of the effect size 
#' and the between-study variance. Profile likelihood confidence intervals around 
#' the estimates are computed by means of inverting the likelihood-ratio test. 
#' Likelihood-ratio tests are used for testing the null hypotheses of no effect 
#' and no between-study variance. The \code{ML} method is the recommended method 
#' for applying p-uniform*. 
#' 
#' The two other methods (\code{P} and \code{LNP}) are moment based estimators. 
#' The method \code{P} is based on the distribution of the sum of independent 
#' uniformly distributed random variables (Irwin-Hall distribution) and the 
#' \code{LNP} method refers to Fisher's method (1950, Chapter 4). For these implementations, 
#' a p-value for testing the null hypothesis of no between-study variance can also be 
#' obtained by means of a parametric bootstrap. This is necessary since the data 
#' are otherwise first used for estimating the effect size in the procedure for testing 
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
#' \item{...}{a number of additional elements}
#'
#' @note The \code{control} argument in the \code{puni_star} function is an optional 
#' argument that gives the user more control over the optimization and root-finding 
#' algorithms. This can be especially useful if estimation of the method does not 
#' converge and NAs are returned by the function. The \code{control} argument should 
#' be specified as a list containing one or more elements. For example, 
#' \code{control = list(verbose = TRUE)} Default values are used if an element is 
#' not specified. The following elements can be specified by the user:
#' 
#' \describe{
#' \item{\code{proc.ml:}}{ A character indicating with optimization procedure should 
#' be used for the method \code{ML}. The initial implementation of p-uniform* iteratively
#' optimized the profile log-likelihood functions. As of version 0.2.6 of this 
#' package, the default optimization routine estimates both parameters at the 
#' same time. The old optimization procedure can be used by specifying \code{proc.ml = "prof"}}.
#' \item{\code{par:}}{ Starting values for the optimization procedure in case of
#' method \code{ML}. The default values are zeros.}
#' \item{\code{bounds.int:}}{ A vector of length two that is used for determining the 
#' bounds for estimating the effect size with \code{P} and \code{LNP}. The default 
#' values are a function of the \code{yi}. The lower bound is the minimum \code{yi} 
#' minus 1 and the upper bound is the maximum \code{yi} plus 1. The effect size 
#' has to be between the lower and upper bound.}
#' \item{\code{tau.int:}}{ A vector of length two that indicates the lower and upper 
#' bound of the interval that is used for estimating the between-study variance. 
#' The estimate of the between-study variance should be included in this interval. 
#' This interval is used for the methods \code{P} and \code{LNP} 
#' and its default values are (0, 2).}
#' \item{\code{est.ci:}}{ A vector of length two indicating the values that are 
#' subtracted from and added to the estimate of the effect size for computing the 
#' 95\% confidence intervals. This vector is used for the methods \code{ML}, 
#' \code{P}, and \code{LNP} and its default values are (3, 3). To give an example, 
#' estimates for the lower and upper bound around the effect size estimate are 
#' searched on the interval (est-3, est) and (est, est+3), respectively.}
#' \item{\code{tau2.ci:}}{ A vector of length two indicating the values that are 
#' added to the estimate of the between-study variance for computing the 95\% confidence 
#' intervals. This vector is used for the methods \code{ML}, \code{P}, and \code{LNP} 
#' and its default values are (0.5, 2).}
#' \item{\code{tol:}}{ A number indicating the desired accuracy of the estimates. 
#' This number is used for the methods \code{P} and \code{LNP} and its 
#' default value is 0.001.} 
#' \item{\code{maxit:}}{ An integer indicating the maximum number of iterations 
#' that is used for estimating the effect size and between-study variance. This 
#' number is used for the methods \code{P} and \code{LNP} and its default 
#' value is 300.}
#' \item{\code{verbose:}}{ A logical indicating whether information should be printed 
#' about the algorithm for estimating the effect size and between-study variance. 
#' This logical is used for the methods \code{ML}, \code{P}, and \code{LNP} and 
#' its default value is FALSE.}
#' \item{\code{reps:}}{ An integer indicating the number of bootstrap replications
#' for computing the bootstrapped p-value for the test of no between-study variance.
#' This integer is used for the methods \code{P} and \code{LNP} and its default value 
#' is 1,000.}
#' \item{\code{type:}}{ A character vector indicating whether Wald-based hypothesis
#' tests and confidence intervals are preferred (\code{type = "Wald"}) or 
#' likelihood-ratio tests and profile likelihood confidence intervals 
#' (\code{type = "profile"}). This character vector is used for method \code{ML}.
#' The default is "profile". There is also the option \code{type = "Wald/profile"} 
#' which implies that Wald tests and confidence intervals are computed for the 
#' fixed effects and likelihood-ratio tests and profile likelihood confidence 
#' intervals for the between-study variance.}
#' \item{\code{optimizer:}}{ A character indicating the optimizer that is used 
#' for method \code{ML}. The default value is "Nelder-Mead". The \code{optim} 
#' function is used for optimization, so the optimization methods implemented in 
#' the \code{optim} function can be used. See the documentation of the \code{optim} 
#' function for more information.}
#' } 
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references Fisher, R.A. (1950). Statistical methods for research workers (11th ed.).
#' London: Oliver & Boyd.
#' @references van Aert, R.C.M., & van Assen, M.A.L.M. (2023). Correcting for 
#' publication bias in a meta-analysis with the p-uniform* method. Manuscript submitted  
#' for publication. Preprint: https://osf.io/preprints/bitss/zqjr9/
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
#' puni_star(mi = mi, ni = ni, sdi = sdi, side = "right")
#'
#' ### Apply p-uniform* method using t-values
#' puni_star(tobs = tobs, ni = ni, side = "right")
#' 
#' ### Generate data in case of a continuous moderator variable 
#' set.seed(12345)
#' 
#' vi <- rep(0.04, 50) # Within-study variance
#' tau2 <- 0.04 # Between-study variance
#' xi <- rnorm(50) # Continuous moderator variable
#' yi <- rnorm(50, mean = 0.5*xi, sd = sqrt(vi+tau2))
#' 
#' ### Apply p-uniform* method with xi as moderator
#' puni_star(yi = yi, vi = vi, mods = ~ xi, sdi = sdi, side = "right")
#' @export

puni_star <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, yi, vi, 
                      mods = NULL, alpha = 0.05, side, method = "ML", 
                      boot = FALSE, control)
{
  
  ##### COMPUTE EFFECT SIZE, VARIANCE, AND Z-VALUES PER STUDY #####
  if (!missing("mi") & !missing("ni") & !missing("sdi")) 
  { # Mean unknown sigma
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side,
                    measure = measure)
  } else if (!missing("ni") & !missing("tobs")) 
  {
    measure <- "MT"
    es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
             !missing("sd1i") & !missing("sd2i")) 
  { # Mean difference unknown sigma
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                    sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) 
  { # Mean difference unknown sigma with observed t-value
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2, side = side,
                    measure = measure)
  } else if (!missing("ri") & !missing("ni")) 
  { # Correlation
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side, measure = measure)
  } else if (!missing("yi") & !missing("vi")) 
  { # User-specified standardized effect sizes
    measure <- "SPE"
    es <- escompute(yi = yi, vi = vi, alpha = alpha/2, side = side, measure = measure)
  }
  
  ### Number of fixed effects parameters to be estimated
  n_bs <- ifelse(is.null(mods), 1, ncol(model.matrix(mods, data = es)))
  
  ##############################################################################
  
  ### Default values for optimizing (ML) and root-finding procedures (P and LNP)
  con <- list(proc.ml = "", # Whether both parameters need to be estimated at the same time (default) or profile likelihoods need to be optimized
              par = rep(0, n_bs+1), # Starting values for ML estimation
              int = c(-2, 2), # Interval that is used for estimating ES with ML when profile likelihoods need to be optimized
              bounds.int = c(min(es$yi)-1,max(es$yi)+1), # Interval that is used for determining bounds for estimating ES (P, LNP)            
              tau.int = c(0, 2), # Interval that is used for estimating tau (ML when iteratively optimizing the profile likelihoods, P, LNP)
              ### Values that are added to the estimates of the ES and tau for estimating 
              # CIs. For example, for CIs around ES estimate lb is searched for on the 
              # interval c(est-3, est) and ub c(est, est+3)
              est.ci = c(3, 3),
              tau2.ci = c(0.5, 2),
              tol = 0.001, # Desired accuracy for the optimizing (ML) and root-finding procedures (P, LNP)
              maxit = 300, # Maximum number of iterations for the optimizing (ML) and root-finding procedures (P, LNP)
              verbose = FALSE, # If verbose = TRUE output is printed about estimation procedures for ES and tau (ML, P, LNP)
              reps = 1000, # Number of bootstrap replications for computing bootstrapped p-value test of heterogeneity (P, LNP)
              type = "profile", # Profile likelihood CIs are computed
              optimizer = "Nelder-Mead") # Optimizer that is used for ML estimation 
              
  ### Check if user has specified values in control and if yes replace values in con
  if (missing(control) == FALSE)
  {
    con.pos <- pmatch(names(control), names(con))
    con[con.pos] <- control[1:length(con.pos)]
  }
  
  ##############################################################################
  
  ### In the absence of moderators, fit an intercept-only model
  if (is.null(mods)) 
  { 
    mods <- ~ 1
    var_names <- ""
  } else
  { # Add data of moderators to es data frame
    es <- cbind(es, model.frame(mods))
    
    ### Extract variable names for the output
    var_names <- colnames(model.matrix(mods, data = es))
  }
  
  ##### EFFECT SIZE ESTIMATION, TESTS OF NO EFFECT, AND TEST OF NO BETWEEN-STUDY
  # VARIANCE #####
  res.es <- esest_nsig(es = es, mods = mods, n_bs = n_bs, method = method, 
                       boot = boot, con = con)
  
  # ##### PUBLICATION BIAS TEST #####
  # Commented out for now. More research is needed to develop a publication bias 
  # test for p-uniform* and to study its properties.
  # res.pub <- pubbias_nsig(yi = es$yi, vi = es$vi, ycv = es$zcv*sqrt(es$vi),
  #                         est = res.es$est, tau.est = res.es$tau.est, method = method)
  res.pub <- data.frame(L.pb = NA, pval.pb = NA)
  
  ##### MIRROR OR TRANSFORM RESULTS #####
  res.trans <- transform_nsig(res.es = res.es, side = side)
  
  ##### CREATE OUTPUT #####
  x <- list(con = con, var_names = var_names, method = method, k = length(es$yi), 
            ksig = sum(es$pval < alpha/2), est = res.trans$est, 
            ci.lb = res.trans$ci.lb, ci.ub = res.trans$ci.ub, 
            L.0 = res.es$L.0, pval.0 = res.es$pval.0, tau2 = res.es$tau2, 
            tau2.lb = res.es$tau2.lb, tau2.ub = res.es$tau2.ub, se = res.es$se,
            L.het = res.es$L.het, pval.het = res.es$pval.het, 
            L.pb = res.pub$L.pb, pval.pb = res.pub$pval.pb, 
            optim.info = res.es$optim.info)
  
  class(x) <- "puni_staroutput"
  
  return(x)
  
}