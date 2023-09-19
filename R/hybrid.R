#' hybrid
#'
#' Function to statistically combine conventional and preregistered/replications 
#' studies by means of the hybrid methods as described in van Aert and van Assen 
#' (2018) and van Aert (2023).
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
#' @param yi A vector of standardized effect sizes (see Details)
#' @param vi A vector of sampling variances belonging to the standardized effect
#' sizes (see Details)
#' @param conventional A vector indicating whether a study is a conventional study 
#' (indicated with a 1) and therefore was susceptible to bias. Studies not susceptible
#' to bias are indicated with a 0
#' @param side A character indicating whether the observed effect size of the
#' conventional studies are in the right-tail of the distribution (i.e., positive) or
#' in the left-tail of the distribution (i.e., negative) (either \code{"right"}
#' or \code{"left"})
#' @param alpha A numerical value specifying the alpha level as used in the conventional
#' study (default is 0.05, see Details)
#' @param mods A one-sided formula to specify the moderators to include. For 
#' example \code{x1} can be included as moderator by specifying \code{mods = ~ x1}
#' @param control An optional list of elements that give the user more control 
#' over the estimation procedures (see Note) 
#'
#' @details Three different effect sizes can be used as input for the
#' \code{hybrid} function: one-sample means, two-independent means, and raw
#' correlation coefficients. Analyzing one-sample means and two-independent
#' means can be done by either providing the function group means (\code{mi} or
#' \code{m1i} and \code{m2i}), standard deviations (\code{sdi} or \code{sd1i}
#' and \code{sd2i}), and sample sizes (\code{ni} or \code{n1i} and \code{n2i})
#' or t-values (\code{tobs}) and sample sizes (\code{ni} or \code{n1i}
#' and \code{n2i}). Pearson correlation coefficients can be analyzed by supplying 
#' \code{ri} and \code{ni} to the \code{hybrid} function. These correlation 
#' coefficients are internally transformed to Fisher's z correlations before 
#' analyzing the data. The results in the output are of the Fisher's z transformed
#' correlations. It is also possible to specify the standardized effect sizes and 
#' its sampling variances directly via the \code{yi} and \code{vi} arguments. 
#'  
#' Two other arguments that need to be specified are \code{side} and \code{conventional}. 
#' \code{side} indicates whether the effect size in the conventional study was expected
#' to be in the right-tail (\code{side = "right"} or in the left-tail 
#' (\code{side = "left"}) of the distribution. The argument \code{conventional} has 
#' to be used to indicate which studies are conventional studies and expected to be 
#' susceptible to bias. A 1 indicates that a study was susceptible to bias and 
#' a 0 indicates that a study was not susceptible to bias.
#' 
#' It is assumed that two-tailed hypothesis tests were conducted in the conventional 
#' studies. In case one-tailed hypothesis tests were conducted in the conventional 
#' studies, the alpha level has to be multiplied by two. For example, if one-tailed 
#' hypothesis tests were conducted with an alpha level of .05, an alpha of 0.1 
#' has to be supplied to the \code{hybrid} function.
#' 
#' \strong{Previous version}
#' 
#' The usage of a previous version of the \code{hybrid} function was more restricted.
#' Users could only apply the method to a single conventional study and replication. 
#' Before the addition of the extra functionality to also analyze multiple conventional 
#' studies and replications, data of the conventional study and replication were 
#' specified in vectors containing two elements with the first element being 
#' the data of the conventional study and the second one the data of the replication. 
#' In order to maintain backwards compatibility, it is still possible to analyze 
#' data like this by using the arguments \code{m1i, m2i, mi, ri, sd1i, sd2i, sdi, 
#' n1i, n2i, ni, tobs}. However, using the \code{hybrid} function in this way is 
#' now deprecated.
#'
#' @return
#' \item{k}{total number of effect sizes}
#' \item{k.conventional}{number of effect sizes of conventional studies}
#' \item{est}{parameter estimates of the fixed effects of the hybrid method}
#' \item{tau2}{estimate of the between-study variance in true effect size of hybrid 
#' method}
#' \item{se}{standard error of the fixed effects and of the estimated between-study
#' variance}
#' \item{ci.lb}{lower bound of hybrid method's confidence interval of the average
#' effect size}
#' \item{ci.ub}{upper bound of hybrid method's confidence interval of the average
#' effect size}
#' \item{L.0}{test statistic of hybrid method's test of null-hypothesis of no
#' effect. This is either a z-value (\code{zval} in the output) or a 
#' chi-square value (\code{LR} in the output) in case of a likelihood-ratio test}
#' \item{pval.0}{p-value of hybrid method's test of null-hypothesis of no effect}
#'  \item{tau2}{estimate of the between-study variance in true effect size of hybrid 
#' method}
#' \item{tau2.lb}{lower bound of hybrid method's confidence interval of the
#' between-study variance}
#' \item{tau2.ub}{upper bound of hybrid method's confidence interval of the
#' between-study variance}
#' \item{L.het}{test statistic of hybrid method's test of null-hypothesis of no
#' heterogeneity This is either a z-value (\code{zval} in the output) or a 
#' chi-square value (\code{LR} in the output) in case of a likelihood-ratio test}
#' \item{pval.het}{p-value of hybrid method's test of null-hypothesis of no
#' heterogeneity}
#' \item{optim.info}{model fitting results if the implementation of van Aert 
#' (2023) is used}
#' 
#' The elements below are only returned if the deprecated implementation 
#' of van Aert and van Assen (2018) is used:
#' 
#' \item{est.hyr}{effect size estimate of hybridR method}
#' \item{ci.lb.hyr}{lower bound of hybridR method's confidence interval of the 
#' effect size}
#' \item{ci.ub.hyr}{upper bound of hybridR method's confidence interval of the 
#'  effect size}
#' \item{L.0.hyr}{test statistic of hybridR method's test of null-hypothesis of
#'  no effect}
#' \item{pval.0.hyr}{p-value of hybridR method's test of null-hypothesis of no 
#' effect} 
#' \item{pval.o}{two-tailed p-value of conventional study}
#' \item{est.hy0}{effect size estimate of hybrid0 method}
#' \item{ci.lb.hy0}{lower bound of hybrid0 method's confidence interval of the 
#' effect size}
#' \item{ci.ub.hy0}{upper bound of hybrid0 method's confidence interval of the 
#' effect size}
#' \item{L.0.hy0}{test statistic of hybrid0 method's test of null-hypothesis of no
#' effect}
#' \item{pval.0.hy0}{two-tailed p-value of hybrid0 method's test of
#' null-hypothesis of no effect}
#' \item{est.fe}{effect size estimate based on fixed-effect
#' meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional
#' fixed-effect meta-analysis}
#' \item{zval.fe}{test statistic of the null-hypothesis of no effect based on
#' fixed-effect meta-analysis}
#' \item{pval.fe}{two-tailed p-value of the null-hypothesis of no effect based
#' on fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional
#' fixed-effect meta-analysis}
#' \item{ci.ub.fe}{upper bound of confidence interval based on
#' fixed-effect meta-analysis}
#' \item{est.repl}{effect size estimate of replication}
#' \item{se.repl}{standard error of replication's effect size estimate}
#' \item{ci.lb.repl}{lower bound of replication's confidence interval}
#' \item{ci.ub.repl}{upper bound of replication's confidence interval}
#' \item{stat.repl}{test statistic of replication for testing null-hypothesis of
#'  no effect}
#' \item{pval.repl}{two-tailed p-value of replication for testing
#' null-hypothesis of no effect}
#' 
#' @note The \code{control} argument in the \code{hybrid} function is an optional 
#' argument that gives the user more control over the estimation procedures. 
#' This can be especially useful if estimation of the method does not converge 
#' and NAs are returned by the function. The \code{control} argument should 
#' be specified as a list containing one or more elements. For example, 
#' \code{control = list(verbose = TRUE)} Default values are used if an element is 
#' not specified. The following elements can be specified by the user:
#' 
#' \itemize{
#' \item{\code{int:}}{ A vector of length two that indicates the lower and upper 
#' bound of the interval that is used for estimating the effect size. The effect 
#' size estimate should be included in this interval. Its default values are -10
#' for the first element and the maximum effect size of a study included in the 
#' analysis + 1 as the second element. This control argument is only applicable 
#' to the implementation of van Aert and van Assen (2018).}
#' \item{\code{est.ci:}}{ A vector of length two indicating the values that are 
#' subtracted from and added to the estimate of the effect size for computing the 
#' 95\% confidence intervals. Its default values are (50, 1). To give an example, 
#' estimates for the lower and upper bound around the effect size estimate are 
#' searched on the interval (est-50, est) and (est, est+1), respectively.}
#' \item{\code{tau2.ci}}{ A vector of length two indicating the values that are 
#' subtracted from and added to the estimate of the between-study variance for 
#' computing the 95\% confidence intervals. Its default values are (0.5, 0.5). 
#' To give an example, estimates for the lower and upper bound around the effect 
#' size estimate are searched on the interval (tau2-0.5, tau2) and (tau2, tau2+0.5), 
#' respectively.}
#' \item{\code{tol:}}{ A number indicating the desired accuracy of the estimates. 
#' Its default value is .Machine$double.eps^0.25. This control argument is only 
#' applicable to the implementation of van Aert and van Assen (2018).} 
#' \item{\code{verbose:}}{ A logical indicating whether information should be printed 
#' about the estimation procedure. Its default value is FALSE.}
#' \item{\code{par:}}{ Starting values for the optimization procedure in case of
#' the implementation of van Aert (2023). The default values are zeros.}
#' \item{\code{implementation:}}{ A character indicating whether the implementation
#' of van Aert and van Assen (2018) based on a single conventional and single replication
#' should be used (\code{implementation = "two"}) or whether the implementation 
#' of van Aert (2023) should be used that allows for including more than two studies
#' (\code{implementation = "multiple"}).}
#' \item{\code{optimizer:}}{ A character indicating the optimizer that is used if
#' the implementation of van Aert (2023) is used. The default value is "Nelder-Mead".
#' The \code{optim} function is used for optimization, so the optimization 
#' methods implemented in the \code{optim} function can be used. See the 
#' documentation of the \code{optim} function for more information.}
#' \item{\code{type:}}{ A character vector indicating whether Wald-based hypothesis
#' tests and confidence intervals are preferred (\code{type = "Wald"}) or 
#' likelihood-ratio tests and profile likelihood confidence intervals 
#' (\code{type = "profile"}). The default is "Wald/profile" which implies that 
#' Wald tests and confidence intervals are computed for the fixed effects and 
#' likelihood-ratio tests and profile likelihood confidence intervals for the 
#' between-study variance.}
#' } 
#' 
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R. C. M., & van Assen, M. A. L. M. (2018). Examining
#' reproducibility in psychology: A hybrid method for statistically combining a
#' biased original study and replication. Behavior Research Methods, 50(4): 1515-1539.
#' doi:10.3758/s13428-017-0967-6
#' 
#' @references van Aert, R.C.M. (2023). Empowering meta-analysis by taking advantage 
#' of preregistered studies: The extended hybrid meta-analysis method. Manuscript 
#' in preparation.
#'
#' @examples
#' 
#' ### Application using the implementation of van Aert (2023) with more than two studies.
#' # Note that these data come from the "metadat" R package that needs to be loaded.
#' # The data come from the meta-analyses by Lehmann et al. (2018). See the help 
#' # file of dat.lehmann2018 for more information about this meta-analysis.
#' library(metadat)
#' dat <- dat.lehmann2018
#'
#' ### Create a new object to indicate which studies were conventional studies and 
#' # which ones were preregistered
#' dat$conventional <- ifelse(dat$Preregistered == "Not Pre-Registered", 1, 0)
#'
#' ### Lehmann et al. split the analyses for males and females. We only use the data
#' # of females for this example.
#' red_romance_femalep <- dat[dat$Gender == "Females", ]
#'
#' ### Prepare data for the analysis
#' yi <- red_romance_femalep$yi
#' vi <- red_romance_femalep$vi
#' conventional <- red_romance_femalep$conventional
#'
#' ### Apply the hybrid method with Wald-type hypothesis tests and confidence intervals
#' hybrid(yi = yi, vi = vi, conventional = conventional, side = "right")
#'
#' ### Apply the hybrid method with likelihood ratio hypothesis tests and profile 
#' # likelihood confidence intervals
#' hybrid(yi = yi, vi = vi, conventional = conventional, side = "right", 
#' control = list(type = "profile"))
#'
#' ### Include Color_Match as moderator
#' Color_Match <- red_romance_femalep$Color_Match
#'
#' ### Apply the hybrid method with Wald-type hypothesis tests and confidence intervals
#' hybrid(yi = yi, vi = vi, conventional = conventional, side = "right", mods = ~ Color_Match)
#'
#' ### Apply the hybrid method with likelihood ratio hypothesis tests and profile 
#' # likelihood confidence intervals
#' hybrid(yi = yi, vi = vi, conventional = conventional, side = "right", mods = ~ Color_Match, 
#' control = list(type = "profile"))
#' 
#' ### Application using the implementation of van Aert and van Assen (2018). 
#' # The hybrid method is applied to example on page 5 of van Aert and van Assen 
#' # (2018).
#'
#' pval <- c(0.03, 0.3) # p-value conventional study and replication
#' n1i <- n2i <- c(40, 80) # Sample sizes per group 
#' tobs <- qt(pval/2, df = n1i+n2i-2, lower.tail = FALSE) # Observed t-values 
#' 
#' ### Apply hybrid method using the implementation of van Aert and van Assen (2018)
#' hybrid(tobs = tobs, n1i = n1i, n2i = n2i, side = "right", conventional = c(1,0), 
#' control = list(implementation = "two"))
#'  
#' @export

hybrid <- function(m1i, m2i, mi, ri, sd1i, sd2i, sdi, n1i, n2i, ni, tobs, yi, vi, 
                   conventional, mods = NULL, alpha = 0.05, side, control) 
{
  
  ### Necessary to pass R CMD check otherwise the following warning will be present:
  # no visible binding for global variable est.fe etc.
  est.fe <- ci.ub.fe <- ci.lb.fe <- zval.fe <- NULL 
  
  if (!missing("mi") | !missing("tobs") | !missing("m1i") | !missing("ri"))
  { # If only two studies are supplied as input
    
    if (!missing("mi") & !missing("ni") & !missing("sdi")) 
    {
      measure <- "M"
      es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$conventional <- c(1, 0) # Create variable to indicate that first study is an conventional study
      
      res.repl <- repl(es = es, mi = mi, sdi = sdi, ni = ni, measure = measure,
                       side = side)
      
    } else if (!missing("ni") & !missing("tobs")) 
    {
      measure <- "MT"
      es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$conventional <- c(1, 0) # Create variable to indicate that first study is an conventional study
      
      res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
      
    } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
               !missing("sd1i") & !missing("sd2i")) 
    {
      measure <- "MD"
      es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                      sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
      
      es$conventional <- c(1, 0) # Create variable to indicate that first study is an conventional study
      
      res.repl <- repl(es = es, m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i,
                       sd1i = sd1i, sd2i = sd2i, measure = measure, side = side)
      
    } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) 
    {
      measure <- "MDT"
      es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2,
                      side = side, measure = measure)
      
      es$conventional <- c(1, 0) # Create variable to indicate that first study is an conventional study
      
      res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
      
    } else if (!missing("ri") & !missing("ni")) 
    {
      measure <- "COR"
      es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$conventional <- c(1, 0) # Create variable to indicate that first study is an conventional study
      
      res.repl <- repl(es = es, measure = measure, side = side)
    }
    
    if (es$pval[1] > alpha/2) 
    {
      stop("Conventional study is not statistically significant")
    }
    
    ### If data of more than one conventional study and replication are provided
  } else
  {
    if (!missing("mi") & !missing("ni") & !missing("sdi")) 
    {
      measure <- "M"
      
      ### Compute effect size
      es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side, 
                      measure = measure)
      es$conventional <- conventional
      
    } else if (!missing("ni") & !missing("tobs")) 
    {
      measure <- "MT"
      
      ### Compute effect size
      es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
      es$conventional <- conventional
      
    } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & 
               !missing("n2i") & !missing("sd1i") & !missing("sd2i")) 
    {
      measure <- "MD"
      
      ### Compute effect size
      es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                      sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
      es$conventional <- conventional
      
    } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) 
    {
      measure <- "MDT"
      
      ### Compute effect size
      es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2,
                      side = side, measure = measure)
      es$conventional <- conventional
      
    } else if (!missing("ri") & !missing("ni")) 
    {
      measure <- "COR"
      
      ### Compute effect size
      es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side,
                      measure = measure)
      es$conventional <- conventional
      
    } else if (!missing("yi") & !missing("vi"))
    {
      measure <- "SPE"
      
      ### Compute effect size
      es <- escompute(yi = yi, vi = vi, alpha = alpha/2, side = side, 
                      measure = measure)
      es$conventional <- conventional
      
    }
  }
  
  k <- nrow(es) # Total number of effect sizes
  k.conventional <- sum(es$conventional == 1) # Number of conventional studies in the meta-analysis
  
  ### Critical value
  es$ycv <- es$zcv*sqrt(es$vi)
  
  ### Number of fixed effects parameters to be estimated
  n_bs <- ifelse(is.null(mods), 1, ncol(model.matrix(mods, data = es)))
  
  ##############################################################################
  
  ### Default values for estimation procedures
  
  # - int            = interval that is used for estimating ES and computing CI
  # - est.ci         = values that are subtracted from (first element) and added to 
  #                    (second element) of the estimate to determine the search interval 
  #                    for computing the CI. For example, the lower bound of the CI is
  #                    searched for on the interval c(est-50, est)
  # - tau2.ci        = values that are subtracted from (first element) and added to 
  #                    (second element) of the estimate to determine the search interval 
  #                    for computing the CI. For example, the upper bound of the CI is
  #                    searched for on the interval c(tau2, tau2+0.5)
  # - tol            = the convergence tolerance used for implementation == "two"
  # - verbose        = if verbose = TRUE output is printed about the estimation procedures
  # - par            = starting values for ML estimation
  # - implementation = "multiple" or "two" indicating whether the implementation
  #                    for two studies as outlined in van Aert and van Assen
  #                    (2018) should be used or not
  # - optimizer      = optimizer that is used for ML estimation in case of 
  #                    implementation == "multiple"
  # - type           = type of hypothesis testing procedure and procedure for 
  #                    creating confidence intervals. Options are "Wald" or "profile"
  
  con <- list(int = c(-10, max(es$yi + 1)),
              est.ci = c(50, 1),
              tau2.ci = c(0.5, 0.5),
              tol = .Machine$double.eps^0.25,
              verbose = FALSE,
              par = rep(0, n_bs+1), # Starting values for ML estimation. +1 for estimating tau^2
              implementation = "multiple",
              optimizer = "Nelder-Mead",
              type = "Wald/profile") # If Wald is used for fixed effects and profile for tau^2 
  
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
  
  ### Apply hybrid method
  res1 <- hy(es = es, measure = measure, side = side, mods = mods, n_bs = n_bs,
             con = con)
  
  # If the implementation of van Aert and van Assen (2018) is used with only two
  # studies (one conventional and one replication study), compute Hybrid^R,
  # Hybrid^0, and fixed-effect meta-analysis
  if (con$implementation == "two")
  {
    
    ### Conduct fixed-effect meta-analysis
    res.repl <- repl(es = es, measure = measure, side = side)
    
    ### Hybrid^R method
    if (res.repl$pval.o < alpha/2) 
    { # Use results of hybrid if two-tailed p-value of conventional studies < alpha/2
      res2 <- data.frame(est.hyr = res1$est, ci.lb.hyr = res1$ci.lb,
                         ci.ub.hyr = res1$ci.ub, L.0.hyr = res1$L.0, 
                         pval.0.hyr = res1$pval.0, pval.o = res.repl$pval.o)
    } else 
    { # Use results of only replications if two-tailed p-value of conventional studies 
      # > alpha/2
      res2 <- data.frame(est.hyr = res.repl$est.repl, ci.lb.hyr = res.repl$ci.lb.repl,
                         ci.ub.hyr = res.repl$ci.ub.repl, L.0.hyr = res.repl$stat.repl,
                         pval.0.hyr = res.repl$pval.repl, pval.o = res.repl$pval.o)
    }
    
    ### Hybrid^0 method
    res3 <- hy0(es = es, res1 = res1, alpha = alpha/2)
    
    ### Apply fixed-effect meta-analysis
    res4 <- fe_ma(yi = es$yi, vi = es$vi)
    
    ### Transform results of fixed-effect meta-analysis
    if (measure == "COR") 
    { # Back transform Fisher z to correlation
      res4$est.fe <- (exp(2*res4$est.fe) - 1)/(exp(2*res4$est.fe) + 1)
      res4$ci.lb.fe <- (exp(2*res4$ci.lb.fe) - 1)/(exp(2*res4$ci.lb.fe) + 1)
      res4$ci.ub.fe <- (exp(2*res4$ci.ub.fe) - 1)/(exp(2*res4$ci.ub.fe) + 1)
      res4$zval.fe <- res4$zval.fe
      
      if (side == "left") 
      { # Re-mirror estimates
        res4$est.fe <- est.fe*-1
        tmp <- ci.ub.fe
        res4$ci.ub.fe <- ci.lb.fe*-1
        res4$ci.lb.fe <- tmp*-1
        res4$zval.fe <- zval.fe*-1
      }
      
    } else if (side == "left" & measure != "COR") 
    { # Re-mirror estimates
      res4$est.fe <- res4$est.fe*-1
      tmp <- res4$ci.ub.fe
      res4$ci.ub.fe <- res4$ci.lb.fe*-1
      res4$ci.lb.fe <- tmp*-1
      res4$zval.fe <- res4$zval.fe*-1
    } else 
    {
      res4$est.fe <- res4$est.fe
      res4$ci.ub.fe <- res4$ci.ub.fe
      res4$ci.lb.fe <- res4$ci.lb.fe
      res4$zval.fe <- res4$zval.fe
    }
  } else
  { # Assign NAs to objects if implementation == "multiple"
    res2 <- data.frame(est.hyr = NA, ci.lb.hyr = NA, ci.ub.hyr = NA, L.0.hyr = NA,
                       pval.0.hyr = NA, pval.o = NA)
    
    res3 <- data.frame(est.hy0 = NA, ci.lb.hy0 = NA, ci.ub.hy0 = NA, 
                       L.0.hy0 = NA, pval.0.hy0 = NA, ave = NA)
    
    res4 <- data.frame(est.fe = NA, se.fe = NA, ci.lb.fe = NA, ci.ub.fe = NA, 
                       zval.fe = NA, pval.fe = NA, pval.fe.one = NA, 
                       Qstat = NA, Qpval = NA)
    
    res.repl <- data.frame(est.repl = NA, se.repl = NA, ci.lb.repl = NA, 
               ci.ub.repl = NA, stat.repl = NA, pval.repl = NA, pval.o = NA)
  }
  
  x <- list(con = con, var_names = var_names, est = res1$est, tau2 = res1$tau2, 
            se = res1$se, L.0 = res1$L.0, pval.0 = res1$pval.0, 
            L.het = res1$L.het, pval.het = res1$pval.het, 
            ci.lb = res1$ci.lb, ci.ub = res1$ci.ub, 
            tau2.lb = res1$tau2.lb, tau2.ub = res1$tau2.ub, 
            optim.info = res1$optim.info, est.hyr = res2$est.hyr, 
            ci.lb.hyr = res2$ci.lb.hyr, ci.ub.hyr = res2$ci.ub.hyr, 
            L.0.hyr = res2$L.0.hyr, pval.0.hyr = res2$pval.0.hyr, pval.o = res2$pval.o, 
            est.hy0 = res3$est.hy0, ci.lb.hy0 = res3$ci.lb.hy0, 
            ci.ub.hy0 = res3$ci.ub.hy0, L.0.hy0 = res3$L.0.hy0, 
            pval.0.hy0 = res3$pval.0.hy0, est.fe = res4$est.fe, 
            se.fe = res4$se.fe, ci.lb.fe = res4$ci.lb.fe, ci.ub.fe = res4$ci.ub.fe, 
            zval.fe = res4$zval.fe, pval.fe = res4$pval.fe, 
            est.repl = res.repl$est.repl, se.repl = res.repl$se.repl,
            ci.lb.repl = res.repl$ci.lb.repl, ci.ub.repl = res.repl$ci.ub.repl,
            stat.repl = res.repl$stat.repl, pval.repl = res.repl$pval.repl, 
            k = k, k.conventional = k.conventional, optim.info = res1$optim.info)
  class(x) <- "hybridoutput"
  return(x)
}