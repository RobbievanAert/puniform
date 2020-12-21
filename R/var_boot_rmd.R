#' var_boot_rmd
#'
#' Function for parametric bootstrapping procedure to estimate the variability in 
#' outcomes' effect size in case of raw mean difference as effect size measure.
#'
#' @param sd1i A vector of standard deviations of the outcomes in group 
#' 1 (see Details)
#' @param sd2i A vector of standard deviations of the outcomes in group 
#' 2 (see Details)
#' @param n1i An integer specifying the sample size of group 1
#' @param n2i An integer specifying the sample size of group 2
#' @param r A numerical value specifying the Pearson correlation coefficient 
#' between participants' scores on the different outcomes
#' @param dv An integer specifying the total number of outcomes (default 
#' is 10, see Details) 
#' @param reps An integer specifying the number of bootstrap replications (default 
#' is 1,000)
#'
#' @details Multiple raw mean differences can be computed in case of two groups 
#' and multiple outcomes. The function estimates the variance of raw 
#' mean differences given a correlation among the outcomes using a 
#' parametric bootstrap procedure. For more information see van Aert & 
#' Wicherts (2020).
#' 
#' The vectors \code{sd1i} and \code{sd2i} can contain a single standard deviation 
#' or multiple standard deviations if information on more than one outcome
#' is available. The integer \code{dv} is an optional argument to specify 
#' the expected number of outcomes used in a primary study. This argument 
#' can be any value between 2 and infinity. Larger values yield more accurate 
#' estimates of the variance but slow down the bootstrap procedure.
#'  
#' The variance that is estimated with this function can be used to correct for 
#' outcome reporting bias by including the variance as a moderator in a 
#' (multivariate) meta-analysis. Please see van Aert & Wicherts (2020) for 
#' more information.
#'
#' @return The \code{var_boot_rmd} function returns a numerical value that is an 
#' estimate of the variance of multiple correlated raw mean differences.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & Wicherts, J.M. (2021). Correcting for outcome 
#' reporting bias in a meta-analysis: A meta-regression approach. Manuscript 
#' submitted for publication.
#'
#' @examples ### Compute variance for an artificial example
#' var_boot_rmd(sd1i = c(0.8, 1.2), sd2i = c(0.85, 1.15), n1i = 100, n2i = 95, r = 0.3)
#'
#' @export
#'
var_boot_rmd <- function(sd1i, sd2i, n1i, n2i, r, dv = 10, reps = 1000)
{
  
  pool <- (sd1i^2*(n1i-1) + sd2i^2*(n2i-1))/(n1i+n2i-2) # Pooled variances
  vi <- pool*(1/n1i+1/n2i) # Variances of raw mean differences
  
  ### Create variance-covariance matrix. Use the mean of the variances, because 
  # this is the best estimate of the common underlying variance
  Sigma <-  matrix(r*sqrt(mean(vi)*mean(vi)), nrow = dv, ncol = dv)
  diag(Sigma) <- mean(vi)
  
  ### Bootstrapping using C++ function
  boot_var <- get_var_boot_rmd(Sigma, dv, reps)
  
  return(boot_var)
}