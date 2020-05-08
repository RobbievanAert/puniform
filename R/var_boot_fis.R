#' var_boot_fis
#'
#' Function for non-parametric bootstrapping procedure to compute the variance of 
#' two overlapping Fisher-z transformed correlation coefficients.
#'
#' @param ri A vector with Pearson correlation coefficients in a primary study 
#' (see Details)
#' @param n A numerical value specifying the total sample size of a primary study
#' @param r A numerical value specifying the Pearson correlation coefficient 
#' between variables h and m (see Details)
#' @param dv An integer specifying the total number of dependent measures (default 
#' is 10, see Details) 
#' @param reps An integer specifying the number of boostrap replications (default 
#' is 1,000)
#'
#' @details In case of three variables (l, h, and m), overlapping Fisher-z 
#' transformed correlation coefficients can be computed between variables l and h 
#' and variables l and m. The function computes the variance of the two overlapping 
#' Fisher-z transformed correlation coefficients using a non-parametric bootstrap 
#' procedure. For more information see van Aert & Wicherts (2020).
#' 
#' The vector \code{ri} can contain a single Pearson correlation coefficient or 
#' multiple coefficients if information on more than one dependent measure is 
#' available. The integer \code{dv} is an optional argument to specify the expected 
#' number of dependent measures used in a primary study. This argument can be any 
#' value between 2 and infinity. Larger values yield more accurate estimates of 
#' the variance but slow down the bootstrap procedure.
#'  
#' The variance that is computed with this function can be used to correct for 
#' outcome reporting bias by including the variance as a moderator in a 
#' (multivariate) meta-analysis. Please see van Aert & Wicherts (2020) for 
#' more information.
#'
#' @return The \code{var_boot_fis} function returns a numerical value that is the 
#' variance of two overlapping Fisher-z transformed correlations.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & Wicherts, J.M. (2020). Correcting for outcome 
#' reporting bias in a meta-analysis: A meta-regression approach. In preparation.
#'
#' @examples ### Compute variance for an artificial example
#' var_boot_fis(ri = 0, sigma2 = 1, n = 100, r = 0.3)
#'
#' @export
#'
var_boot_fis <- function(ri, n, r, dv = 10, reps = 1000)
{
  
  ### Mean of the correlation coefficients (best estimate of the effect)
  r_thetai <- mean(ri)
  
  ### Create variance-covariance matrix
  ### Covariance among dependent variables equal r, because variance of scores 
  # in the population (sigma2) is assumed to be 1 which does not affect the 
  # results
  cov <- r 
  Sigma <- matrix(cov, nrow = dv+1, ncol = dv+1)
  
  ### Relationship with variable that will always be included when computing a
  # correlation
  Sigma[1, ] <- r_thetai
  Sigma[ ,1] <- r_thetai
  
  diag(Sigma) <- 1
  
  ### Bootstrapping using C++ function
  boot_var <- get_var_boot_fis(Sigma, n, dv, reps)
  
  return(boot_var)
}