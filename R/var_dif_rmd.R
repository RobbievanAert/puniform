#' var_dif_rmd
#'
#' Function for computing the variance of the difference between two raw mean 
#' differences.
#'
#' @param sd1i A vector of standard deviations of the outcomes in group 
#' 1 (see Details)
#' @param sd2i A vector of standard deviations of the outcomes in group 
#' 2 (see Details)
#' @param n1i An integer specifying the sample size of group 1
#' @param n2i An integer specifying the sample size of group 2
#' @param r A numerical value specifying the Pearson correlation coefficient 
#' between participants' scores on the different outcomes
#'
#' @details Multiple raw mean differences can be computed in case of two groups 
#' and multiple outcomes. The function computes the variance of the 
#' difference of two raw mean differences given a correlation between the outcomes. 
#' For a derivation of this variance see the supplemental materials of 
#' van Aert & Wicherts (2020).
#' 
#' The vectors \code{sd1i} and \code{sd2i} can contain a single standard deviation 
#' or multiple standard deviations if information on more than one outcome 
#' measure is available.
#' 
#' The variance that is computed with this function can be used to correct for 
#' outcome reporting bias by including the variance as a moderator in a 
#' (multivariate) meta-analysis. Please see van Aert & Wicherts (2020) for 
#' more information.
#'
#' @return The \code{var_dif_rmd} function returns a numerical value that is the 
#' variance of the difference of two raw mean differences given r.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & Wicherts, J.M. (2020). Correcting for outcome 
#' reporting bias in a meta-analysis: A meta-regression approach. In preparation.
#'
#' @examples ### Compute variance for an artificial example
#' var_dif_rmd(sd1i = c(0.8, 1.2), sd2i = c(0.85, 1.15), n1i = 100, n2i = 95, r = 0.3)
#'
#' @export
#'
var_dif_rmd <- function(sd1i, sd2i, n1i, n2i, r)
{
  ### Compute pooled variance per dependent measure
  pool_var <- (sd1i^2*(n1i-1)+sd2i^2*(n2i-1))/(n1i+n2i-2)
  
  ### Mean pooled variance
  pool <- mean(pool_var)
  
  ### Compute variance of the difference (see supplemental materials of van Aert 
  # & Wicherts)
  var_dif <- (2*pool-2*pool*r)/n1i+(2*pool-2*pool*r)/n2i
  
  return(var_dif)
}