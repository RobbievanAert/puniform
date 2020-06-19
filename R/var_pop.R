#' var_pop
#'
#' Function for estimating the population variance of correlated outcomes' 
#' effect size within a primary study.
#'
#' @param v A numerical value specifying the sampling variance of the effect size
#' (see Details)
#' @param r A numerical value specifying the Pearson correlation coefficient 
#' between the outcomes (see Details)
#'
#' @details This function estimates the population variance of the effect size of 
#' correlated outcomes within a study. That is, it estimates the population 
#' variance from a single draw of a multivariate normal distribution. The function 
#' assumes equal true effect size of all outcomes, equal sampling variances of 
#' the outcomes' effect size, and equal correlation (i.e., \code{r}) among the 
#' outcomes. 
#' 
#' For a derivation of this estimator see van Aert & Wicherts (2020).
#' 
#' The variance that is computed with this function can be used to correct for 
#' outcome reporting bias by including the variance as a moderator in a 
#' (multivariate) meta-analysis. Please see van Aert & Wicherts (2020) for 
#' more information.
#'
#' @return The \code{var_pop} function returns a numerical value that is the 
#' estimate of the population variance of correlated outcomes' effect size 
#' given v and r.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & Wicherts, J.M. (2020). Correcting for outcome 
#' reporting bias in a meta-analysis: A meta-regression approach. Manuscript 
#' submitted for publication.
#'
#' @examples ### Compute variance for an artificial example
#' var_pop(v = 0.1, r = 0.3)
#'
#' @export
#'
var_pop <- function(v, r)
{
  v - r*sqrt(v*v)
}
