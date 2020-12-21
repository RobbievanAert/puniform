#' var_dif_fis
#'
#' Function for computing the variance of the difference between two overlapping 
#' Fisher-z transformed correlation coefficients.
#'
#' @param n A numerical value specifying the total sample size of a primary study
#' @param r A numerical value specifying the Pearson correlation coefficient 
#' between variables h and m (see Details)
#' @param rho A numerical value specifying the Pearson correlation coefficient 
#' between variables l and h and variables h and m (see Details)
#'
#' @details In case of three variables (l, h, and m), overlapping Fisher-z 
#' transformed correlation coefficients can be computed between variables l and h 
#' and variables l and m. The function computes the variance of the difference 
#' between these two overlapping Fisher-z transformed correlations. For a derivation 
#' of this variance see van Aert & Wicherts (2020).
#' 
#' The variance that is computed with this function can be used to correct for 
#' outcome reporting bias by including the variance as a moderator in a 
#' (multivariate) meta-analysis. Please see van Aert & Wicherts (2020) for 
#' more information.
#'
#' @return The \code{var_dif_fis} function returns a numerical value that is the 
#' variance of the difference of two overlapping Fisher-z transformed correlations 
#' given n, r, and rho.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & Wicherts, J.M. (2021). Correcting for outcome 
#' reporting bias in a meta-analysis: A meta-regression approach. Manuscript 
#' submitted for publication.
#'
#' @examples ### Compute variance for an artificial example
#' var_dif_fis(n = 100, r = 0.3, rho = 0.5)
#'
#' @export
#'
var_dif_fis <- function(n, r, rho)
{
  2/(n-3) - 2*((r*(-2*rho^2+1)-0.5*rho^2*(-r^2+1-2*rho^2))/(1-rho^2)^2/(n-3))
}
