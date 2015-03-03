#' escompute
#'
#' Function that computes Hedges' g for one-sample mean and two-sample mean and computes p-value as in primary studies.
#'
#' @param mi A vector of group means for one-sample mean
#' @param ni A vector of sample sizes for one-sample mean
#' @param sdi A vector of standard deviations for one-sample mean
#' @param m1i A vector of means in group 1 for two-sample mean
#' @param m2i A vector of means in group 2 for two-sample mean
#' @param n1i A vector of sample sizes in group 1 for two-sample mean
#' @param n2i A vector of sample sizes in group 2 for two-sample mean
#' @param sd1i A vector of standard deviations in group 1 for two-sample mean
#' @param sd2i A vector of standard deviations in group 2 for two-sample mean
#' @param side A character indicating the direction of the tested hypothesis in the primary studies (either "right" or "left")
#' @param measure A character indicating whether a Hedges' g should be computed based on one-sample mean or two-sample mean ("M" or "MD")
#'
#' @return Function returns a data frame with Hedges' g effect sizes (yi), Hedges' g variances (vi), z-values (zval), and p-values as computed in primary studies (pval).
#'
#' @export

escompute <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, side, measure) {

  if(measure == "M") {
    di <- mi/sdi
    J <- 1 - 3/(4*(ni-1)-1)
    yi <- J * di
    vi <- J^2 * (1/ni+yi^2/(2*ni*(ni-1)))
    zval <- yi/sqrt(vi)
    tval <- mi/(sdi/sqrt(ni))
    if(side == "right") { pval <- pt(tval, df = ni - 1, lower.tail = FALSE) }
    if(side == "left") { pval <- pt(tval, df = ni - 1) }
  }

  if(measure == "MD") {

    s.pool <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i+n2i-2))
    di <- (m1i-m2i)/s.pool
    J <- 1 - 3/(4*(n1i+n2i-2)-1)
    yi <- J * di
    vi <- 1/n1i+1/n2i+(1-(n1i+n2i-2)/(n1i+n2i*J^2))*yi^2
    zval <- yi/sqrt(vi)
    tval <- (m1i - m2i)/sqrt(s.pool*(1/n1i+1/n2i))
    if(side == "right") { pval <- pt(tval, df = n1i+n2i-2, lower.tail = FALSE) }
    if(side == "left") { pval <- pt(tval, df = n1i+n2i-2) }
  }

  if(measure == "COR") {

    yi <- .5*log((1 + ri) / (1 - ri))
    vi <- 1/(ni-3)
    zval <- yi/sqrt(vi)
    if(side == "right") { pval <- pnorm(zval, lower.tail = FALSE) }
    if(side == "left") { pval <- pnorm(zval) }
  }

  if(side == "left") {
    yi <- yi * -1
    zval <- zval * -1
  }

  return(data.frame(yi, vi, zval, pval))
}
