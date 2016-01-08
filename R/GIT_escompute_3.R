#' escompute
#'
#' Function that computes Hedges' g for an one-sample mean and a two-independent means, Fisher's r-to-z transformed correlation coefficient for a raw correlation coefficient and computes a p-value as in primary studies.
#'
#' @param mi A vector of group means for one-sample mean
#' @param ni A vector of sample sizes for one-sample mean
#' @param sdi A vector of standard deviations for one-sample mean
#' @param m1i A vector of means in group 1 for two-independent means
#' @param m2i A vector of means in group 2 for two-independent means
#' @param n1i A vector of sample sizes in group 1 for two-independent means
#' @param n2i A vector of sample sizes in group 2 for two-independent means
#' @param sd1i A vector of standard deviations in group 1 for two-independent means
#' @param sd2i A vector of standard deviations in group 2 for two-independent means
#' @param ri A vector of raw correlation coefficients
#' @param tobs A vector of t-values
#' @param alpha A integer specifying the alpha level as used in primary studies
#' @param side A character indicating the direction of the tested hypothesis in the primary studies (either "\code{right}" or "\code{left}")
#' @param measure A character indicating what kind of effect size should be computed (Hedges' g or Fisher's r-to-z transformed correlation coefficients) and which arguments are used as input ("\code{M}", "\code{MT}", "\code{MD}", "\code{MDT}", or "\code{COR}"). See Details below.
#'
#' @details The \code{measure} argument has to be used to specify the desired effect size and what input parameters are used. There are five options:
#' \itemize{
#' \item{\code{"M"}}{ for one-sample mean with \code{mi}, \code{ni}, \code{sdi}, \code{alpha}, and \code{side} as input parameters}
#' \item{\code{"MT"}}{ for one-sample mean with \code{tobs}, \code{ni}, \code{alpha}, and \code{side} as input parameters}
#' \item{\code{"MD"}}{ for two-sample mean with \code{m1i}, \code{m2i}, \code{n1i}, \code{n2i}, \code{sd1i}, \code{sd2i}, \code{alpha}, and \code{side} as input parameters}
#' \item{\code{"MDT"}}{ for two-sample mean with \code{tobs}, \code{n1i}, \code{n2i}, \code{alpha}, and \code{side} as input parameters}
#' \item{\code{"COR"}}{ for raw correlation coefficients with \code{ri}, \code{ni}, \code{alpha}, and \code{side} as input parameters}
#' }
#'
#' @return Function returns a data frame with standardized effect sizes (yi), variances of these standardized effect sizes (vi), z-values (zval), p-values as computed in primary studies (pval), and critical z-values (zcv).
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @export

escompute <- function(mi, ri, ni, sdi, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, alpha, side, measure) {

  if (measure == "M" | measure == "MT") {
    if (measure == "M") {
      di <- mi/sdi
      tval <- mi/(sdi/sqrt(ni))
    } else if (measure == "MT") {
      di <- tobs*(1/sqrt(ni))
      tval <- tobs
    }
    J <- 1 - 3/(4*(ni-1)-1)
    yi <- J * di
    vi <- J^2 * (1/ni+di^2/(2*ni))
    zval <- yi/sqrt(vi)
    if (side == "right") { pval <- pt(tval, df = ni - 1, lower.tail = FALSE) }
    if (side == "left") { pval <- pt(tval, df = ni - 1) }
    dcvi <- qt(alpha, df = ni-1, lower.tail = FALSE)*1/sqrt(ni)
    ycvi <- J * dcvi
    zcv <- ycvi/sqrt(vi)

  }

  if (measure == "MD" | measure == "MDT") {
    if (measure == "MD") {
      s.pool <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i+n2i-2))
      di <- (m1i-m2i)/s.pool
      tval <- (m1i - m2i)/sqrt(s.pool^2*(1/n1i+1/n2i))
    } else if (measure == "MDT") {
      di <- tobs*sqrt((n1i+n2i)/(n1i*n2i))
      tval <- tobs
    }
    J <- 1 - 3/(4*(n1i+n2i-2)-1)
    yi <- J * di
    vi <- 1/n1i+1/n2i+(1-(n1i+n2i-2-2)/((n1i+n2i-2)*J^2))*yi^2
    zval <- yi/sqrt(vi)
    if (side == "right") { pval <- pt(tval, df = n1i+n2i-2, lower.tail = FALSE) }
    if (side == "left") { pval <- pt(tval, df = n1i+n2i-2) }
    dcvi <- qt(alpha, df = n1i+n2i-2, lower.tail = FALSE)*sqrt((n1i+n2i)/(n1i*n2i))
    ycvi <- J * dcvi
    zcv <- ycvi/sqrt(vi)
  }

  if (measure == "COR") {
    yi <- .5*log((1 + ri) / (1 - ri))
    vi <- 1/(ni-3)
    zval <- yi/sqrt(vi)
    zcv <- qnorm(alpha, lower.tail = FALSE)
    if (side == "right") { pval <- pnorm(zval, lower.tail = FALSE) }
    if (side == "left") { pval <- pnorm(zval) }
  }

  if (side == "left") {
    yi <- yi * -1
    zval <- zval * -1
  }

  return(data.frame(yi, vi, zval, pval, zcv))
}
