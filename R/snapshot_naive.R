#' snapshot_naive
#'
#' Function for applying Snapshot Bayesian Meta-Analysis Method (snapshot naive)
#' for two-independent means and raw correlation coefficients.
#'
#' @param m1i A vector of length two containing the means in group 1 for the original
#' study and replication for two-independent means
#' @param m2i A vector of length two containing the means in group 2 for the original
#' and replication for two-independent means
#' @param n1i A vector of length two containing the sample sizes in group 1 for
#' the original study and replication for two-independent means
#' @param n2i A vector of length two containing the sample sizes in group 2 for
#' the original study and replication for two-independent means
#' @param sd1i A vector of length two containing the standard deviations in group 1
#' for the original study and replication for two-independent means
#' @param sd2i A vector of length two containing the standard deviations in group 2
#' for the original study and replication for two-independent means
#' @param ri A vector of length two containing the raw correlation coefficients
#' of the original study and replication
#' @param ni A vector of length two containing the sample size of the original 
#' study and replication for the raw correlation coefficient
#' @param tobs A vector of length two containing the t-values of the original
#' study and replication
#'
#' @details The function computes posterior probabilities (assuming a uniform prior
#' distribution) for four true effect sizes (no, small, medium, and large) based
#' on an original study and replication. For more information see van Aert and
#' van Assen (2016).
#'
#' Two different effect size measures can be used as input for the \code{snapshot.naive}
#' function: two-independent means and raw correlation coefficients.
#' Analyzing two-independent means can be done by either providing
#' the function group means (\code{m1i} and \code{m2i}), standard deviations
#' (\code{sd1i} and \code{sd2i}), and sample sizes (\code{n1i} and \code{n2i}) or
#' t-values (\code{tobs}) and sample sizes (\code{n1i} and \code{n2i}).See the Example section for
#' an example. Raw correlation coefficients can be analyzed by supplying \code{ri}
#' and \code{ni} to the \code{snapshot.naive}.
#'
#' @return The \code{snapshot.naive} function returns a data frame with posterior
#' probabilities for no (\code{p.0}), small (\code{p.sm}), medium (\code{p.me}),
#' and large (\code{p.la}) true effect size.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & van Assen, M.A.L.M. (2017). Bayesian evaluation
#' of effect size after replicating an original study. PLoS ONE, 12(4), e0175302. 
#' doi:10.1371/journal.pone.0175302
#'
#' @examples ### Example as presented on page 491 in Maxwell, Lau, and Howard (2015)
#' snapshot_naive(ri = c(0.243, 0.114), ni = c(80, 172))
#'
#' @export


snapshot_naive <- function(ri, ni, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs)
{

  alpha <- 0.05/2 # Compute alpha in two-tailed tests with results reported in predicted direction

  ### Compute standardized effect sizes
  if (!missing(ri) & !missing(ni)) { # Correlation
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha, side = "right", measure = measure)
    true.es <- c(0, fis_trans(r=0.1), fis_trans(r=0.3), fis_trans(r=0.5)) # True effect sizes
  } else if (!missing(m1i) & !missing(m2i) & !missing(n1i) & !missing(n2i) &
             !missing(sd1i) & !missing(sd2i)) { # Mean difference unknown sigma
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                    sd2i = sd2i, alpha = alpha, side = "right", measure = measure)
    true.es <- c(0, 0.2, 0.5, 0.8) # True effect sizes
  } else if (!missing(n1i) & !missing(n2i) & !missing(tobs)) { # Mean difference unknown sigma with observed t-value
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha, side = "right",
                    measure = measure)
    true.es <- c(0, 0.2, 0.5, 0.8) # True effect sizes
  }

  if (es$yi[1] < 0) { # Mirror effect size and z-value if effect size in original study is smaller than zero
    es$yi <- es$yi*-1
  }

  ### Compute probabilities
  f.0 <- dnorm(es$yi[1],true.es[1],sqrt(es$vi[1]))*dnorm(es$yi[2],true.es[1],sqrt(es$vi[2]))
  f.sm <- dnorm(es$yi[1],true.es[2],sqrt(es$vi[1]))*dnorm(es$yi[2],true.es[2],sqrt(es$vi[2]))
  f.me <- dnorm(es$yi[1],true.es[3],sqrt(es$vi[1]))*dnorm(es$yi[2],true.es[3],sqrt(es$vi[2]))
  f.la <- dnorm(es$yi[1],true.es[4],sqrt(es$vi[1]))*dnorm(es$yi[2],true.es[4],sqrt(es$vi[2]))

  ## Posterior probabilities of hypotheses
  p.0 <- f.0/(f.0+f.sm+f.me+f.la)
  p.sm <- f.sm/(f.0+f.sm+f.me+f.la)
  p.me <- f.me/(f.0+f.sm+f.me+f.la)
  p.la <- f.la/(f.0+f.sm+f.me+f.la)

  return(data.frame(p.0 = p.0, p.sm = p.sm, p.me = p.me, p.la = p.la))

}
