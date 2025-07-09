#' snapshot
#'
#' Function for applying Snapshot Bayesian Hybrid Meta-Analysis Method for two-independent
#' means and raw correlation coefficients.
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
#' @param alpha A numerical value specifying the alpha level as used in the original
#' study (default is 0.05, see Details)
#'
#' @details The function computes posterior probabilities (assuming a uniform prior
#' distribution) for four true effect sizes (no, small, medium, and large) based
#' on an original study and replication while taking into account statistical
#' significance in the original study. For more information see van Aert and van
#' Assen (2016).
#'
#' Two different effect size measures can be used as input for the \code{snapshot}
#' function: two-independent means and raw correlation coefficients.
#' Analyzing two-independent means can be done by either providing
#' the function group means (\code{m1i} and \code{m2i}), standard deviations
#' (\code{sd1i} and \code{sd2i}), and sample sizes (\code{n1i} and \code{n2i}) or
#' t-values (\code{tobs}) and sample sizes (\code{n1i} and \code{n2i}). Both options
#' should be accompanied with input for the argument \code{alpha}. See the Example section for
#' an example. Raw correlation coefficients can be analyzed by supplying \code{ri}
#' and \code{ni} to the \code{snapshot} function next to input for the argument
#' \code{alpha}.
#'
#' The \code{snapshot} function assumes that a two-tailed hypothesis test was
#' conducted in the original study. In case a one-tailed hypothesis test was
#' conducted in the original study, the alpha level has to be multiplied by two.
#' For example, if a one-tailed hypothesis test was conducted with an alpha level
#' of .05, an alpha of 0.1 has to be submitted to \code{snapshot}.
#'
#' @return The \code{snapshot} function returns a data frame with posterior probabilities
#' for no (\code{p.0}), small (\code{p.sm}), medium (\code{p.me}), and large (\code{p.la})
#' true effect size.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & van Assen, M.A.L.M. (2017). Bayesian evaluation
#' of effect size after replicating an original study. PLoS ONE, 12(4), e0175302. 
#' doi:10.1371/journal.pone.0175302
#'
#' @examples ### Example as presented on page 491 in Maxwell, Lau, and Howard (2015)
#' snapshot(ri = c(0.243, 0.114), ni = c(80, 172))
#'
#' @export

snapshot <- function(ri, ni, m1i, m2i, n1i, n2i, sd1i, sd2i, tobs, alpha = 0.05)
{

  alpha <- alpha/2 # Compute alpha in two-tailed tests with results reported in predicted direction

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

  ycvi <- qnorm(alpha, lower.tail = FALSE)*sqrt(es$vi) # Critical value

  ### Compute log likelihood given a particular effect size (i.e., snapshot) 
  f.0 <- (dnorm(es$yi[1],true.es[1],sqrt(es$vi[1]), log = TRUE) - 
                pnorm(ycvi[1],true.es[1],sqrt(es$vi[1]), lower.tail = FALSE, log.p = TRUE)) + 
    dnorm(es$yi[2],true.es[1],sqrt(es$vi[2]), log = TRUE)
  
  f.sm <- (dnorm(es$yi[1],true.es[2],sqrt(es$vi[1]), log = TRUE) - 
                 pnorm(ycvi[1],true.es[2],sqrt(es$vi[1]), lower.tail = FALSE, log.p = TRUE)) + 
    dnorm(es$yi[2],true.es[2],sqrt(es$vi[2]), log = TRUE)
  
  f.me <- (dnorm(es$yi[1],true.es[3],sqrt(es$vi[1]), log = TRUE) - 
                 pnorm(ycvi[1],true.es[3],sqrt(es$vi[1]), lower.tail = FALSE, log.p = TRUE)) + 
    dnorm(es$yi[2],true.es[3],sqrt(es$vi[2]), log = TRUE)
  
  f.la <- (dnorm(es$yi[1],true.es[4],sqrt(es$vi[1]), log = TRUE) - 
                 pnorm(ycvi[1],true.es[4],sqrt(es$vi[1]), lower.tail = FALSE, log.p = TRUE)) + 
    dnorm(es$yi[2],true.es[4],sqrt(es$vi[2]), log = TRUE)
  
  ### Compute posterior probabilities of hypotheses using the log sum exponent
  # trick to avoid numerical underflow
  fs <- c(f.0, f.sm, f.me, f.la)
  constant <- max(fs)
  ps <- exp(fs - (constant + log(sum(exp(fs - constant)))))
  
  ### Store posterior probabilities
  p.0 <- ps[1]
  p.sm <- ps[2]
  p.me <- ps[3]
  p.la <- ps[4]

  return(data.frame(p.0 = p.0, p.sm = p.sm, p.me = p.me, p.la = p.la))
}