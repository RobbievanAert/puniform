#' req_ni_r
#'
#' Function for computing the required sample size for a replication based on the
#' Snapshot Bayesian Hybrid Meta-Analysis Method for two-independent
#' means and raw correlation coefficients.
#'
#' @param m1i.o An integer containing the mean in group 1 of the original
#' study for two-independent means
#' @param m2i.o An integer containing the mean in group 2 of the original
#' study for two-independent means
#' @param n1i.o An integer containing the sample size in group 1 of
#' the original study for two-independent means
#' @param n2i.o An integer containing the sample size in group 2 of
#' the original study for two-independent means
#' @param sd1i.o An integer containing the standard deviation in group 1
#' of the original study for two-independent means
#' @param sd2i.o An integer containing the standard deviation in group 2
#' of the original study for two-independent means
#' @param ri.o An integer containing the raw correlation coefficient
#' of the orignal study
#' @param ni.o An integer containing the sample size for the raw correlation 
#' coefficient
#' @param tobs.o An integer containing the t-value of the original
#' study
#' @param alpha An integer specifying the alpha level as used in the original study
#' @param des.pprob An integer specifying the posterior probablity that an user
#' desires to obtain for one of the four true effect sizes
#' @param des.pow An integer specifying the probability of observing a posterior
#' probability larger than des.pprob that an user desires to obtain for one of the
#' four true effect sizes
#' @param lo An integer specifiying the lower bound of the search interval that is
#' used for the optimization procedure (default is 4)
#' @param hi An integer specifiying the upper bound of the search interval that is
#' used for the optimization procedure (default is 100,000)
#'
#' @details The function computes the required sample size for the replication based
#' on the Snapshot Bayesian Hybrid Meta-Analysis Method for four true effect sizes
#' (no, small, medium, and large). The required sample size is computed by optimizing
#' \eqn{P(\pi_x \ge a)=b} with \eqn{\pi_x} being the posterior probability with x
#' referring to no (0), small (S), medium (M), and large (L) true effect size and \eqn{a}
#' the desired posterior probability, and \eqn{b} the desired probability of observing
#' a posterior probability larger than \eqn{a}. The required sample size for the
#' replication is computed with and without including information of the originial
#' study. Computing the required sample size with the Snapshot Bayesian Hybrid
#' Meta-Analysis Method is akin to computing the required sample size with a power
#' analysis in null hypothesis significance testing. For more information see van
#' Aert and van Assen (2016).
#'
#' The \code{req.ni.r} function assumes that a two-tailed hypothesis test was
#' conducted in the original study. In case one-tailed hypothesis tests was
#' conducted in the original study, the alpha level has to be multiplied by two.
#' For example, if a one-tailed hypothesis test was conducted with an alpha level
#' of .05, an alpha of 0.1 has to be submitted to \code{req.ni.r}.
#'
#' @return The \code{req.ni.r} function returns a 4x2 matrix with in the first
#' column the required total sample size of the replication when information of the
#' original study is taken into account and in the second column the required sample
#' size if information of the original study is ignored.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & van Assen, M.A.L.M. (2016). Bayesian evaluation
#' of effect size after replicating an original study. Manuscript submitted for
#' publication.
#'
#' @examples ### Example as presented on page 491 in Maxwell, Lau, and Howard (2015)
#' req_ni_r(ri.o = 0.243, ni.o = 80, alpha = .05, des.pprob = 0.75, des.pow = 0.8)
#'
#' @export

req_ni_r <- function(ri.o, ni.o, m1i.o, m2i.o, n1i.o, n2i.o, sd1i.o, sd2i.o, tobs.o,
                     alpha, des.pprob, des.pow, lo=4, hi=100000) {

  alpha <- alpha/2 # Compute alpha in two-tailed tests with results reported in predicted direction

  ### Compute standardized effect sizes
  if (!missing(ri.o) & !missing(ni.o)) { # Correlation
    measure <- "COR"
    es <- escompute(ri = ri.o, ni = ni.o, alpha = alpha, side = "right", measure = measure)
    true.es <- c(0, fis_trans(r=0.1), fis_trans(r=0.3), fis_trans(r=0.5)) # True effect sizes
  } else if (!missing(m1i.o) & !missing(m2i.o) & !missing(n1i.o) & !missing(n2i.o) &
             !missing(sd1i.o) & !missing(sd2i.o)) { # Mean difference unknown sigma
    measure <- "MD"
    es <- escompute(m1i = m1i.o, m2i = m2i.o, n1i = n1i.o, n2i = n2i.o, sd1i = sd1i.o,
                    sd2i = sd2i.o, alpha = alpha, side = "right", measure = measure)
    true.es <- c(0, 0.2, 0.5, 0.8) # True effect sizes
  } else if (!missing(n1i.o) & !missing(n2i.o) & !missing(tobs.o)) { # Mean difference unknown sigma with observed t-value
    measure <- "MDT"
    es <- escompute(n1i = n1i.o, n2i = n2i.o, tobs = tobs.o, alpha = alpha, side = "right",
                    measure = measure)
    true.es <- c(0, 0.2, 0.5, 0.8) # True effect sizes
  }

  if (es$yi < 0) { # Mirror effect size and z-value if effect size in original study is smaller than zero
    es$yi <- es$yi*-1
    es$zval <- es$zval*-1
  }

  points <- 1000 # Number of points that are used to get distribution of Fisher-z
  perc <- 1:points/(points+1) # Percentiles for getting distribution of Fisher-z

  vi.o <- es$vi # Store sampling variance
  yi.o <- es$yi # Store standardized effect size
  cv.o <- qnorm(1-alpha, mean=0, sd=sqrt(vi.o)) # Critical value in original study

  ni.r.o <- ni.r <- numeric(length(true.es)) # Empty objects for storing results

  for (i in 1:4) { # Compute required sample size for no, small, medium, and large effect

    if (es$zval < es$zcv) { # If original study is not significant return NA
      ni.r.o[i] <- NA
    } else {
      if(optim_ni_r(ni.r=4, perc=perc, true.es=true.es[i], yi.o=yi.o, vi.o=vi.o,
                    cv.o=cv.o, measure=measure, des.pprob=des.pprob, des.pow=des.pow) > 0) {
        ### If required sample size is smaller than 4 return "< 4" (sampling variance cannot be computed for correlation)
        ni.r.o[i] <- "< 4"
      } else { # Apply bisection to search for required sample size with original study
        tmp <- try(ceiling(bisect(optim_ni_r, lo=lo, hi=hi, perc=perc, true.es=true.es[i],
                                  yi.o=yi.o, vi.o=vi.o, cv.o=cv.o, measure=measure,
                                  des.pprob=des.pprob, des.pow=des.pow)))

        ### If required sample size could not be computed, the sample size is larger than 100,000
        if (class(tmp) == "try-error") { ni.r.o[i] <- "> 100,000"
        } else { ni.r.o[i] <- tmp }
      }
    }

    ### Apply bisection to search for required sample size without original study
    tmp <- try(ceiling(bisect(optim_ni_r, lo=lo, hi=hi, perc=perc, true.es=true.es[i],
                              measure=measure, des.pprob=des.pprob,
                              des.pow=des.pow)))

    ### If required sample size could not be computed, the sample size is larger than 100,000
    if (class(tmp) == "try-error") { ni.r[i] <- "> 100,000"
    } else { ni.r[i] <- tmp }
  }

  ### If original study is not significant return notification
  if (is.na(ni.r.o[1])) { message("Original study is not statistically significant")}

  return(noquote(matrix(c(ni.r.o, ni.r), nrow=4, ncol=2,
                        dimnames=list(c("0", "S", "M", "L"), c("ni.r.o", "ni.r")))))

}
