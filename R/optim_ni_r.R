### Function that is optimized to obtain required sample size in replication
optim_ni_r <- function(ni.r, perc, true.es, yi.o, vi.o, cv.o, measure, des.pprob, des.pow) {

  if (measure == "COR") { # If measure is correlation coefficient
    vi.r.0 <- vi.r.sm <- vi.r.me <- vi.r.la <- 1/(ni.r-3) # Compute sampling variance

    es.sm <- fis_trans(r=0.1)
    es.me <- fis_trans(r=0.3)
    es.la <- fis_trans(r=0.5)

  } else if (measure == "MD" | measure == "MDT") { # If measure is mean difference

    es.sm <- 0.2
    es.me <- 0.5
    es.la <- 0.8

    J <- 1 - 3/(4 * (ni.r/2 + ni.r/2 - 2) - 1)
    vi.r.0 <- 1/(ni.r/2) + 1/(ni.r/2) + (1 - (ni.r/2 + ni.r/2 - 2 - 2)/((ni.r/2 + ni.r/2 - 2) * J^2)) * 0^2
    vi.r.sm <- 1/(ni.r/2) + 1/(ni.r/2) + (1 - (ni.r/2 + ni.r/2 - 2 - 2)/((ni.r/2 + ni.r/2 - 2) * J^2)) * es.sm^2
    vi.r.me <- 1/(ni.r/2) + 1/(ni.r/2) + (1 - (ni.r/2 + ni.r/2 - 2 - 2)/((ni.r/2 + ni.r/2 - 2) * J^2)) * es.me^2
    vi.r.la <- 1/(ni.r/2) + 1/(ni.r/2) + (1 - (ni.r/2 + ni.r/2 - 2 - 2)/((ni.r/2 + ni.r/2 - 2) * J^2)) * es.la^2

  }

  ### Distribution of effect sizes with sampling variance of replication
  dist.0 <- qnorm(perc, mean=true.es, sd=sqrt(vi.r.0))
  dist.sm <- qnorm(perc, mean=true.es, sd=sqrt(vi.r.sm))
  dist.me <- qnorm(perc, mean=true.es, sd=sqrt(vi.r.me))
  dist.la <- qnorm(perc, mean=true.es, sd=sqrt(vi.r.la))

  if (missing(yi.o) == FALSE) { # Compute likelihoods with original study
    f.0 <- (dnorm(yi.o, mean=0, sd=sqrt(vi.o))/pnorm(cv.o, mean=0, sd=sqrt(vi.o),lower.tail=FALSE))*
      dnorm(dist.0, mean=0, sd=sqrt(vi.r.0))
    f.sm <- (dnorm(yi.o, mean=es.sm, sd=sqrt(vi.o))/pnorm(cv.o, mean=es.sm, sd=sqrt(vi.o),lower.tail=FALSE))*
      dnorm(dist.sm, mean=es.sm, sd=sqrt(vi.r.sm))
    f.me <- (dnorm(yi.o, mean=es.me, sd=sqrt(vi.o))/pnorm(cv.o, mean=es.me, sd=sqrt(vi.o),lower.tail=FALSE))*
      dnorm(dist.me, mean=es.me, sd=sqrt(vi.r.me))
    f.la <- (dnorm(yi.o, mean=es.la, sd=sqrt(vi.o))/pnorm(cv.o, mean=es.la, sd=sqrt(vi.o),lower.tail=FALSE))*
      dnorm(dist.la, mean=es.la, sd=sqrt(vi.r.la))

  } else if (missing(yi.o) == TRUE) { # Compute likelihoods without original study
    f.0 <- dnorm(dist.0, mean=0, sd=sqrt(vi.r.0))
    f.sm <- dnorm(dist.sm, mean=es.sm, sd=sqrt(vi.r.sm))
    f.me <- dnorm(dist.me, mean=es.me, sd=sqrt(vi.r.me))
    f.la <- dnorm(dist.la, mean=es.la, sd=sqrt(vi.r.la))
  }

  ### Posterior probabilities of hypotheses
  if (true.es == 0) { pprob <- f.0/(f.0+f.sm+f.me+f.la)
  } else if (true.es == fis_trans(r=0.1) | true.es == 0.2) { pprob <- f.sm/(f.0+f.sm+f.me+f.la)
  } else if (true.es == fis_trans(r=0.3) | true.es == 0.5) { pprob <- f.me/(f.0+f.sm+f.me+f.la)
  } else if (true.es == fis_trans(r=0.5) | true.es == 0.8) { pprob <- f.la/(f.0+f.sm+f.me+f.la) }

  mean(pprob > des.pprob)-des.pow
}
