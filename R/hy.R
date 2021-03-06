### Function for applying hybrid method
hy <- function(es, measure, side) 
{
  
  ### Compute bounds for effect size estimation
  bo <- bounds_hy(es = es, ext = FALSE)
  
  ### Apply bisection method for effect size
  est.hy <- try(bisect(func = pdist_hy_helper, lo = bo[1], hi = bo[2], es = es, 
                       val = "est"), silent = TRUE)
  
  if (class(est.hy) == "try-error") 
  { # Apply bisection method with extended search interval
    bo.ext <- bounds_hy(es = es, ext = TRUE)
    est.hy <- try(bisect(func = pdist_hy_helper, lo = bo.ext[1], hi = bo.ext[2], 
                         es = es, val = "est"), silent = TRUE)
    
    if (class(est.hy) == "try-error") 
    { # If estimate cannot be computed, return NA
      est.hy <- NA
    }
  }
  
  ### Apply bisection method for lower bound
  ci.lb.hy <- try(bisect(func = pdist_hy_helper, lo = bo[1], hi = est.hy, es = es, 
                         val = "ci.lb", cv.P = get_cv_P(nrow(es))), silent = TRUE)
  
  if (class(ci.lb.hy) == "try-error") 
  { # Apply bisection method with extended search interval
    bo.ext <- bounds_hy(es = es, ext = TRUE)
    ci.lb.hy <- try(bisect(func = pdist_hy_helper, lo = bo.ext[1], hi = est.hy, 
                           es = es, val = "ci.lb", cv.P = get_cv_P(nrow(es))), 
                    silent = TRUE)
    
    if (class(ci.lb.hy) == "try-error") 
    {
      ci.lb.hy <- NA
    }
  }
  
  ### Apply bisection method for upper bound
  ci.ub.hy <- try(bisect(func = pdist_hy_helper, lo = est.hy, hi = bo[2], es = es,
                         val = "ci.ub", cv.P = nrow(es) - get_cv_P(nrow(es))), 
                  silent = TRUE)
  
  if (class(ci.ub.hy) == "try-error") 
  {
    ci.ub.hy <- NA
  }
  
  if (measure == "COR") 
  { # Back transform Fisher z to correlation
    est.hy <- (exp(2 * est.hy) - 1)/(exp(2 * est.hy) + 1)
    ci.lb.hy <- (exp(2 * ci.lb.hy) - 1)/(exp(2 * ci.lb.hy) + 1)
    ci.ub.hy <- (exp(2 * ci.ub.hy) - 1)/(exp(2 * ci.ub.hy) + 1)
  }
  
  if (side == "left") 
  { # Re-mirror estimates
    est.hy <- est.hy * -1
    tmp <- ci.ub.hy
    ci.ub.hy <- ci.lb.hy * -1
    ci.lb.hy <- tmp * -1
  }
  
  ### Test of H0: mu = 0
  q <- pdist_hy(d = 0, es = es, val = "est")$q
  
  if (length(q) == 2)
  { # If there is only one original study and one replication
    x.hy <- ifelse(sum(q) < 1, sum(q), 2-sum(q))  # Compute probability density  
    pval.hy <- ifelse(sum(q) < 1, 0.5*sum(q)^2, -0.5*sum(q)^2 + 2*sum(q)-1)
    pval.hy <- ifelse(pval.hy > 0.5, (1-pval.hy)*2, pval.hy*2) # Compute two-tailed p-value  
  } else
  { # If there are more than two studies
    x.hy <- (sum(q)-nrow(es)*0.5)/sqrt(nrow(es)/12)
    pval.hy <- pnorm(x.hy)
    pval.hy <- ifelse(pval.hy > 0.5, (1-pval.hy)*2, pval.hy*2) # Compute two-tailed p-value
  }
  
  return(data.frame(est.hy = est.hy, ci.lb.hy = ci.lb.hy, ci.ub.hy = ci.ub.hy, 
                    x.hy = x.hy, pval.hy = pval.hy))
} 