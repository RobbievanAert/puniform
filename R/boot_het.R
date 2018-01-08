### Function to compute p-value with bootstrapping for test of homogeneity
boot_het <- function(k, est0, vi, ycv, method, con) 
{ 
  
  yi.boot <- rnorm(k, mean = est0, sd = sqrt(vi)) # Sample k effect sizes given est0 and vi
  
  ### Estimate effect size with tau=0
  est.boot <- try(bisect_est(pdist_nsig, lo = -5, hi = 5, tol = con$tol, tau.est = 0, 
                             yi = yi.boot, vi = vi, param = "est", ycv = ycv, method = method, 
                             val = "es"), silent = TRUE)
  
  if (class(est.boot) == "try-error") 
  { # If effect size cannot be estimated, return NA
    stat <- NA
  } else 
  {
    ### Compute conditional probabilities at d.boot
    tr.q <- trq(est.boot, 0, yi.boot, vi, ycv)
    
    het.q.boot <- 2*abs(tr.q-0.5) # Compute heterogeneity statistic
    stat <- sum(-log(1-het.q.boot))
  }
  return(stat = stat)
}