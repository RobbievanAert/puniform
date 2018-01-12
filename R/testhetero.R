### Function for testing the null hypothesis of no between-study variance
testhetero <- function(yi, vi, est, tau.est, ycv, method, boot, con)
{
  
  if (method == "ML")
  {
    
    pval.boot <- NA # P-value is not bootstrapped with method == ML
    
    ### Optimize profile likelihood function of delta
    # (suppressWarnings() in order to be able to specify wide search intervals)
    est0 <- suppressWarnings(try(optimize(ml_est, c(-1,1), 0, yi, vi, ycv, 
                                          maximum = TRUE)$maximum, silent = TRUE))
    
    if (class(est0) == "try-error" | any(est0 > 0.95 | est0 < -0.95))
    {
      est0 <- suppressWarnings(try(optimize(ml_est, c(-5,5), 0, yi, vi, ycv, 
                                            maximum = TRUE)$maximum, silent = TRUE))
    }
    
    if (class(est0) == "try-error" | any(est0 > 4.95 | est0 < -4.95))
    { # If effect size could not be estimated, return NA
      L.het <- NA
      pval.het <- NA
    } else
    { # Conduct likelihood-ratio test
      L.het <- -2*(ml_est(est0, 0, yi, vi, ycv)-ml_est(est, tau.est, yi, vi, ycv))
      pval.het <- pchisq(L.het, df = 1, lower.tail = FALSE)
    }
    
  } else if (method == "P" | method == "LNP")
  {
    
    k <- length(yi) # Number of effect sizes
    
    ### Estimate effect size with tau=0
    est0 <- try(bisect_est(pdist_nsig, lo = -4, hi = 4, tol = con$tol, tau.est = 0, 
                           yi = yi, vi = vi, param = "est", ycv = ycv, method = method, 
                           val = "es"), silent = TRUE)
    
    if (class(est0) == "try-error")
    {
      est0 <- try(bisect_est(pdist_nsig, lo = -10, hi = 10, tol = con$tol, tau.est = 0,
                             yi = yi, vi = vi, param = "est", ycv = ycv, method = method, 
                             val = "es"), silent = TRUE)
    }
    
    if (class(est0) == "try-error")
    { # If effect size cannot be estimated, return NA
      L.het <- NA
      pval.het <- NA
      pval.boot <- NA
    } else 
    {
      ### Compute conditional probabilities at est0
      tr.q <- trq(est0, 0, yi, vi, ycv)
      
      het.q <- 2*abs(tr.q-0.5) # Compute heterogeneity statistic
      L.het <- sum(-log(1-het.q))
      pval.het <- pgamma(L.het, k, 1, lower.tail = FALSE)
      
      if (boot == TRUE)
      { # If boot == TRUE, bootstrapped p-value is computed
        
        ### Conduct bootstrapping
        L.het.boot <- replicate(con$reps, expr = boot_het(k = k, est0 = est0, vi = vi, 
                                                      ycv = ycv, method = method,
                                                      con = con))
        
        ### Compute p-value with bootstrapping
        pval.boot <- length(L.het.boot[L.het.boot > L.het & !is.na(L.het.boot)])/con$reps
        
      } else 
      {
        pval.boot <- NA
      }
    }
    
  }
  
  return(data.frame(L.het = L.het, pval.het = pval.het, pval.boot = pval.boot))
  
}