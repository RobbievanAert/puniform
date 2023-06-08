### Function for testing the null hypothesis of no between-study variance
testhetero <- function(yi, vi, est, tau.est, ycv, method, boot, con)
{
  reps <- con$reps
  
  if (method == "ML")
  {
    
    pval.boot <- NA # P-value is not bootstrapped with method == ML
    
    ### Log-likelihood when estimating both parameters
    ll <- ml_star(par = c(est, tau.est), yi = yi, vi = vi, ycv = ycv) 
    
    ### Log-likelihood when tau is constrained to zero
    ll0 <- optimize(f = ml_star_est, interval = con$int, tau = 0, yi = yi, vi = vi,
                    ycv = ycv, maximum = TRUE)$objective
    
    ### Conduct likelihood-ratio test
    L.het <- -2*(ll0-ll)
    
    ### 0.5 x chisq, because the tested null-hypothesis H0: tau2 = 0 is on the 
    # boundary of the parameter space. See Andrews (2001) and Molenberghs and 
    # Verbeke (2012)
    pval.het <- 0.5*pchisq(L.het, df = 1, lower.tail = FALSE)
    
  } else if (method == "P" | method == "LNP")
  {
    
    k <- length(yi) # Number of effect sizes
    
    ### Estimate effect size with tau=0
    est0 <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(-4, 4), tau = 0, 
                                         yi = yi, vi = vi, param = "est", ycv = ycv, 
                                         method = method, val = "es", cv_P = 0)$root, 
                                 silent = TRUE))
    
    if (inherits(est0, what = "try-error")) 
    {
      est0 <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(-10, 10), tau = 0, 
                                           yi = yi, vi = vi, param = "est", ycv = ycv, 
                                           method = method, val = "es", cv_P = 0)$root, 
                                   silent = TRUE))
    }
    
    if (inherits(est0, what = "try-error"))
    { # If effect size cannot be estimated, return NA
      L.het <- NA
      pval.het <- NA
      pval.boot <- NA
    } else 
    {
      ### Compute conditional probabilities at est0
      tr.q <- trq(est = est0, tau = 0, yi = yi, vi = vi, ycv = ycv, param = "est")
      
      het.q <- 2*abs(tr.q-0.5) # Compute heterogeneity statistic
      L.het <- sum(-log(1-het.q))
      pval.het <- pgamma(L.het, k, 1, lower.tail = FALSE)
      
      if (boot == TRUE)
      { # If boot == TRUE, bootstrapped p-value is computed
        
        ### Conduct bootstrapping
        L.het.boot <- replicate(reps, expr = boot_het(k = k, est0 = est0, vi = vi, 
                                                      ycv = ycv, method = method,
                                                      con = con))
        
        ### Compute p-value with bootstrapping
        pval.boot <- length(L.het.boot[L.het.boot > L.het & !is.na(L.het.boot)])/reps
        
      } else 
      {
        pval.boot <- NA
      }
    }
    
  }
  
  return(data.frame(L.het = L.het, pval.het = pval.het, pval.boot = pval.boot))
}