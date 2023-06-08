### Function to test null hypothesis of no effect
testeffect_nsig <- function(yi, vi, est, tau.est, ycv, method, con)
{
  
  if (method == "ML")
  { # Conduct likelihood-ratio test
    
    ### Log-likelihood when estimating both parameters
    ll <- ml_star(par = c(est, tau.est), yi = yi, vi = vi, ycv = ycv) 
    
    ### Log-likelihood when average effect size is constrained to zero
    ll0 <- optimize(f = ml_star_tau, interval = con$tau.int, d = 0, yi = yi, vi = vi,
                    ycv = ycv, maximum = TRUE)$objective
    
    ### Conduct likelihood-ratio test
    L.0 <- -2*(ll0-ll)
    pval.0 <- pchisq(L.0, df = 1, lower.tail = FALSE)
  } else if (method == "P" | method == "LNP")
  {
    L.0 <- pval.0 <- NA 
  }
  
  return(data.frame(L.0 = L.0, pval.0 = pval.0))
}