### Function to compute the bounds for estimating effect size with P and LNP 
# in puni_star in such a way that the observed effect sizes are smaller (or larger) than 
# the lower bound and upper bound for computing the conditional probabilities 
# for estimating effect size and estimating tau
bounds_nsig <- function(yi, vi, ycv, method, tau.est, bounds.int)
{
  for (m in 1:(length(bounds.int)-1))
  {
    int <- c(bounds.int[m], bounds.int[m+1])
    
    ### Check if d can be estimated
    est.l <- pdist_nsig(est = int[1], tau = tau.est, yi = yi, vi = vi, 
                        param = "est", ycv = ycv, method = method, val = "es", 
                        cv_P = 0)
    est.u <- pdist_nsig(est = int[2], tau = tau.est, yi = yi, vi = vi, 
                        param = "est", ycv = ycv, method = method, val = "es", 
                        cv_P = 0)
    
    if (method == "LNP")
    {
      if (est.l > 0 & est.u < 0)
      { # Estimate of effect size is within the bounds
        break
      } 
    } else if (method == "P")
    {
      if (est.l < 0 & est.u > 0)
      { # Estimate of effect size is within the bounds
        break
      } 
    }
    
  }
  return(int)
}