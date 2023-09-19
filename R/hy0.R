### Function for applying hybrid^0 method
hy0 <- function(es, res1, alpha) 
{
  ### Compute average conditional probability
  ave <- mean(pdist_hy(d = 0, es = es, val = "est")$q)
  
  if (ave > 0.5) 
  { # If average conditional probability is larger than 0.5
    res1$est <- 0
    res1$ci.lb <- NA
    res1$ci.ub <- NA
    res1$L.0 <- 1
    res1$pval.0 <- 0.5
  }   
  
  return(data.frame(est.hy0 = res1$est, ci.lb.hy0 = res1$ci.lb, 
                    ci.ub.hy0 = res1$ci.ub, L.0.hy0 = res1$L.0, 
                    pval.0.hy0 = res1$pval.0, ave = ave))
}