### Function for publication bias test p-uniform*
pubbias_nsig <- function(yi, vi, ycv, est, tau.est, method) 
{
  
  ### Random-effects meta-analysis
  res.re <- metafor::rma(yi = yi, vi = vi, method = "PM")
  
  if (method == "ML")
  {
    
    L.pb <- (2*(ml_est(est, tau.est, yi, vi, ycv) - 
                  ml_est(res.re$b[1], sqrt(res.re$tau2), yi, vi, ycv)))
    
    pval.pb <- pchisq(L.pb, df = 2, lower.tail = FALSE)
    
  } else if (method == "P" | method == "LNP")
  {
    L.pb <- pval.pb <- NA
  }
  
  return(list(L.pb = L.pb, pval.pb = pval.pb))
}