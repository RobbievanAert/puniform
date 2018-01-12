### Function for applying bisection method for estimating effects size with P or LNP
bisect_est <- function(func, lo, hi, tol, verbose = FALSE, tau.est, yi, vi, 
                       param, ycv, method, val) 
{ 
  
  if (method == "P")
  {
    cv_P <- get_cv_P(length(yi))
  } else if (method == "LNP")
  {
    cv_P <- 0
  }
  
  flo <- func(lo, tau.est, yi, vi, param, ycv, method, val, cv_P)
  fhi <- func(hi, tau.est, yi, vi, param, ycv, method, val, cv_P)
  
  if (flo * fhi > 0)
  {
    stop("root is not included in the specified interval") 
  }
  
  chg <- hi - lo 
  while (abs(chg) > tol) { 
    root <- (lo + hi)/2 
    
    froot <- func(root, tau.est, yi, vi, param, ycv, method, val, cv_P) 
    
    if(verbose == TRUE) { cat("root = ", root, "froot =", froot,  "chg = ", chg, fill = TRUE) }
    
    if (abs(froot) <= tol) break 
    if (flo * froot < 0) hi <- root 
    if (fhi * froot < 0) lo <- root 
    chg <- hi - lo 
  } 
  
  return(root) 
} 
