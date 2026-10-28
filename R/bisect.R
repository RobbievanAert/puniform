### Function for applying the bisection method
bisect <- function(func, lo, hi, tol = 0.0001, verbose = FALSE, ...) 
{ 
  
  flo <- func(lo, ...) 
  fhi <- func(hi, ...) 
  
  if (flo * fhi > 0) 
    stop("root is not included in the specified interval") 
  
  chg <- hi - lo 
  while (abs(chg) > tol) 
  { 
    root <- (lo + hi)/2 
    froot <- func(root, ...) 
    
    if(verbose == TRUE) 
    { 
      cat("root = ", root, "froot =", froot,  "chg = ", chg, fill = TRUE) 
    }
    
    if (abs(froot) <= tol) break 
    if (flo * froot < 0) hi <- root 
    if (fhi * froot < 0) lo <- root 
    chg <- hi - lo 
  } 
  
  return(root) 
} 