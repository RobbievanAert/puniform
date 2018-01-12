### Function for computing search interval
bounds_hy <- function(yi, vi, zval, zcv, ext) {
  
  ub <- max(yi + 1)
  
  if(ext == FALSE) {
    lb <- max(yi - 38*sqrt(vi))
  } else {
    ### If search interval for bisection function is too small extend interval 
    ### to the lowest possible lower bound. Since we do not evaluate the fraction 
    ### for the replication study the extended lower bound is manually selected
    lb <- -600 
  }
  
  return(c(lb, ub))
}