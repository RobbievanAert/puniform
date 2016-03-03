bounds.hy <- function(yi, vi, zval, zcv, ext) {
  
  ub <- max(yi + 1)
  
  if(ext == FALSE) {
    lb <- max(yi - 38*sqrt(vi))
  } else {
    lb <- -600
  }
  
  return(c(lb, ub))
}