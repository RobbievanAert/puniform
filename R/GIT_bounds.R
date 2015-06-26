bounds <- function(yi, vi, zval, zcv, ext) {

  ub <- max(yi + 1)

  if(ext == FALSE) {
    lb <- max(yi - 38*sqrt(vi))
  } else {
    lb <- max((-(700-.5*(zval-zcv)^2)/(zval-zcv) + zcv)*sqrt(vi))
  }

  return(c(lb, ub))
}
