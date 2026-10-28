### Function for obtaining search interval for bisection function
bounds <- function(yi, vi, zval, zcv, ext) {

    ub <- max(yi + 3)

    if(ext == FALSE) {
      lb <- max(yi - 38*sqrt(vi))
    } else {
      ### If search interval for bisection function is too small extend interval
      # to the lowest possible lower bound
      lb <- max((-(700-0.5*(zval-zcv)^2)/(zval-zcv)+zcv)*sqrt(vi))
    }

  return(c(lb, ub))
}
