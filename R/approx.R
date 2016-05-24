### Function for z-value approximation
approx <- function(zd, zval, zcv) {
  a <- zval - zcv
  oz <- 1-1/((zcv-zd)^2+2)+1/(((zcv-zd)^2+2)*((zcv-zd)^2+4))-
    5/(((zcv-zd)^2+2)*((zcv-zd)^2+4)*((zcv-zd)^2+6))
  oz.a <- 1-1/((zval-zd)^2+2)+1/(((zval-zd)^2+2)*((zval-zd)^2+4))-
    5/(((zval-zd)^2+2)*((zval-zd)^2+4)*((zval-zd)^2+6))
  exp(-(zcv-zd)*a-0.5*a^2)*(zcv-zd)/(zval-zd)*oz.a/oz
}