### Function for z-value approximation
approx_puni <- function(zd, zval, zcv, method) {
  
  a <- zval - zcv
  
  if (any(method == c("LNP", "LN1MINP", "P", "KS", "AD")))
  { # If method is not ML
    oz <- 1-1/((zcv-zd)^2+2)+1/(((zcv-zd)^2+2)*((zcv-zd)^2+4))-
      5/(((zcv-zd)^2+2)*((zcv-zd)^2+4)*((zcv-zd)^2+6))
    oz.a <- 1-1/((zval-zd)^2+2)+1/(((zval-zd)^2+2)*((zval-zd)^2+4))-
      5/(((zval-zd)^2+2)*((zval-zd)^2+4)*((zval-zd)^2+6))
    out <- exp(-(zcv-zd)*a-0.5*a^2)*(zcv-zd)/(zval-zd)*oz.a/oz
    
  } else if (method == "ML") { # If method is ML
    ozcv <- 1 - 1/((zcv-zd)^2+2) + 1/(((zcv-zd)^2+2)*((zcv-zd)^2+4)) - 5/
      (((zcv-zd)^2+2)*((zcv-zd)^2+4)*((zcv-zd)^2+6))
    out <- exp(-(zcv-zd)*a)*exp(-0.5*a^2)*((zcv-zd)/ozcv)
  }
  return(out=out)
}
