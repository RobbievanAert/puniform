### Function used as input for the bisection method
pdist_hy <- function(d, yi, vi, zval, zcv, k, val, cv.P) { 
  
  ### Transform d to zd for approximation
  zd <- d/sqrt(vi)
  
  q <- numeric(k)    ### Empty object for storing transformed p-values
  
  if (zcv[1] - zd[1] <= 38) { # If no extreme probability has to be computed for the original study  
    pmarg <- exp(pnorm(zcv[1]*sqrt(vi[1]), d, sqrt(vi[1]), lower.tail = FALSE, log.p = TRUE))
    ph1 <- exp(pnorm(yi[1], d, sqrt(vi[1]), lower.tail = FALSE, log.p = TRUE))
    q[1] <- ph1/pmarg
  } else { q[1] <- approx(zd[1], zval[1], zcv[1]) } # Use approximation to compute extreme probability for original study
  
  q[2] <- exp(pnorm(yi[2], d, sqrt(vi[2]), lower.tail = FALSE, log.p = TRUE))  
  
  stat <- sum(q)
  
  if (val == "est") { out <- stat - k }
  else if (val == "ci.ub") { out <- stat - cv.P }
  else if (val == "ci.lb") { out <- stat - cv.P }      
  
  out
}