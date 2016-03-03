pdist.hy <- function(d, yi, vi, zval, zcv, k, val, cv.P) { 
  zd <- d/sqrt(vi)
  q <- numeric(k)    
  
  if (zcv[1] - zd[1] <= 38) { 
    pmarg <- exp(pnorm(zcv[1]*sqrt(vi[1]), d, sqrt(vi[1]), lower.tail = FALSE, log.p = TRUE))
    ph1 <- exp(pnorm(yi[1], d, sqrt(vi[1]), lower.tail = FALSE, log.p = TRUE))
    q[1] <- ph1/pmarg
  } else { q[1] <- approx(zd[1], zval[1], zcv[1]) } 
  
  q[2] <- exp(pnorm(yi[2], d, sqrt(vi[2]), lower.tail = FALSE, log.p = TRUE))  
  
  stat <- sum(q)
  
  if (val == "est") { out <- stat - k }
  else if (val == "ci.ub") { out <- stat - cv.P }
  else if (val == "ci.lb") { out <- stat - cv.P }      
  
  out
}