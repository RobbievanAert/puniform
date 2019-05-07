### Function used as input for the bisection method
pdist_hy <- function(d, es, val, cv.P) 
{ 
  
  ### Add objects to workspace
  yi <- es$yi
  vi <- es$vi
  zval <- es$zval
  zcv <- es$zcv
  ori <- es$ori
  k <- nrow(es)
  
  zd <- d/sqrt(vi) # Transform d to zd for approximation
  
  q <- numeric(k) # Empty object for storing transformed p-values
  
  ### Loop to compute conditional probabilities
  for (i in 1:k)
  {
    if (ori[i] == 1)
    { # In case of an original study
      
      if (zcv[i] - zd[i] <= 38) 
      { # If no extreme probability has to be computed for the original study  
        pmarg <- exp(pnorm(zcv[i]*sqrt(vi[i]), d, sqrt(vi[i]), lower.tail = FALSE, 
                           log.p = TRUE))
        ph1 <- exp(pnorm(yi[i], d, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
        q[i] <- ph1/pmarg
        
      } else 
      { # Use approximation to compute extreme probability for original study
        q[i] <- approx(zd[i], zval[i], zcv[i]) 
      } 
      
    } else if (ori[i] == 0)
    { # In case of a replication
      q[i] <- exp(pnorm(yi[i], d, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))  
    }
    
  } 

  stat <- sum(q)
  
  if (val == "est") 
  { 
    out <- stat - k/2 
  } else if (val == "ci.ub") 
  { 
    out <- stat - cv.P 
  } else if (val == "ci.lb") 
  { 
    out <- stat - cv.P 
  }      
  
  return(list(out = out, q = q))
}

### Helper function for root-finding
pdist_hy_helper <- function(d, es, val, cv.P)
{
  pdist_hy(d = d, es = es, val = val, cv.P = cv.P)$out
}