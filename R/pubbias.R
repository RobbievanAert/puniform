### Function for publication bias test p-uniform
pubbias <- function(yi, vi, zval, zcv, ksig, alpha, method, est.fe, est) {
  
  ### Check if there are significant studies
  if (ksig == 0)
  {
    L.pb <- NA
    pval.pb <- NA
    approx.pb <- NA
  } else {
    if (method == "KS" | method == "AD")
    {
      ### If method is KS or AD return NAs for publication bias test
      L.pb <- NA
      pval.pb <- NA
      approx.pb <- 0 # Create object for notification in output
    } else if (method == "ML")
    { # Likelihood-ratio test
      L.pb <- LR_test(d.alt = est, d.null = est.fe, yi = yi, vi = vi, zcv = zcv)
      pval.pb <- pchisq(L.pb, df = 1, lower.tail = FALSE)
      approx.pb <- 0 # Create object for notification in output
    } else 
    { 
      zd <- est.fe/sqrt(vi) # Transform est.fe to zd for approximation
      q <- numeric(ksig)  # Empty object for storing transformed p-values
      ### Loop for computing transformed p-values
      for (i in 1:length(yi))
      {
        if (zcv[i] - zd[i] <= 38)
        { # Do not use approximation
          approx.pb <- 0 # Create object for notification in output
          pmarg <- exp(pnorm(zcv[i] * sqrt(vi[i]), est.fe, sqrt(vi[i]),
                             lower.tail = FALSE, log.p = TRUE))
          ph1 <- exp(pnorm(yi[i], est.fe, sqrt(vi[i]), lower.tail = FALSE,
                           log.p = TRUE))
          q[i] <- ph1/pmarg
        } else if (zd[i] > -(700-0.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i])+zcv[i])
        {
          ### If statement above is TRUE: approximate P(Z>=zval) and P(Z>=zcv)
          approx.pb <- 1 # Create object for notification in output
          q[i] <- approx_puni(zd[i], zval[i], zcv[i], method = method)
        } else { # Use maximum bound to calculate q
          approx.pb <- 2 # Create object for notification in output
          zx <- -(700-0.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i])+zcv[i]
          q[i] <- approx_puni(zx, zval[i], zcv[i])
        }
        ### Calculate test statistic and p-value
        if (method == "LNP")
        {
          L.pb <- sum(-log(q))
          pval.pb <- pgamma(L.pb, ksig, 1)
        } else if (method == "LN1MINP")
        {
          L.pb <- sum(-log(1-q))
          pval.pb <- exp(pgamma(L.pb, ksig, 1, lower.tail = FALSE, log.p = TRUE))
        } else if (method == "P")
        {
          L.pb <- (sum(q)-ksig*0.5)/sqrt(ksig/12)
          pval.pb <- pnorm(L.pb, lower.tail = FALSE)
        } 
      }
    }
  }
  
  return(list(L.pb = L.pb, pval.pb = pval.pb, ksig = ksig, approx.pb = max(approx.pb)))
}
