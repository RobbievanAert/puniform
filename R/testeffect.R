### Function for test of an effect p-uniform
testeffect <- function(yi, vi, zval, zcv, ksig, method, est) {

  if (method == "KS" | method == "AD")
  { # Return NA for test of an effect
    L.0 <- NA
    pval.0 <- NA
    approx.0.imp <- 0 # Create object for notification in output
  } else if (method == "ML")
    { # Likelihood-ratio test
      L.0 <- LR_test(d.alt = est, d.null = 0, yi = yi, vi = vi, zcv = zcv)
      pval.0 <- pchisq(L.0, df = 1, lower.tail = FALSE)
      approx.0.imp <- 0 # Create object for notification in output
    } else
    {
    q <- numeric(ksig) # Empty object for storing transformed p-values
    ### Loop for computing transformed p-values
    for(i in 1:length(zval))
    {
      if (zval[i] <= 38)
      { # Do not use approximation
        approx.0.imp <- 0 # Create object for notification in output
        pmarg <- exp(pnorm(zcv[i], 0, 1, lower.tail = FALSE, log.p = TRUE))
        ph1 <- exp(pnorm(zval[i], 0, 1, lower.tail = FALSE, log.p = TRUE))
        q[i] <- ph1/pmarg
      } else { # Use maximum bound to calculate q
        approx.0.imp <- 1 # Create object for notification in output
        zx <- sqrt(1400 + zcv[i]^2)
        q[i] <- approx_puni(0, zx, zcv[i], method = method)
      }
      ### Calculate test statistic and p-value
      if (method == "LNP")
      {
        L.0 <- sum(-log(q))
        pval.0 <- exp(pgamma(L.0, ksig, 1, lower.tail = FALSE, log.p = TRUE))
      } else if (method == "LN1MINP")
      {
        L.0 <- sum(-log(1 - q))
        pval.0 <- pgamma(L.0, ksig, 1)
      } else if (method == "P")
      {
        L.0 <- (sum(q)-ksig*0.5)/sqrt(ksig/12)
        pval.0 <- pnorm(L.0)
      } 
    }
  }

  return(data.frame(L.0 = L.0, pval.0 = pval.0, approx.0.imp = max(approx.0.imp)))
}