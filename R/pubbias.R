###################################################################
##### FUNCTION FOR PUBLICATION BIAS TEST P-UNIFORM            #####
###################################################################

pubbias <- function(es, alpha, method, est.fe) {

  # ### Conduct fixed-effect meta-analysis
  # wi <- 1/es$vi  # Weight per study
  # est.fe <- sum(es$yi * wi)/sum(wi)  # FE meta-analytical estimate
  # se.fe <- sqrt(1/sum(wi))  # Standard error of meta-analytical estimate
  # zval.fe <- est.fe/se.fe  # Z-value for test of no effect
  # if (zval.fe > 0) {
  #   pval.fe <- pnorm(zval.fe, lower.tail = FALSE) # Compute one-sided p-value
  # } else {
  #   pval.fe <- pnorm(zval.fe)
  # }
  # ci.lb.fe <- est.fe - qnorm(0.975) * se.fe # Lower bound CI meta-analytical estimate
  # ci.ub.fe <- est.fe + qnorm(0.975) * se.fe # Upper bound CI meta-analytical estimate
  # Qstat <- sum(wi * (es$yi - est.fe)^2) # Q-statistic
  # Qpval <- pchisq(Qstat, df = length(es$yi) - 1, lower.tail = FALSE) # p-value of Q-statistic

  ### Create a subset of significant studies
  sub <- subset(es, es$pval < alpha)
  yi <- sub$yi
  vi <- sub$vi
  zval <- sub$zval
  zcv <- sub$zcv
  ksig <- nrow(sub)

  ### Check if there are significant studies
  if (ksig == 0) {
    L.pb <- NA
    pval.pb <- NA
    approx.pb <- NA
  } else {
    if (method == "KS" | method == "AD") {
      ### If method is KS or AD return NAs for publication bias test
      L.pb <- NA
      pval.pb <- NA
      approx.pb <- 0 # Create object for notification in output
    } else { # Transform est.fe to zd for approximation
      zd <- est.fe/sqrt(vi)
      q <- numeric(ksig)  # Empty object for storing transformed p-values
      ### Loop for computing transformed p-values
      for (i in 1:length(yi)) {
        if (zcv[i] - zd[i] <= 38) { # Do not use approximation
          approx.pb <- 0 # Create object for notification in output
          pmarg <- exp(pnorm(zcv[i] * sqrt(vi[i]), est.fe, sqrt(vi[i]),
                             lower.tail = FALSE, log.p = TRUE))
          ph1 <- exp(pnorm(yi[i], est.fe, sqrt(vi[i]), lower.tail = FALSE,
                           log.p = TRUE))
          q[i] <- ph1/pmarg
        } else if (zd[i] > -(700-0.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i])+zcv[i]) {
          ### If statement above is TRUE: approximate P(Z>=zval) and P(Z>=zcv)
          approx.pb <- 1 # Create object for notification in output
          q[i] <- approx(zd[i], zval[i], zcv[i])
        } else { # Use maximum bound to calculate q
          approx.pb <- 2 # Create object for notification in output
          zx <- -(700-0.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i])+zcv[i]
          q[i] <- approx(zx, zval[i], zcv[i])
        }
        ### Calculate test statistic and p-value
        if (method == "LNP") {
          L.pb <- sum(-log(q))
          pval.pb <- pgamma(L.pb, ksig, 1)
        } else if (method == "LN1MINP") {
          L.pb <- sum(-log(1-q))
          pval.pb <- exp(pgamma(L.pb, ksig, 1, lower.tail = FALSE, log.p = TRUE))
        } else if (method == "P") {
          L.pb <- (sum(q)-ksig*0.5)/sqrt(ksig/12)
          pval.pb <- pnorm(L.pb, lower.tail = FALSE)
        }
      }
    }
  }

  return(list(data = data.frame(yi, vi, zval, zcv), L.pb = L.pb, pval.pb = pval.pb,
              ksig = ksig, approx.pb = max(approx.pb)))
}
