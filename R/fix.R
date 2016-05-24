### Function for fixed-effect meta-analysis
fix <- function(yi, vi, measure, side) {
  wi <- 1/vi # Weight per study
  est.fe <- sum(yi*wi)/sum(wi) # FE meta-analytic estimate
  se.fe <- sqrt(1/sum(wi)) # Standard error of meta-analytic estimate
  ci.lb.fe <- est.fe-qnorm(0.975)*se.fe # Lower bound CI meta-analytical estimate
  ci.ub.fe <- est.fe+qnorm(0.975)*se.fe # Upper bound CI meta-analytical estimate
  zval.fe <- est.fe/se.fe # Z-value for test of no effect
  pval.fe <- pnorm(zval.fe, lower.tail = FALSE) # Compute one-sided p-value
  pval.fe <- ifelse(pval.fe > 0.5, (1-pval.fe)*2, pval.fe*2) # Compute two-tailed p-value
  
  if (measure == "COR") { # Back transform Fisher z to correlation
    est.fe <- (exp(2*est.fe) - 1)/(exp(2*est.fe) + 1)
    se.fe <- (exp(2*se.fe) - 1)/(exp(2*se.fe) + 1)
    ci.lb.fe <- (exp(2*ci.lb.fe) - 1)/(exp(2*ci.lb.fe) + 1)
    ci.ub.fe <- (exp(2*ci.ub.fe) - 1)/(exp(2*ci.ub.fe) + 1)    
  }
  
  if (side == "left") { # Re-mirror estimates
    est.fe <- est.fe*-1
    tmp <- ci.ub.fe
    ci.ub.fe <- ci.lb.fe*-1
    ci.lb.fe <- tmp*-1
    zval.fe <- zval.fe*-1
  }
  
  return(data.frame(est.fe = est.fe, se.fe = se.fe, ci.lb.fe = ci.lb.fe, 
                    ci.ub.fe = ci.ub.fe, zval.fe = zval.fe, pval.fe = pval.fe))
}