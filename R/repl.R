### Function for estimation based on only replication
repl <- function(es, m1i, m2i, mi, sd1i, sd2i, sdi, n1i, n2i, ni, tobs, measure, side) {
  
  if (measure == "M") {
    stat.repl <- mi[2]/(sdi[2]/sqrt(ni[2]))
  } else if (measure == "MT") {
    stat.repl <- tobs[2]
  } else if (measure == "MD") {
    s.pool <- sqrt(((n1i[2]-1)*sd1i[2]^2 + (n2i[2]-1)*sd2i[2]^2)/(n1i[2]+n2i[2]-2))
    stat.repl <- (m1i[2]-m2i[2])/sqrt(s.pool^2*(1/n1i[2]+1/n2i[2]))
  } else if (measure == "MDT") {
    stat.repl <- tobs[2]
  } else if (measure == "COR") {
    stat.repl <- es$yi[2]/sqrt(es$vi[2]) # z-value for test of no effect
  }
  
  est.repl <- es$yi[2] # Effect size estimate replication
  ci.lb.repl <- est.repl-qnorm(.975)*sqrt(es$vi[2]) # Lower bound of CI replication
  ci.ub.repl <- est.repl+qnorm(.975)*sqrt(es$vi[2]) # Upper bound of CI replication
  pval.repl <- pnorm(es$zval[2], lower.tail = FALSE) # Compute p-value replication
  pval.repl <- ifelse(pval.repl > 0.5, (1-pval.repl)*2, pval.repl*2) # Compute two-tailed p-value of replication
  pval.o <- ifelse(es$pval[1] > 0.5, (1-es$pval[1])*2, es$pval[1]*2) # Compute two-tailed p-value of original study
  
  if (measure == "COR") { # Back transform Fisher z to correlation
    est.repl <- (exp(2*es$yi[2]) - 1)/(exp(2*es$yi[2]) + 1)
    ci.lb.repl <- (exp(2*ci.lb.repl) - 1)/(exp(2*ci.lb.repl) + 1)
    ci.ub.repl <- (exp(2*ci.ub.repl) - 1)/(exp(2*ci.ub.repl) + 1)
  }
  
  if (side == "left") { # Re-mirror estimates
    est.repl <- est.repl*-1
    tmp <- ci.ub.repl
    ci.ub.repl <- ci.lb.repl*-1
    ci.lb.repl <- tmp*-1    
  }
  
  return(data.frame(est.repl = est.repl, se.repl = sqrt(es$vi[2]), ci.lb.repl = ci.lb.repl, 
                    ci.ub.repl = ci.ub.repl, stat.repl = stat.repl, pval.repl = pval.repl, 
                    pval.o = pval.o))
}