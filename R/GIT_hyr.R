hyr <- function(es, m1i, m2i, mi, sd1i, sd2i, sdi, n1i, n2i, ni, tobs, measure, side) {
  
  if (measure == "M") {
    stat.hyr <- mi[2]/(sdi[2]/sqrt(ni[2]))
  } else if (measure == "MT") {
    stat.hyr <- tobs[2]
  } else if (measure == "MD") {
    s.pool <- sqrt(((n1i[2]-1)*sd1i[2]^2 + (n2i[2]-1)*sd2i[2]^2)/(n1i[2]+n2i[2]-2))
    stat.hyr <- (m1i[2]-m2i[2])/sqrt(s.pool^2*(1/n1i[2]+1/n2i[2]))
  } else if (measure == "MDT") {
    stat.hyr <- tobs[2]
  } else if (measure == "COR") {
    stat.hyr <- es$yi[2]/sqrt(es$vi[2]) 
  }
  
  est.hyr <- es$yi[2] 
  ci.lb.hyr <- est.hyr-qnorm(.975)*sqrt(es$vi[2]) 
  ci.ub.hyr <- est.hyr+qnorm(.975)*sqrt(es$vi[2]) 
  pval.hyr <- ifelse(es$pval[2] > 0.5, (1-es$pval[2])*2, es$pval[2]*2) 
  pval.o <- ifelse(es$pval[1] > 0.5, (1-es$pval[1])*2, es$pval[1]*2) 
  
  if (measure == "COR") { 
    est.hyr <- (exp(2*es$yi[2]) - 1)/(exp(2*es$yi[2]) + 1)
    ci.lb.hyr <- (exp(2*ci.lb.hyr) - 1)/(exp(2*ci.lb.hyr) + 1)
    ci.ub.hyr <- (exp(2*ci.ub.hyr) - 1)/(exp(2*ci.ub.hyr) + 1)
  }
  
  if (side == "left") { 
    est.hyr <- est.hyr*-1
    tmp <- ci.ub.hyr
    ci.ub.hyr <- ci.lb.hyr*-1
    ci.lb.hyr <- tmp*-1    
  }
  
  return(data.frame(est.hyr = est.hyr, ci.lb.hyr = ci.lb.hyr, ci.ub.hyr = ci.ub.hyr, stat.hyr = stat.hyr, pval.hyr = pval.hyr, pval.o = pval.o))
}