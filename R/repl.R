### Function for estimation based on only replication
repl <- function(es, m1i, m2i, mi, ri, sd1i, sd2i, sdi, n1i, n2i, ni, tobs, 
                 measure, side) 
{
  
  ### Conduct FE meta-analysis based on replications
  sub <- subset(es, es$ori == 0)
  tmp <- fe_ma(yi = sub$yi, vi = sub$vi)
  est.repl <- tmp$est.fe
  se.repl <- tmp$se.fe
  ci.lb.repl <- tmp$ci.lb.fe
  ci.ub.repl <- tmp$ci.ub.fe
  
  if (!missing("mi") | !missing("tobs") | !missing("m1i") | !missing("ri"))
  { # If only two studies are supplied as input
    
    if (measure == "M") 
    {
      stat.repl <- mi[2]/(sdi[2]/sqrt(ni[2]))
    } else if (measure == "MT") 
    {
      stat.repl <- tobs[2]
    } else if (measure == "MD") 
    {
      s.pool <- sqrt(((n1i[2]-1)*sd1i[2]^2 + (n2i[2]-1)*sd2i[2]^2)/(n1i[2]+n2i[2]-2))
      stat.repl <- (m1i[2]-m2i[2])/sqrt(s.pool^2*(1/n1i[2]+1/n2i[2]))
    } else if (measure == "MDT") 
    {
      stat.repl <- tobs[2]
    } else if (measure == "COR") 
    {
      stat.repl <- es$yi[2]/sqrt(es$vi[2]) # z-value for test of no effect
    }
    
    pval.repl <- pnorm(stat.repl, lower.tail = FALSE) # One-tailed p-value
    pval.repl <- ifelse(pval.repl > 0.5, (1-pval.repl)*2, pval.repl*2) # Two-tailed p-value
    
  } else
  {
    stat.repl <- tmp$zval.fe
    pval.repl <- tmp$pval.fe
  }
  
  ### Conduct FE meta-analysis based on original studies
  sub <- subset(es, es$ori == 1)
  tmp <- fe_ma(yi = sub$yi, vi = sub$vi)
  pval.o <- tmp$pval.fe
  
  if (measure == "COR") 
  { # Back transform Fisher z to correlation
    est.repl <- (exp(2*est.repl) - 1)/(exp(2*est.repl) + 1)
    ci.lb.repl <- (exp(2*ci.lb.repl) - 1)/(exp(2*ci.lb.repl) + 1)
    ci.ub.repl <- (exp(2*ci.ub.repl) - 1)/(exp(2*ci.ub.repl) + 1)
  }
  
  if (side == "left") 
  { # Re-mirror estimates
    est.repl <- est.repl*-1
    tmp <- ci.ub.repl
    ci.ub.repl <- ci.lb.repl*-1
    ci.lb.repl <- tmp*-1    
  }
  
  return(data.frame(est.repl = est.repl, se.repl = se.repl, ci.lb.repl = ci.lb.repl, 
                    ci.ub.repl = ci.ub.repl, stat.repl = stat.repl, pval.repl = pval.repl, 
                    pval.o = pval.o))
}