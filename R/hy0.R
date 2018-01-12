### Function for applying hybrid^0 method
hy0 <- function(es, res1, alpha) {
  ### Compute average conditional probability
  ave <- mean(c(pnorm(es$zval[1], lower.tail = FALSE)/alpha, pnorm(es$zval[2], lower.tail = FALSE))) 
  if (ave > 0.5) { # If average conditional probability is larger than 0.5
    res1$est.hy <- 0
    res1$ci.lb.hy <- NA
    res1$ci.ub.hy <- NA
    res1$x.hy <- 1
    res1$pval.hy <- 0.5
  }   
  return(data.frame(est.hy0 = res1$est, ci.lb.hy0 = res1$ci.lb, ci.ub.hy0 = res1$ci.ub, 
                    x.hy0 = res1$x, pval.hy0 = res1$pval, ave = ave))
}