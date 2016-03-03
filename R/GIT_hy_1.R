hy <- function(es, measure, side, alpha) {
    
  bo <- bounds.hy(yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, ext = FALSE)
    
  est.hy <- try(bisect(func = pdist.hy, lo = bo[1], hi = bo[2], yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, k = 1, val = "est"), silent = TRUE)
  if(class(est.hy) == "try-error") { 
    
    bo.ext <- bounds.hy(yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, ext = TRUE)
    est.hy <- try(bisect(func = pdist.hy, lo = bo.ext[1], hi = bo.ext[2], yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, k = 1, val = "est"), silent = TRUE)
    
    if(class(est.hy) == "try-error") { est.hy <- NA }
  }  
  
  ci.lb.hy <- try(bisect(func = pdist.hy, lo = bo[1], hi = est.hy, yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, k = nrow(es), val = "ci.lb", cv.P = 0.2236068), silent = TRUE)
  if(class(ci.lb.hy) == "try-error") {
    
    bo.ext <- bounds.hy(yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, ext = TRUE)
    ci.lb.hy <- try(bisect(func = pdist.hy, lo = bo.ext[1], hi = est.hy, yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, k = nrow(es), val = "ci.lb", cv.P = 0.2236068), silent = TRUE)
    if(class(ci.lb.hy) == "try-error") { ci.lb.hy <- NA }
  }    
  
  ci.ub.hy <- try(bisect(func = pdist.hy, lo = est.hy, hi = bo[2], yi = es$yi, vi = es$vi, zval = es$zval, zcv = es$zcv, k = nrow(es), val = "ci.ub", cv.P = 2-0.2236068), silent = TRUE)
  if(class(ci.ub.hy) == "try-error") { ci.ub.hy <- NA }
  
  if (measure == "COR") { 
    est.hy <- (exp(2*est.hy) - 1)/(exp(2*est.hy) + 1)
    ci.lb.hy <- (exp(2*ci.lb.hy) - 1)/(exp(2*ci.lb.hy) + 1)
    ci.ub.hy <- (exp(2*ci.ub.hy) - 1)/(exp(2*ci.ub.hy) + 1)
  }
  
  if (side == "left") { 
    est.hy <- est.hy*-1
    tmp <- ci.ub.hy
    ci.ub.hy <- ci.lb.hy*-1
    ci.lb.hy <- tmp*-1
  }  
  
  q <- numeric(nrow(es))    
  pmarg <- exp(pnorm(es$zcv[1], mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE))
  ph1 <- exp(pnorm(es$zval[1], mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE))
  q[1] <- ph1/pmarg
  q[2] <- exp(pnorm(es$zval[2], mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE))
  x.hy <- ifelse(sum(q) < 1, sum(q), 2-sum(q)) 
  pval.hy <- ifelse(sum(q) < 1, 0.5*sum(q)^2, -0.5*sum(q)^2+2*sum(q)-1)   
  pval.hy <- ifelse(pval.hy > 0.5, (1-pval.hy)*2, pval.hy*2) 
  
  return(data.frame(est.hy = est.hy, ci.lb.hy = ci.lb.hy, ci.ub.hy = ci.ub.hy, x.hy = x.hy, pval.hy = pval.hy))
}