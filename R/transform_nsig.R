### Function for re-mirroring in case of left-tailed tests with p-uniform*
transform_nsig <- function(res.es, side) 
{
  
  if(side == "left") 
  {
    ### Re-mirror effect sizes
    est <- res.es$est * -1
    tmp <- res.es$ci.ub
    ci.ub <- res.es$ci.lb * -1
    ci.lb <- tmp * -1
  } else {
    est <- res.es$est
    ci.lb <- res.es$ci.lb
    ci.ub <- res.es$ci.ub
  }
  
  return(data.frame(est = est, ci.lb = ci.lb, ci.ub = ci.ub))
}