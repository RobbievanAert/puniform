### Function for re-mirroring in case of left-tailed tests with p-uniform*
transform_nsig <- function(res.es, side) 
{
  
  if(side == "left") 
  {
    ### Re-mirror effect sizes
    est <- res.es$est * -1
    tmp <- res.es$ub
    ub <- res.es$lb * -1
    lb <- tmp * -1
  } else {
    est <- res.es$est
    lb <- res.es$lb
    ub <- res.es$ub
  }
  
  return(data.frame(est = est, lb = lb, ub = ub))
}