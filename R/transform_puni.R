### Function for re-mirroring in case of left-tailed tests with p-uniform
transform_puni <- function(res.fe, res.es, side) 
{
  
  if(side == "left") 
  {
    ### Re-mirror effect sizes
    est <- res.es$est * -1
    tmp <- res.es$ci.ub
    ci.ub <- res.es$ci.lb * -1
    ci.lb <- tmp * -1
    est.fe <- res.fe$est.fe * -1
    tmp <- res.fe$ci.ub.fe
    ci.ub.fe <- res.fe$ci.lb.fe * -1
    ci.lb.fe <- tmp * -1
    se.fe <- res.fe$se.fe
    zval.fe <- res.fe$zval.fe * -1
  } else 
  {
    est <- res.es$est
    ci.lb <- res.es$ci.lb
    ci.ub <- res.es$ci.ub
    est.fe <- res.fe$est.fe
    ci.lb.fe <- res.fe$ci.lb.fe
    ci.ub.fe <- res.fe$ci.ub.fe
    se.fe <- res.fe$se.fe
    zval.fe <- res.fe$zval.fe
  }
  
  return(data.frame(est = est, ci.lb = ci.lb, ci.ub = ci.ub, est.fe = est.fe,
                    zval.fe = zval.fe, ci.lb.fe = ci.lb.fe, ci.ub.fe = ci.ub.fe,
                    se.fe = se.fe))
}