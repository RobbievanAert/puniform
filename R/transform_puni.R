#' transform
#' 
#' @keywords internal

##############################################################
##### FUNCTION FOR TRANSFORMING FISHER Z CORRELATIONS TO #####
##### RAW CORRELATIONS AND RE-MIRRORING IN CASE OF LEFT  #####
##### TAILED TESTS                                       #####
##############################################################

transform_puni <- function(res.fe, res.es, side, measure) {

  ### Tranform to raw correlations
  if(measure == "COR") {
    est <- (exp(2*res.es$est) - 1)/(exp(2*res.es$est) + 1)
    ci.lb <- (exp(2*res.es$ci.lb) - 1)/(exp(2*res.es$ci.lb) + 1)
    ci.ub <- (exp(2*res.es$ci.ub) - 1)/(exp(2*res.es$ci.ub) + 1)
    est.fe <- (exp(2*res.fe$est.fe) - 1)/(exp(2*res.fe$est.fe) + 1)
    ci.lb.fe <- (exp(2*res.fe$ci.lb.fe) - 1)/(exp(2*res.fe$ci.lb.fe) + 1)
    ci.ub.fe <- (exp(2*res.fe$ci.ub.fe) - 1)/(exp(2*res.fe$ci.ub.fe) + 1)
    se.fe <- (exp(2*res.fe$se.fe) - 1)/(exp(2*res.fe$se.fe) + 1)
    zval.fe <- res.fe$zval.fe

    ### Re-mirror effect sizes
    if(side == "left") {
      est <- est * -1
      tmp <- ci.ub
      ci.ub <- ci.lb * -1
      ci.lb <- tmp * -1
      est.fe <- est.fe * -1
      tmp <- ci.ub.fe
      ci.ub.fe <- ci.lb.fe * -1
      ci.lb.fe <- tmp * -1
      zval.fe <- res.fe$zval.fe * -1
    }
  } else if(measure != "COR" & side == "left") {
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
  } else {
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
