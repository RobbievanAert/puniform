transform <- function(res1, res3, side, measure) {

  if(measure == "COR") {
    est <- (exp(2*res3$est) - 1)/(exp(2*res3$est) + 1)
    ci.lb <- (exp(2*res3$ci.lb) - 1)/(exp(2*res3$ci.lb) + 1)
    ci.ub <- (exp(2*res3$ci.ub) - 1)/(exp(2*res3$ci.ub) + 1)

    est.fe <- (exp(2*res1$res$b[1]) - 1)/(exp(2*res1$res$b[1]) + 1)
    ci.lb.fe <- (exp(2*res1$res$ci.lb) - 1)/(exp(2*res1$res$ci.lb) + 1)
    ci.ub.fe <- (exp(2*res1$res$ci.ub) - 1)/(exp(2*res1$res$ci.ub) + 1)
    se.fe <- (exp(2*res1$res$se) - 1)/(exp(2*res1$res$se) + 1)

    if(side == "left") {
      est <- est * -1
      tmp <- ci.ub
      ci.ub <- ci.lb * -1
      ci.lb <- tmp * -1
      est.fe <- est.fe * -1
      tmp <- ci.ub.fe
      ci.ub.fe <- ci.lb.fe * -1
      ci.lb.fe <- tmp * -1
    }
  } else if(measure != "COR" & side == "left") {
    est <- res3$est * -1
    tmp <- res3$ci.ub
    ci.ub <- res3$ci.lb * -1
    ci.lb <- tmp * -1
    est.fe <- res1$res$b[1] * -1
    tmp <- res1$res$ci.ub
    ci.ub.fe <- res1$res$ci.lb * -1
    ci.lb.fe <- tmp * -1
    se.fe <- res1$res$se
  } else {
    est <- res3$est
    ci.lb <- res3$ci.lb
    ci.ub <- res3$ci.ub
    est.fe <- res1$res$b[1]
    ci.lb.fe <- res1$res$ci.lb
    ci.ub.fe <- res1$res$ci.ub
    se.fe <- res1$res$se
  }

  return(data.frame(est = est, ci.lb = ci.lb, ci.ub = ci.ub, est.fe = est.fe, ci.lb.fe = ci.lb.fe, ci.ub.fe = ci.ub.fe, se.fe = se.fe))
}
