esest <- function(yi, vi, zval, zcv, ksig, method) {

  bo <- bounds(yi = yi, vi = vi, zval = zval, zcv = zcv, ext = FALSE)
  bo.ext <- bounds(yi = yi, vi = vi, zval = zval, zcv = zcv, ext = TRUE)

  pdist <- function(d, yi, vi, zval, zcv, ksig, val, method, cv.P) {
    zd <- d/sqrt(vi)
    q <- numeric(ksig)
    for(i in 1:length(yi)) {
      if(zcv[i] - zd[i] <= 38) {
        pmarg <- exp(pnorm(zcv[i]*sqrt(vi[i]), mean = d, sd = sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
        ph1 <- exp(pnorm(yi[i], mean = d, sd = sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
        q[i] <- ph1/pmarg
      } else {
        q[i] <- approx(zd[i], zval[i], zcv[i])
      }
    }
    if(val == "est") { tr.q <<- q }
    if(method == "LNP") { stat <- sum(-log(q)) }
    else if(method == "LN1MINP") { stat <- sum(-log(1 - q)) }
    else if(method == "P") { stat <- sum(q) }
    else if(method == "KS" & val == "est") { out <- ks.test(x = q, y = punif)$statistic }
    else if(method == "AD" & val == "est") { out <- as.numeric(ADGofTest::ad.test(x = q, distr.fun = punif)$statistic) }
    if(val == "est" & (method == "LNP" | method == "LN1MINP" | method == "P")) { out <- stat - ksig }
    else if(val == "ci.lb" & (method == "LNP" | method == "LN1MINP")) { out <- stat - qgamma(.975, ksig, 1) }
    else if(val == "ci.ub" & (method == "LNP" | method == "LN1MINP")) { out <- stat - qgamma(.025, ksig, 1) }
    else if(val == "ci.ub" & method == "P") { out <- stat - cv.P }
    else if(val == "ci.lb" & method == "P") { out <- stat - cv.P }
    out
  }

  if(method == "KS" | method == "AD") {
    ci.lb <- NA
    ci.ub <- NA
    approx.ci.lb <- 0
    approx.ci.ub <- 0
    if(method == "KS") {
      est.bo <- try(bisect(func = pdist, lo = bo[1], hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig/2, val = "est", method = "P"), silent = TRUE)
      if(class(est.bo) == "try-error") {
        est.bo <- bisect(func = pdist, lo = bo.ext[1], hi = bo.ext[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig/2, val = "est", method = "P")
      }
      est <- try(as.numeric(optimize(pdist, interval = c(est.bo-1.5, est.bo+1.5), yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "est", method = method)$minimum), silent = TRUE)
      if(class(est) == "try-error") { est <- NA }
    }
    if(method == "AD") {
      est.AD <- nlm(pdist, p = 0, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "est", method = method)
      if(est.AD$gradient < 0.1) { est <- est.AD$estimate
      } else { est <- NA }
    }
    if(any(is.na(est) == FALSE & any(zcv - est/sqrt(vi) > 38))) { approx.est <- 1
    } else { approx.est <- 0 }

  } else {
    if(method == "P") {
      ksig.est <- ksig/2
    } else { ksig.est <- ksig }
    approx.est <- 0
    est <- try(bisect(func = pdist, lo = bo[1], hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig.est, val = "est", method = method), silent = TRUE)
    if(class(est) == "try-error") {
      approx.est <- 1
      est <- try(bisect(func = pdist, lo = bo.ext[1], hi = bo.ext[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig.est, val = "est", method = method), silent = TRUE)
      if(class(est) == "try-error") { est <- NA }
    }

    if(method == "P") {
      approx.ci.lb <- 0
      ci.lb <- try(bisect(func = pdist, lo = bo[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method, cv.P = get.cv.P(ksig)), silent = TRUE)
      if(class(ci.lb) == "try-error") {
        approx.ci.lb <- 1
        ci.lb <- try(bisect(func = pdist, lo = bo.ext[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method, cv.P = get.cv.P(ksig)), silent = TRUE)
        if(class(ci.lb) == "try-error") { ci.lb <- NA }
      }
    } else if(method == "LN1MINP") {
      approx.ci.lb <- 0
      ci.lb <- try(bisect(func = pdist, lo = bo[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.ub", method = method), silent = TRUE)
      if(class(ci.lb) == "try-error") {
        approx.ci.lb <- 1
        ci.lb <- try(bisect(func = pdist, lo = bo.ext[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.ub", method = method), silent = TRUE)
        if(class(ci.lb) == "try-error") { ci.lb <- NA }
      }
    } else {
      approx.ci.lb <- 0
      ci.lb <- try(bisect(func = pdist, lo = bo[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method), silent = TRUE)
      if(class(ci.lb) == "try-error") {
        approx.ci.lb <- 1
        ci.lb <- try(bisect(func = pdist, lo = bo.ext[1], hi = est, yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method), silent = TRUE)
        if(class(ci.lb) == "try-error") { ci.lb <- NA }
      }
    }

    if(method == "P") {
      ci.ub <- try(bisect(func = pdist, lo = est, hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.ub", method = method, cv.P = ksig-get.cv.P(ksig)), silent = TRUE)
      if(class(ci.ub) == "try-error") { ci.ub <- NA }
    } else if(method == "LN1MINP") {
      if(is.na(est) == TRUE) {
        ci.ub <- try(bisect(func = pdist, lo = bo[1], hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method), silent = TRUE)
      } else {
        ci.ub <- try(bisect(func = pdist, lo = est, hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.lb", method = method), silent = TRUE)
      }
    } else if(is.na(est) == TRUE) {
      ci.ub <- try(bisect(func = pdist, lo = bo[1], hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.ub", method = method), silent = TRUE)
    } else {
      ci.ub <- try(bisect(func = pdist, lo = est, hi = bo[2], yi = yi, vi = vi, zval = zval, zcv = zcv, ksig = ksig, val = "ci.ub", method = method), silent = TRUE)
    }
    if(class(ci.ub) == "try-error") { ci.ub <- NA }
  }

  return(list(est = est, ci.lb = ci.lb, ci.ub = ci.ub, approx.est = max(approx.est), approx.ci.lb = max(approx.ci.lb), ext.lb = bo.ext[1], tr.q = tr.q))
}
