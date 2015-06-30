pubbias <- function(es, alpha, method) {

  wi <- 1/es$vi
  est.fe <- sum(es$yi*wi)/sum(wi)
  se.fe <- sqrt(1/sum(wi))
  zval.fe <- est.fe/se.fe
  if (zval.fe > 0) {
    pval.fe <- pnorm(zval.fe, lower.tail = FALSE)
  } else { pval.fe <- pnorm(zval.fe) }
  ci.lb.fe <- est.fe-qnorm(0.975)*se.fe
  ci.ub.fe <- est.fe+qnorm(0.975)*se.fe
  Qstat <- sum(wi*(es$yi-est.fe)^2)
  Qpval <- pchisq(Qstat, df = length(es$yi)-1, lower.tail = FALSE)

  sub <- subset(es, es$pval < alpha)
  yi <- sub$yi
  vi <- sub$vi
  zval <- sub$zval
  zcv <- sub$zcv
  ksig <- nrow(sub)

  if(ksig == 0) {
    L.pb <- NA
    pval.pb <- NA
    approx.pb <- NA
  } else {
    if(method == "KS" | method == "AD") {
      L.pb <- NA
      pval.pb <- NA
      approx.pb <- 0
    } else {
      zd <- est.fe/sqrt(vi)
      q <- numeric(ksig)
      for(i in 1:length(yi)) {
        if(zcv[i] - zd[i] <= 38) {
          approx.pb <- 0
          pmarg <- exp(pnorm(zcv[i]*sqrt(vi[i]), est.fe, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
          ph1 <- exp(pnorm(yi[i], est.fe, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
          q[i] <- ph1/pmarg
        } else if(zd[i] > -(700-.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i]) + zcv[i]) {
          approx.pb <- 1
          q[i] <- approx(zd[i], zval[i], zcv[i])
        } else {
          approx.pb <- 2
          zx <- -(700-.5*(zval[i]-zcv[i])^2)/(zval[i]-zcv[i]) + zcv[i]
          q[i] <- approx(zx, zval[i], zcv[i])
        }
        if(method == "LNP") {
          L.pb <- sum(-log(q))
          pval.pb <- pgamma(L.pb, ksig, 1)
        } else if(method == "LN1MINP") {
          L.pb <- sum(-log(1 - q))
          pval.pb <- exp(pgamma(L.pb, ksig, 1, lower.tail = FALSE, log.p = TRUE))
        } else if(method == "P") {
          L.pb <- (sum(q)-ksig*0.5)/sqrt(ksig/12)
          pval.pb <- pnorm(L.pb, lower.tail = FALSE)
        }
      }
    }
  }
  return(list(data = data.frame(yi, vi, zval, zcv), L.pb = L.pb, pval.pb = pval.pb, ksig = ksig, est.fe = est.fe,
              se.fe = se.fe, zval.fe = zval.fe, pval.fe = pval.fe, ci.lb.fe = ci.lb.fe, ci.ub.fe = ci.ub.fe,
              Qstat = Qstat, Qpval = Qpval, approx.pb = max(approx.pb)))
}
