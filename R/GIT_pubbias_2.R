pubbias <- function(es, alpha, method) {

  res <- rma(yi = es$yi, vi = es$vi, method = "FE")
  est.fe <- res$b[1]

  sub <- subset(es, es$pval < alpha)
  yi <- sub$yi
  vi <- sub$vi
  zval <- sub$zval
  zcv <- qnorm(alpha, lower.tail = FALSE)
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
      q <- numeric()
      for(i in 1:length(yi)) {
        if(zcv - zd[i] <= 38) {
          approx.pb <- 0
          pmarg <- exp(pnorm(zcv*sqrt(vi[i]), est.fe, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
          ph1 <- exp(pnorm(yi[i], est.fe, sqrt(vi[i]), lower.tail = FALSE, log.p = TRUE))
          q[i] <- ph1/pmarg
        } else if(zd[i] > -(700-.5*(zval[i]-zcv)^2)/(zval[i]-zcv) + zcv) {
          approx.pb <- 1
          q[i] <- approx(zd[i], zval[i], zcv)
        } else {
          approx.pb <- 2
          zx <- -(700-.5*(zval[i]-zcv)^2)/(zval[i]-zcv) + zcv
          q[i] <- approx(zx, zval[i], zcv)
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
  return(list(data = data.frame(yi, vi, zval), L.pb = L.pb, pval.pb = pval.pb, ksig = ksig, res = res, approx.pb = max(approx.pb)))
}
