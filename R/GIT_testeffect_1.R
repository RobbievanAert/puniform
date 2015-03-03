testeffect <- function(zval, ksig, alpha, method) {

  zcv <- qnorm(alpha, lower.tail = FALSE)

  if(method == "KS" | method == "AD") {
    L.0 <- NA
    pval.0 <- NA
    approx.0.imp <- 0
  } else {
    q <- numeric()
    for(i in 1:length(zval)) {
      if(zval[i] <= 38) {
        approx.0.imp <- 0
        ph1 <- exp(pnorm(zval[i], 0, 1, lower.tail = FALSE, log.p = TRUE))
        q[i] <- ph1/alpha
      } else {
        approx.0.imp <- 1
        zx <- sqrt(1400 + zcv^2)
        q[i] <- approx(0, zx, zcv)
      }
      if(method == "LNP") {
        L.0 <- sum(-log(q))
        pval.0 <- exp(pgamma(L.0, ksig, 1, lower.tail = FALSE, log.p = TRUE))
      } else if(method == "LN1MINP") {
        L.0 <- sum(-log(1 - q))
        pval.0 <- pgamma(L.0, ksig, 1)
      } else if(method == "P") {
        L.0 <- (sum(q)-ksig*0.5)/sqrt(ksig/12)
        pval.0 <- pnorm(L.0)
      }
    }
  }
  return(data.frame(L.0 = L.0, pval.0 = pval.0, approx.0.imp = max(approx.0.imp)))
}
