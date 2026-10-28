### Log likelihood function for puniform with ML method
ml <- function(d, yi, vi, zcv)
{

  q <- mapply(function(d, yi, vi, zcv)
  {
    ### In case of extreme conditional density, use approximation
    ifelse(yi/sqrt(vi)-d/sqrt(vi) < 36,
           dnorm(yi/sqrt(vi), mean = d/sqrt(vi), sd = 1)/
             exp(pnorm(zcv, mean = d/sqrt(vi), sd = 1, lower.tail = FALSE, log.p = TRUE)),
           approx_puni(zd = d/sqrt(vi), zval = yi/sqrt(vi), zcv = zcv, method = "ML"))
  }, yi = yi, zcv = zcv, vi = vi, MoreArgs = list(d = d))

  log(prod(q))

}
