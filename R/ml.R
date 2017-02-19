### Log likelihood function for puniform with ML method
ml <- function(d, yi, vi, zcv)
{

  q <- mapply(function(d, yi, vi, zcv)
  {
    ### In case of extreme conditional density, use approximation
    ifelse(yi/sqrt(vi)-d/sqrt(vi) < 36,
           dnorm(yi, mean = d, sd = sqrt(vi))/
             exp(pnorm(zcv * sqrt(vi), mean = d, sd = sqrt(vi), lower.tail = FALSE, log.p = TRUE)),
           approx(zd = d/sqrt(vi), zval = yi/sqrt(vi), zcv = zcv, method = "ML"))
  }, yi = yi, zcv = zcv, vi = vi, MoreArgs = list(d = d))

  log(prod(q))

}
