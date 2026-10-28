### Function for likelihood ratio test with puniform and method ML
LR_test <- function (d.alt, d.null, yi, vi, zcv)
{
  2*(ml(d = d.alt, yi = yi, vi = vi, zcv = zcv)-ml(d = d.null, yi = yi, vi = vi, zcv = zcv))
}
