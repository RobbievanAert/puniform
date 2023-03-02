### Function for estimating profile likelihood CI for d
get_LR_est <- function(prof.est, yi, vi, est, tau.est, ycv)
{
  (2*(ml_star(par = c(est, tau.est), yi = yi, vi = vi, ycv = ycv)-
        ml_star(par = c(prof.est, tau.est), yi = yi, vi = vi, ycv = ycv)))-qchisq(.95, df = 1)
}