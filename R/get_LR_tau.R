### Function for estimating profile likelihood CI for tau
get_LR_tau <- function(prof.tau, yi, vi, est, tau.est, ycv)
{
  (2*(ml_star(par = c(est, tau.est), yi = yi, vi = vi, ycv = ycv)-
        ml_star(par = c(est, prof.tau), yi = yi, vi = vi, ycv = ycv)))-qchisq(.95, df = 1)
}