### Function for estimating profile likelihood CI for tau
get_LR_tau <- function(prof.tau, yi, vi, est, tau.est, ycv)
{
  (2*(ml_tau(tau.est, est, yi, vi, ycv)-
        ml_tau(prof.tau, est, yi, vi, ycv)))-qchisq(.95, df = 1)
}