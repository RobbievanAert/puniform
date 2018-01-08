### Function for estimating profile likelihood CI for d
get_LR_est <- function(prof.est, yi, vi, est, tau.est, ycv)
{
  (2*(ml_est(est, tau.est, yi, vi, ycv)-
        ml_est(prof.est, tau.est, yi, vi, ycv)))-qchisq(.95, df = 1)
}