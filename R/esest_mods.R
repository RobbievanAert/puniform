esest_mods <- function(yi, vi, xi, ycv, con)
{
  par <- con$par
  est <- nlminb(start = par, objective = ml_mods, yi = yi, vi = vi, xi = xi, 
                ycv = ycv)
  
  # return(data.frame(b0 = est$par[1], b1 = est$par[2], tau = est$par[3]))
  return(est)
}