### Function for computing log-likelihood in case of p-uniform* with moderators
ml_mods <- function(par, es, mods, n_bs)
{
  bs <- par[1:n_bs] # Get values for bs
  tau <- par[n_bs+1] # Get value for tau
  
  yi <- es$yi
  vi <- es$vi
  ycv <- es$ycv
  
  ### Evaluate the regression equation at the estimated bs
  X <- model.matrix(mods, data = es)
  M <- X %*% bs
  
  ### Compute the log-likelihood of the truncated densities
  q <- mapply(function(M, yi, vi, ycv, tau)
  {
    ifelse(yi > ycv,
           dnorm(yi, mean = M, sd = sqrt(vi+tau^2), log = TRUE) -
             pnorm(ycv, mean = M, sd = sqrt(vi+tau^2), lower.tail = FALSE, log.p = TRUE),
           dnorm(yi, mean = M, sd = sqrt(vi+tau^2), log = TRUE) - 
             pnorm(ycv, mean = M, sd = sqrt(vi+tau^2), log.p = TRUE))
  }, M = M, yi = yi, vi = vi, ycv = ycv, MoreArgs = list(tau = tau))
  
  return(-sum(q))
}