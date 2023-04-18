### Function for computing log-likelihood in case of p-uniform* with moderators
ml_mods <- function(par, es, n_bs, par_fixed, mods, transf)
{
  ### Add fixed parameters to vector of optimized parameters to get the correct
  # log-likelihood
  par_est <- par_fixed
  par_est[which(is.na(par_est))] <- par
  
  # cat("b1 = ", par_est[1], "b2 = ", par_est[2], "tau2 = ", par_est[3], fill = TRUE)
  
  bs <- par_est[1:n_bs] # Get values for bs
  tau2 <- par_est[n_bs+1] # Get value for tau
  
  ### If the log of tau2 is optimized, take the exponent of tau2  
  if (transf == TRUE) tau2 <- exp(tau2)
  
  yi <- es$yi
  vi <- es$vi
  ycv <- es$ycv
  
  ### Evaluate the regression equation at the estimated bs
  X <- model.matrix(mods, data = es)
  
  ### Compute the means
  M <- X %*% bs
  
  ### Compute the log-likelihood of the truncated densities
  q <- mapply(function(M, yi, vi, ycv, tau2)
  {
    ifelse(yi > ycv,
           dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) -
             pnorm(ycv, mean = M, sd = sqrt(vi+tau2), lower.tail = FALSE, log.p = TRUE),
           dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) - 
             pnorm(ycv, mean = M, sd = sqrt(vi+tau2), log.p = TRUE))
           
           
  }, M = M, yi = yi, vi = vi, ycv = ycv, MoreArgs = list(tau2 = tau2))
  
  return(-sum(q))
}