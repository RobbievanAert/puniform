### Function for computing log-likelihood in case of hybrid method with moderators
ml_hy <- function(par, es, mods, n_bs, par_fixed, transf, verbose)
{
  yi <- es$yi
  vi <- es$vi
  ycv <- es$ycv
  conventional <- es$conventional
  
  ### Add fixed parameters to vector of optimized parameters to get the correct
  # log-likelihood
  par_est <- par_fixed
  par_est[which(is.na(par_est))] <- par
  
  bs <- par_est[1:n_bs] # Get values for bs
  tau2 <- par_est[n_bs+1] # Get value for tau
  
  ### If the log of tau2 is optimized, take the exponent of tau2  
  if (transf == TRUE) tau2 <- exp(tau2)
  
  ### If progress of model fitting needs to be shown
  if (verbose == TRUE) cat("bs = ", bs, "tau2 = ", tau2 , fill = TRUE)
  
  ### Evaluate the regression equation at the estimated bs
  X <- model.matrix(mods, data = es)
  
  ### Compute the means
  M <- X %*% bs
  
  ### Compute the log-likelihood of the truncated densities
  q <- mapply(function(M, yi, vi, ycv, conventional, tau2)
  {
    ifelse(conventional == 1, 
           ifelse(yi > ycv,
                  dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) -
                    pnorm(ycv, mean = M, sd = sqrt(vi+tau2), lower.tail = FALSE, log.p = TRUE),
                  dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) - 
                    pnorm(ycv, mean = M, sd = sqrt(vi+tau2), log.p = TRUE)),
           dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE))
    
  }, M = M, yi = yi, vi = vi, ycv = ycv, conventional = conventional, MoreArgs = list(tau2 = tau2))
  
  return(-sum(q))
}