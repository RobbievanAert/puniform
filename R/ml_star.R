# ### Function for computing log likelihood of p-uniform*
# ml_star <- function(par, yi, vi, ycv, verbose = FALSE) 
# {
#   
#   d <- par[1]
#   tau <- par[2]
#   
#   q <- mapply(function(d, tau, yi, vi, ycv)
#   { # Compute conditional probabilities for significant and nonsignificant 
#     # effect sizes
#     ifelse(yi > ycv,
#            dnorm(yi, mean = d, sd = sqrt(vi+tau^2), log = TRUE)-
#              pnorm(ycv, mean = d, sd = sqrt(vi+tau^2), lower.tail = FALSE,
#                    log.p = TRUE), 
#            dnorm(yi, mean = d, sd = sqrt(vi+tau^2), log = TRUE)-
#              pnorm(ycv, mean = d, sd = sqrt(vi+tau^2), log.p = TRUE))
#   }, yi = yi, vi = vi, ycv = ycv, MoreArgs = list(d = d, tau = tau))
#   
#   ll <- sum(q) # Compute log-likelihood
#   
#   if (verbose == TRUE)
#   {
#     cat("d = ", d, "; tau^2 = ", tau^2, "; log-lik = ", ll, fill = TRUE, sep = "")
#   }
#   
#   return(ll)
# }

### Function for computing log-likelihood of p-uniform* with moderators
ml_star <- function(par, es, mods, n_bs, par_fixed, transf, verbose)
{
  yi <- es$yi
  vi <- es$vi
  ycv <- es$zcv*sqrt(vi)
  
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
    ifelse(yi > ycv,
           dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) -
             pnorm(ycv, mean = M, sd = sqrt(vi+tau2), lower.tail = FALSE, log.p = TRUE),
           dnorm(yi, mean = M, sd = sqrt(vi+tau2), log = TRUE) - 
             pnorm(ycv, mean = M, sd = sqrt(vi+tau2), log.p = TRUE))
  }, M = M, yi = yi, vi = vi, ycv = ycv, MoreArgs = list(tau2 = tau2))
  
  return(-sum(q))
}
