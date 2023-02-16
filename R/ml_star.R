### Function for computing log likelihood of p-uniform*
ml_star <- function(par, yi, vi, ycv) 
{
  
  d <- par[1]
  tau <- par[2]
  
  q <- mapply(function(d, tau, yi, vi, ycv)
  { # Compute conditional probabilities for significant and nonsignificant 
    # effect sizes
    ifelse(yi > ycv,
           dnorm(yi, mean = d, sd = sqrt(vi+tau^2), log = TRUE)-
             pnorm(ycv, mean = d, sd = sqrt(vi+tau^2), lower.tail = FALSE,
                   log.p = TRUE), 
           dnorm(yi, mean = d, sd = sqrt(vi+tau^2), log = TRUE)-
             pnorm(ycv, mean = d, sd = sqrt(vi+tau^2), log.p = TRUE))
  }, yi = yi, vi = vi, ycv = ycv, MoreArgs = list(d = d, tau = tau))
  
  -1*sum(q) # Compute log likelihood
  
}