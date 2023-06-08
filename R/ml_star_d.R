### Function for computing the profile log likelihood of p-uniform* for the
# average effect size
ml_star_d <- function(d, tau, yi, vi, ycv) 
{
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
  
  ll <- sum(q) # Compute log-likelihood
  
  return(ll)
}