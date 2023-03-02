### Function for estimating the parameters with p-uniform* in case of moderators
esest_mods <- function(es, mods, con)
{
  par <- con$par # Load starting values
  
  n_bs <- ncol(model.matrix(mods), data = es) # Number of regression parameters
    
  lower <- c(rep(-Inf, n_bs), 0) # Lower bound for optimizing
  
  out <- optim(par = par, fn = ml_mods, method = "L-BFGS-B", lower = lower, 
        es = es, mods = mods, n_bs = n_bs)
  
  return(list(list(est = out$par[1:n_bs], tau.est = out$par[n_bs+1]),
              optim.info = out))
}