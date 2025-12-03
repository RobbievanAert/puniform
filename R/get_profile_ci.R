### Function used for estimating profile likelihood confidence intervals hybrid 
# method
get_profile_ci <- function(x, es, n_bs, par_fixed, mods, est, tau2, ind, 
                           chi_cv, ll, tau2.fixed)
{
  par_fixed[ind] <- x
  
  ### Remove the fixed parameters from par
  par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
  
  ### Re-estimate model without the fixed parameter. Multiply with -1 to 
  # get the log-likelihood, because ml_hy() returns the negative log-likelihood
  if (n_bs == 1 & is.na(tau2.fixed) == FALSE)
  { # If only one parameter is estimated and tau2 is set to a particular 
    # value, the null model is the log-likelihood where the effect size 
    # is zero
    ll0 <- -1*ml_hy(par = x, es = es, mods = mods, n_bs = n_bs, 
                    par_fixed = c(NA, 0), transf = FALSE, verbose = FALSE)
  } else if (length(par_transf) == 1)
  { # If only one parameter is estimated, use optimize() instead of optim()
    ll0 <- -1*optimize(ml_hy, interval = c(-10,10), es = es, mods = mods, 
                       n_bs = n_bs, par_fixed = par_fixed, transf = TRUE, 
                       verbose = FALSE)$objective
  } else
  {
    ll0 <- -1*optim(par = par_transf, fn = ml_hy, method = "Nelder-Mead", 
                    es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                    transf = TRUE, verbose = FALSE)$value
  }
  
  return(-2*(ll0-ll)-chi_cv)
}