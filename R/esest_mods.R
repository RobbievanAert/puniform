### Function for estimating the parameters with p-uniform* in case of moderators
esest_mods <- function(es, mods, n_bs, par_fixed = rep(NA, n_bs+1), type, con)
{
  par <- con$par # Load starting values
  lower <- c(rep(-Inf, n_bs), 0) # Lower bounds if "L-BFGS-B" is the optimizer
  
  ### Set lower bounds for optimization if "L-BFGS-B" is the optimizer
  if (con$optimizer == "L-BFGS-B") lower <- c(rep(-Inf, n_bs), 0)
  
  if (any(is.na(par_fixed) == FALSE))
  { # If there are parameters fixed for hypothesis testing
    
    ### Remove the fixed parameters from par
    par <- par[is.na(par_fixed) == TRUE]
    
    ### Remove lower bound in case of a fixed parameter and "L-BFGS-B" as optimizer
    if (con$optimizer == "L-BFGS-B") lower <- lower[is.na(par_fixed) == TRUE]
  }
  
  ### If unconstrained optimization is used, optimize exp(tau2) rather than tau2
  # on the interval 0 to Inf (default settings)
  transf <- ifelse(con$optimizer != "L-BFGS-B", TRUE, FALSE)
  
  if (con$optimizer != "L-BFGS-B")
  {
    out <- optim(par = par, fn = ml_mods, es = es, n_bs = n_bs, par_fixed = par_fixed, 
                 mods = mods, transf = transf)
  } else if (con$optimizer == "L-BFGS-B")
  {
    out <- optim(par = par, fn = ml_mods, method = "L-BFGS-B", lower = lower,
                 es = es, n_bs = n_bs, par_fixed = par_fixed, mods = mods, 
                 transf = transf)
  }
  
  ### Take into account any fixed parameters
  par_est <- par_fixed
  par_est[which(is.na(par_est))] <- out$par
  
  ### Store the estimated parameters
  est <- par_est[1:n_bs]
  tau2.est <- ifelse(is.na(par_fixed[length(par_fixed)]) == TRUE, par_est[n_bs+1],
                     par_fixed[length(par_fixed)])
  
  ### Take the exponent of tau2 if unconstrained optimization was used
  tau2.est <- ifelse(transf == TRUE, exp(tau2.est), tau2.est)
  
  ### -1*ll was used, because optim() minimizes the log-likelihood. Multiply with
  # -1 to get the log-likelihood
  ll <- -1*out$value
  
  ### If parameters are fixed, do not compute confidence intervals for now
  if (all(is.na(par_fixed) == TRUE))
  {
    ### Estimate the standard errors based on the inverse of the Hessian. Note that
    # we are minimizing the negative log-likelihood function, so the computed Hessian
    # is actually the negative Hessian.
    H <- numDeriv::hessian(func = ml_mods, x = c(est, tau2.est), es = es, n_bs = n_bs,
                           par_fixed = par_fixed, mods = mods, transf = FALSE)
    inv_H <- try(solve(H), silent = TRUE)
    
    if (inherits(inv_H, what = "try-error"))
    {
      se <- rep(NA, n_bs+1)
      
      warning("Error when inverting Hessian", call. = FALSE)
    } else if (any(diag(inv_H) < 0))
    {
      se <- rep(NA, n_bs+1)
      
      warning("Error when inverting Hessian", call. = FALSE)
    } else 
    {
      se <- sqrt(diag(H))  
    }
    
    if (type == "profile")
    { # Compute profile likelihood confidence intervals
      
      message("Profile likelihood confidence intervals are computed")
      
      get_profile_ci <- function(x, out, uni_lower, es, n_bs, par_fixed, 
                                 mods, ind, chi_cv, ll, transf, con)
      {
        par_fixed[ind] <- x
        
        ### Remove the fixed parameters from par
        par <- out$par[is.na(par_fixed) == TRUE]
        
        ### Remove lower bound in case of a fixed parameter and "L-BFGS-B" as optimizer
        if (con$optimizer == "L-BFGS-B") lower <- uni_lower[is.na(par_fixed) == TRUE]
        
        ### Re-estimate model without the fixed parameter. Multiply with -1 to 
        # get the log-likelihood, because ml_mods() returns the negative log-likelihood
        if (con$optimizer != "L-BFGS-B")
        {
          ll0 <- -1*optim(par = par, fn = ml_mods, es = es, n_bs = n_bs, 
                          par_fixed = par_fixed, mods = mods, transf = transf)$value
        } else if (con$optimizer == "L-BFGS-B")
        {
          ll0 <- -1*optim(par = par, fn = ml_mods, method = "L-BFGS-B", lower = lower,
                          es = es, n_bs = n_bs, par_fixed = par_fixed, mods = mods, 
                          transf = transf)$value
        }
        
        return(-2*(ll0-ll)-chi_cv)
      }
      
      ### Compute lower bound of confidence interval for fixed effects
      ci.lb <- sapply(1:n_bs, FUN = function(ind)
      {
        tmp <- try(uniroot(f = get_profile_ci, interval = c(con$int[1], out$par[ind]),
                           out = out, uni_lower = lower, es = es, n_bs = n_bs,
                           par_fixed = par_fixed, mods = mods, ind = ind, 
                           chi_cv = qchisq(.95, df = 1), ll = ll, transf = transf,
                           con = con)$root, silent = TRUE)
        
        if (inherits(tmp, what = "try-error"))
        {
          tmp <- NA
        }
        
        return(tmp)
      })
      
      ### Compute upper bound of confidence interval for fixed effects
      ci.ub <- sapply(1:n_bs, FUN = function(ind)
      {
        tmp <- try(uniroot(f = get_profile_ci, interval = c(out$par[ind], con$int[2]),
                           out = out, uni_lower = lower, es = es, n_bs = n_bs,
                           par_fixed = par_fixed, mods = mods, ind = ind, 
                           chi_cv = qchisq(.95, df = 1), ll = ll, transf = transf,
                           con = con)$root, silent = TRUE)
        
        if (inherits(tmp, what = "try-error"))
        {
          tmp <- NA
        }
        
        return(tmp)
      })
      
      ### Check if lower bound of CI of tau2 is negative
      ll_at_zero <- ifelse(con$optimizer == "L-BFGS-B", 
                           get_profile_ci(x = 0, out = out, uni_lower = lower, 
                                          es = es, n_bs = n_bs, par_fixed = par_fixed, 
                                          mods = mods, ind = n_bs+1, 
                                          chi_cv = qchisq(.95, df = 1), ll = ll, 
                                          transf = transf, con = con),
                           get_profile_ci(x = log(0), out = out, uni_lower = lower, 
                                          es = es, n_bs = n_bs, par_fixed = par_fixed, 
                                          mods = mods, ind = n_bs+1, 
                                          chi_cv = qchisq(.95, df = 1), ll = ll, 
                                          transf = transf, con = con))
      
      if (ll_at_zero < 0)
      {
        ci.lb.tau2.est <- 0
      } else
      {
        if (con$optimizer == "L-BFGS-B")
        {
          ci.lb.tau2.est <- try(uniroot(f = get_profile_ci,
                                        interval = c(max(0,tau2.est-con$tau.ci[1]^2), tau2.est),
                                        out = out, uni_lower = lower, es = es, 
                                        n_bs = n_bs, par_fixed = par_fixed, 
                                        mods = mods, ind = n_bs+1, 
                                        chi_cv = qchisq(.95, df = 1), ll = ll, 
                                        transf = transf, con = con)$root, silent = TRUE)
        } else if (con$optimizer != "L-BFGS-B")
        {
          ci.lb.tau2.est <- try(uniroot(f = get_profile_ci,
                                        interval = log(c(max(1e-50,tau2.est-con$tau.ci[1]^2), tau2.est)),
                                        out = out, uni_lower = lower, es = es, 
                                        n_bs = n_bs, par_fixed = par_fixed, 
                                        mods = mods, ind = n_bs+1, 
                                        chi_cv = qchisq(.95, df = 1), ll = ll, 
                                        transf = transf, con = con)$root, silent = TRUE)
          
          if (!inherits(ci.lb.tau2.est, what = "try-error"))
          { # If lower bound could be computed transform to tau2 scale
            ci.lb.tau2.est <- exp(ci.lb.tau2.est)
          }
        }
        
        if (inherits(ci.lb.tau2.est, what = "try-error"))
        {
          ci.lb.tau2.est <- NA
        }
        
      }
      
      if (con$optimizer == "L-BFGS-B")
      {
        ci.ub.tau2.est <- try(uniroot(f = get_profile_ci,
                                      interval = c(tau2.est, tau2.est+con$tau.ci[2]^2),
                                      out = out, uni_lower = lower, es = es, 
                                      n_bs = n_bs, par_fixed = par_fixed, 
                                      mods = mods, ind = n_bs+1, 
                                      chi_cv = qchisq(.95, df = 1), ll = ll, 
                                      transf = transf, con = con)$root, silent = TRUE)
      } else if (con$optimizer != "L-BFGS-B")
      {
        ci.ub.tau2.est <- try(uniroot(f = get_profile_ci,
                                      interval = log(c(tau2.est, tau2.est+con$tau.ci[2]^2)),
                                      out = out, uni_lower = lower, es = es, 
                                      n_bs = n_bs, par_fixed = par_fixed, 
                                      mods = mods, ind = n_bs+1, 
                                      chi_cv = qchisq(.95, df = 1), ll = ll, 
                                      transf = transf, con = con)$root, silent = TRUE)
        
        if (!inherits(ci.ub.tau2.est, what = "try-error"))
        { # If upper bound could be computed transform to tau2 scale
          ci.ub.tau2.est <- exp(ci.ub.tau2.est)
        }
      }
      
      if (inherits(ci.ub.tau2.est, what = "try-error"))
      {
        ci.ub.tau2.est <- NA
      }
      
    } else if (type == "Wald")
    { # Wald confidence interval
      
      ### Only compute Wald confidence intervals if se could be computed
      if (all(is.na(se) == FALSE))
      {
        ci.lb <- est - qnorm(.975)*se[1:n_bs]
        ci.ub <- est + qnorm(.975)*se[1:n_bs]
        
        ci.lb.tau2.est <- tau2.est - qnorm(.975)*se[length(se)]
        ci.ub.tau2.est <- tau2.est + qnorm(.975)*se[length(se)]
        
        ci.lb.tau2.est <- ifelse(ci.lb.tau2.est < 0, 0, ci.lb.tau2.est)
        ci.ub.tau2.est <- ifelse(ci.ub.tau2.est < 0, 0, ci.ub.tau2.est)
      } else
      {
        ci.lb <- ci.ub <- rep(NA, n_bs)
        ci.lb.tau2.est <- ci.ub.tau2.est <- NA
      }
      
    }
  } else
  {
    ci.lb <- ci.ub <- rep(NA, n_bs)
    ci.lb.tau2.est <- ci.ub.tau2.est <- NA
    se <- rep(NA, n_bs+1)
  }
  
  return(list(est = est, tau2.est = tau2.est, se = se, ll = ll, ci.lb = ci.lb, 
              ci.ub = ci.ub, ci.lb.tau2.est = ci.lb.tau2.est, 
              ci.ub.tau2.est = ci.ub.tau2.est, optim.info = out))
}