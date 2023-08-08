### Function for applying hybrid method
hy <- function(es, measure, side, mods, n_bs, par_fixed = rep(NA, n_bs+1), con) 
{
  
  int <- con$int
  est.ci <- con$est.ci
  tau2.ci <- con$tau2.ci
  verbose <- con$verbose
  tol <- con$tol 
  par <- con$par
  implementation <- con$implementation
  type <- con$type
  optimizer <- con$optimizer
  
  if (implementation == "two")
  { # If the implementation of van Aert and van Assen (2018) is used with only two
    # studies (one original and one replication study)
    
    ### Apply bisection method for effect size
    est <- try(bisect(func = pdist_hy_helper, lo = int[1], hi = int[2], es = es, 
                         val = "est", tol = tol, verbose = verbose), silent = TRUE)
    
    if (inherits(est, what = "try-error"))
    { # If estimate cannot be computed, return NA
      est <- NA
    }
    
    ### Apply bisection method for lower bound
    ci.lb <- try(bisect(func = pdist_hy_helper, lo = est-est.ci[1], 
                           hi = est, es = es, 
                           val = "ci.lb", cv.P = get_cv_P(nrow(es)), 
                           tol = tol, verbose = verbose), silent = TRUE)
    
    if (inherits(ci.lb, what = "try-error")) 
    {
      ci.lb <- NA
    }
    
    ### Apply bisection method for upper bound
    ci.ub <- try(bisect(func = pdist_hy_helper, lo = est, hi = est+est.ci[2], 
                           es = es, val = "ci.ub", cv.P = nrow(es) - get_cv_P(nrow(es)), 
                           tol = tol, verbose = verbose), silent = TRUE)
    
    if (inherits(ci.ub, what = "try-error")) 
    {
      ci.ub <- NA
    }
    
    ### Test of H0: mu = 0
    q <- pdist_hy(d = 0, es = es, val = "est")$q
    
    if (length(q) == 2)
    { # If there is only one original study and one replication
      L.0 <- ifelse(sum(q) < 1, sum(q), 2-sum(q))  # Compute probability density  
      pval.0 <- ifelse(sum(q) < 1, 0.5*sum(q)^2, -0.5*sum(q)^2 + 2*sum(q)-1)
      pval.0 <- ifelse(pval.0 > 0.5, (1-pval.0)*2, pval.0*2) # Compute two-tailed p-value  
    } else
    { # If there are more than two studies
      L.0 <- (sum(q)-nrow(es)*0.5)/sqrt(nrow(es)/12)
      pval.0 <- pnorm(L.0)
      pval.0 <- ifelse(pval.0 > 0.5, (1-pval.0)*2, pval.0*2) # Compute two-tailed p-value
    }
    
    ### Objects that are not obtained in case of only one original and replication study
    tau2 <- se <- L.het <- pval.het <- tau2.lb <- tau2.ub <- 
      out <- NA
    
    ############################################################################
    ############################################################################
    ############################################################################
    
  } else if (implementation == "multiple")
  { # If the implementation is used that allows for more than two studies
    
    ### Set lower bounds for optimization if "L-BFGS-B" is the optimizer
    if (optimizer == "L-BFGS-B") lower <- c(rep(-Inf, n_bs), 0)
    
    if (any(is.na(par_fixed) == FALSE))
    { # If there are parameters fixed for hypothesis testing
      
      ### Remove the fixed parameters from par
      par <- par[is.na(par_fixed) == TRUE]
      
      ### Remove lower bound in case of a fixed parameter and "L-BFGS-B" as optimizer
      if (optimizer == "L-BFGS-B") lower <- lower[is.na(par_fixed) == TRUE]
    }
    
    ### If unconstrained optimization is used, optimize exp(tau2) rather than tau2
    # on the interval 0 to Inf (default settings)
    transf <- ifelse(optimizer != "L-BFGS-B", TRUE, FALSE)
    
    ############################################################################
    
    ##### Estimate parameters #####
    
    if (optimizer != "L-BFGS-B")
    {
      out <- optim(par = par, fn = ml_hy, method = optimizer, es = es, 
                   mods = mods, n_bs = n_bs, par_fixed = par_fixed, transf = transf, 
                   verbose = verbose)
    } else if (optimizer == "L-BFGS-B")
    {
      out <- optim(par = par, fn = ml_hy, method = optimizer, lower = lower, 
                   es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                   transf = transf, verbose = verbose)
    }
    
    if (out$convergence != 0)
    { # Return warning message if optim returns a nonzero convergence code
      warning("Convergence code is nonzero suggesting that nonconvergence has occured")
    }
    
    ### Store the estimated parameters
    est <- out$par[1:n_bs]
    tau2 <- out$par[n_bs+1]
    
    ### Take the exponent of tau2 if unconstrained optimization was used
    tau2 <- ifelse(transf == TRUE, exp(tau2), tau2)
    
    ### Store log-likelihood
    ll <- -1*out$value
    
    ############################################################################
    
    ##### Compute standard errors #####
    
    ### Estimate the standard errors based on the inverse of the Hessian. Note that
    # we are minimizing the negative log-likelihood function, so the computed Hessian
    # is actually the negative Hessian.
    H <- numDeriv::hessian(func = ml_hy, x = c(est, tau2), es = es, mods = mods, 
                           n_bs = n_bs, par_fixed = par_fixed, transf = FALSE,
                           verbose = FALSE)
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
      se <- sqrt(diag(inv_H))  
    }
    
    ############################################################################
    
    ##### Test whether the fixed effects are different from zero #####
    
    if (type == "profile")
    { # Likelihood-ratio test
      ll0 <- numeric(n_bs)
      
      for (b in 1:n_bs)
      {
        
        par_fixed <- rep(NA, n_bs+1)
        par_fixed[b] <- 0
        
        if (n_bs == 1)
        { # If only one parameter is estimated, use optimize() instead of optim()
          ### Multiplied by minus 1, because log-likelihood is minimized
          ll0[b] <- -1*optimize(ml_hy, interval = c(-10,10), es = es, mods = mods, 
                                n_bs = n_bs, par_fixed = par_fixed, transf = TRUE, 
                                verbose = FALSE)$objective
        } else
        { 
          ### Remove the fixed parameters from par
          par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
          
          ll0[b]<- -1*optim(par = par_transf, fn = ml_hy, method = "Nelder-Mead", 
                            es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                            transf = TRUE, verbose = FALSE)$value
        }
      }
      
      ### Conduct likelihood-ratio test
      L.0 <- -2*(ll0-ll)
      pval.0 <- pchisq(L.0, df = 1, lower.tail = FALSE)
      
    } else if (type == "Wald")
    { # Wald test
      L.0 <- est/se[1:n_bs]
      pval.0 <- 2*pnorm(abs(L.0), lower.tail = FALSE)
    }
    
    ##############################################################################
    
    ##### Test whether there is (residual) between-study variance #####
    
    if (type == "profile")
    { # Likelihood-ratio test
      
      par_fixed <- c(rep(NA, n_bs), 0)
      
      ### Remove the fixed parameters from par
      par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
      
      if (length(par_transf) == 1)
      { # If only one parameter is estimated, use optimize() instead of optim()
        ### Multiplied by minus 1, because log-likelihood is minimized
        ll0 <- -1*optimize(ml_hy, interval = c(-10,10), es = es, mods = mods, 
                           n_bs = n_bs, par_fixed = par_fixed, transf = FALSE, 
                           verbose = FALSE)$objective
      } else
      {
        ll0 <- -1*optim(par = par_transf, fn = ml_hy, method = "Nelder-Mead", 
                        es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                        transf = FALSE, verbose = FALSE)$value
      }
      
      ### Conduct likelihood-ratio test
      L.het <- -2*(ll0-ll)
      
      ### 0.5 x chisq, because the tested null-hypothesis H0: tau2 = 0 is on the 
      # boundary of the parameter space. See Andrews (2001) and Molenberghs and 
      # Verbeke (2012)
      pval.het <- 0.5*pchisq(L.het, df = 1, lower.tail = FALSE)
      
    } else if (type == "Wald")
    { # Wald test
      L.het <- tau2/se[length(se)]
      pval.het <- 2*pnorm(abs(L.het), lower.tail = FALSE)
    }
    
    ##############################################################################
    
    ##### Compute 95% confidence intervals #####
    
    if (type == "profile")
    { # Compute profile likelihood confidence intervals
      
      par_fixed <- rep(NA, n_bs+1)
      
      ### Function used for estimating profile likelihood confidence intervals
      get_profile_ci <- function(x, es, n_bs, par_fixed, mods, est, tau2, ind, 
                                 chi_cv, ll)
      {
        par_fixed[ind] <- x
        
        ### Remove the fixed parameters from par
        par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
        
        ### Re-estimate model without the fixed parameter. Multiply with -1 to 
        # get the log-likelihood, because ml_mods() returns the negative log-likelihood
        if (length(par_transf) == 1)
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
      
      message("Profile likelihood confidence intervals are computed")
      
      ### Compute lower bound of confidence interval for fixed effects
      ci.lb <- sapply(1:n_bs, FUN = function(ind)
      {
        tmp <- try(uniroot(f = get_profile_ci, 
                           interval = c(est[ind]-est.ci[1], est[ind]),
                           es = es, n_bs = n_bs, par_fixed = par_fixed, mods = mods, 
                           est = est, tau2 = tau2, ind = ind, 
                           chi_cv = qchisq(.95, df = 1), ll = ll)$root, 
                   silent = TRUE)
        
        if (inherits(tmp, what = "try-error"))
        {
          tmp <- NA
        }
        
        return(tmp)
      })
      
      ### Compute upper bound of confidence interval for fixed effects
      ci.ub <- sapply(1:n_bs, FUN = function(ind)
      {
        tmp <- try(uniroot(f = get_profile_ci, 
                           interval = c(est[ind], est.ci[2]+est[ind]),
                           es = es, n_bs = n_bs, par_fixed = par_fixed, mods = mods, 
                           est = est, tau2 = tau2, ind = ind, 
                           chi_cv = qchisq(.95, df = 1), ll = ll)$root, 
                   silent = TRUE)
        
        if (inherits(tmp, what = "try-error"))
        {
          tmp <- NA
        }
        
        return(tmp)
      })
      
      ### Check if lower bound of CI of tau2 is negative
      ll_at_zero <- get_profile_ci(x = log(0), es = es, n_bs = n_bs, 
                                   par_fixed = par_fixed, mods = mods, est = est, 
                                   tau2 = tau2, ind = n_bs+1, 
                                   chi_cv = qchisq(.95, df = 1), ll = ll)
      
      if (ll_at_zero < 0)
      {
        tau2.lb <- 0
      } else
      {
        tau2.lb <- try(uniroot(f = get_profile_ci,
                                     interval = log(c(max(1e-50,tau2-tau2.ci[1]), 
                                                      tau2)),
                                     es = es, n_bs = n_bs, par_fixed = par_fixed, 
                                     mods = mods, est = est, tau2 = tau2,
                                     ind = n_bs+1, chi_cv = qchisq(.95, df = 1), 
                                     ll = ll)$root, silent = TRUE)
        
        if (!inherits(tau2.lb, what = "try-error"))
        { # If lower bound could be computed transform to tau2 scale
          tau2.lb <- exp(tau2.lb)
        }
      }
      
      if (inherits(tau2.lb, what = "try-error"))
      {
        tau2.lb <- NA
      }
      
      tau2.ub <- try(uniroot(f = get_profile_ci,
                                   interval = log(c(tau2, tau2+tau2.ci[2])),
                                   es = es, n_bs = n_bs, par_fixed = par_fixed, 
                                   mods = mods, est = est, tau2 = tau2,
                                   ind = n_bs+1, chi_cv = qchisq(.95, df = 1), 
                                   ll = ll)$root, silent = TRUE)
      
      if (!inherits(tau2.ub, what = "try-error"))
      { # If upper bound could be computed transform to tau2 scale
        tau2.ub <- exp(tau2.ub)
      }
      
      if (inherits(tau2.ub, what = "try-error"))
      {
        tau2.ub <- NA
      }
      
    } else if (type == "Wald")
    { # Wald confidence interval
      
      ### Only compute Wald confidence intervals if se could be computed
      if (all(is.na(se) == FALSE))
      {
        ci.lb <- est - qnorm(.975)*se[1:n_bs]
        ci.ub <- est + qnorm(.975)*se[1:n_bs]
        
        tau2.lb <- tau2 - qnorm(.975)*se[length(se)]
        tau2.ub <- tau2 + qnorm(.975)*se[length(se)]
        
        tau2.lb <- ifelse(tau2.lb < 0, 0, tau2.lb)
        tau2.ub <- ifelse(tau2.ub < 0, 0, tau2.ub)
      } else
      {
        ci.lb <- ci.ub <- rep(NA, n_bs)
        tau2.lb <- tau2.ub <- NA
      }
    }
  }
  
  if (side == "left") 
  { # Re-mirror estimates
    est <- est * -1
    tmp <- ci.ub
    ci.ub <- ci.lb * -1
    ci.lb <- tmp * -1
    
    if (type == "Wald")
    {
      L.0 <- L.0 * -1
    }
  }
  
  return(list(est = est, tau2 = tau2, se = se, 
              L.0 = L.0, pval.0 = pval.0, L.het = L.het, 
              pval.het = pval.het, ci.lb = ci.lb, ci.ub = ci.ub, 
              tau2.lb = tau2.lb, tau2.ub = tau2.ub, 
              optim.info = out))
} 