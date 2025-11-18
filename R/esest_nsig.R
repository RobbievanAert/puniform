### Function for applying p-uniform*
esest_nsig <- function(es, mods, n_bs, par_fixed = rep(NA, n_bs+1), method, boot, con) 
{
  
  yi <- es$yi
  vi <- es$vi
  ycv <- es$zcv*sqrt(vi)
  
  est.ci <- con$est.ci
  tau2.ci <- con$tau2.ci
  par <- con$par
  optimizer <- con$optimizer
  type <- con$type
  int <- con$int
  bounds.int <- con$bounds.int
  tau.int <- con$tau.int
  tol <- con$tol
  maxit <- con$maxit
  verbose <- con$verbose
  reps <- con$reps
  
  if (method == "ML")
  {
    
    if(con$proc.ml == "prof")
    { # Old estimation procedure with ML where the profile log-likelihood functions
      # are iteratively optimized
      
      if (mods != ~1)
      { # Return an error if not the default optimization procedure is requested
        stop("Moderators cannot be included when the profile log-likelihood functions 
        are iteratively optimized. Please use the default estimation procedure.")
      }
      
      ### Starting values for optimization
      tau.est <- sqrt(par[2])
      int <- con$int
      tau.int <- con$tau.int
      tol <- con$tol
      maxit <- con$maxit
      verbose <- con$verbose
      
      stay <- TRUE # In order to stay in while loop
      est <- -999 # Unreasonable estimate to force procedure to use at least two iterations
      i <- 0 # Counter for number of iterations
      
      ### While loop for optimizing profile likelihood functions
      while(stay) {
        
        i <- i+1 # Counter for number of iterations
        
        ### For next iteration new estimate becomes old estimate
        old <- est
        tau.old <- tau.est
        
        ### Optimize profile likelihood function of delta 
        # (suppressWarnings() in order to be able to specify wide search intervals)
        est <- suppressWarnings(optimize(ml_est, int, tau.est, yi, vi, ycv, 
                                         maximum = TRUE)$maximum)
        
        ### Optimize profile likelihood function of tau
        # (suppressWarnings() in order to be able to specify wide search intervals)
        tau.est <- suppressWarnings(optimize(ml_tau, tau.int, est, yi, vi, ycv, 
                                             maximum = TRUE)$maximum)
        
        ### Print intermediate steps if requested
        if (verbose == TRUE)
        {
          cat("est = ", est, "tau.est = ", tau.est, fill = TRUE)
        }
        
        ### Stay in while loop till difference between previous and new estimates
        # is less than tol or maxit equals maximum number of iterations
        stay <- ifelse(abs(tau.est-tau.old) < tol & abs(est-old) < tol | 
                         i == maxit, FALSE, TRUE)
        
      }
      
      if (i == maxit | any(round(est, 3) %in% int | round(tau.est, 3) %in% 
                           tau.int) & round(tau.est, 3) != 0) 
      { # If maximum number of iterations is reached or if estimates are equal to 
        # bounds of interval that was used for optimization return NA except if 
        # estimate of tau is equal to zero
        est <- NA
        tau.est <- NA
        
      }
      
      ### Return NA, because optimization information is only returned if both
      # parameters are estimated at the same time
      out <- NA
      
    } else
    {
      # ### Starting values
      # par <- c(con$stval.d, con$stval.tau)
      # 
      # ### Control arguments of optim(). -1 times con$fnscale to maximize the 
      # # log-likelihood function
      # control.optim <- list(fnscale = -1*con$fnscale, maxit = con$maxit)
      # 
      # ### Optimize log likelihood function
      # out <- optim(par = par, fn = ml_star, yi = yi, vi = vi, ycv = ycv, 
      #              lower = c(-Inf, 0), method = "L-BFGS-B", 
      #              verbose = con$verbose, control = control.optim)
      # 
      # ### Store estimates
      # est <- out$par[1]
      # tau.est <- out$par[2]
      # 
      # ### Return warning message if there are indications for non-convergence
      # if (out$convergence != 0)
      # {
      #   warning("Convergence code is non-zero indicating convergence issues. Try changing the control parameters 'fnscale' and 'maxit' to reach convergence.")
      # }
      
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
        out <- optim(par = par, fn = ml_star, method = optimizer, es = es, 
                     mods = mods, n_bs = n_bs, par_fixed = par_fixed, transf = transf, 
                     verbose = verbose)
      } else if (optimizer == "L-BFGS-B")
      {
        out <- optim(par = par, fn = ml_star, method = optimizer, lower = lower, 
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
      
      ##########################################################################
      
      ##### Compute standard errors #####
      
      ### Estimate the standard errors based on the inverse of the Hessian. Note that
      # we are minimizing the negative log-likelihood function, so the computed Hessian
      # is actually the negative Hessian.
      H <- numDeriv::hessian(func = ml_star, x = c(est, tau2), es = es, mods = mods, 
                             n_bs = n_bs, par_fixed = par_fixed, transf = FALSE,
                             verbose = FALSE)
      inv_H <- try(solve(H), silent = TRUE)
      
      if (inherits(inv_H, what = "try-error"))
      {
        se <- rep(NA, n_bs+1)
        
        warning("Error when inverting Hessian", call. = FALSE)
      } else 
      {
        ### Suppress warning in case of taking the square root of a negative value
        se <- suppressWarnings(sqrt(diag(inv_H)))
      }
      
      ############################################################################
      
      if (any(is.na(c(est, tau2))))
      {
        ci.lb <- NA
        ci.ub <- NA
        tau2.lb <- NA
        tau2.ub <- NA
        L.0 <- NA
        pval.0 <- NA
        L.het <- NA
        pval.het <- NA
      } else 
      {
        
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
              ll0[b] <- -1*optimize(ml_star, interval = c(-10,10), es = es, mods = mods, 
                                    n_bs = n_bs, par_fixed = par_fixed, transf = TRUE, 
                                    verbose = FALSE)$objective
            } else
            { 
              ### Remove the fixed parameters from par
              par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
              
              ll0[b]<- -1*optim(par = par_transf, fn = ml_star, method = "Nelder-Mead", 
                                es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                                transf = TRUE, verbose = FALSE)$value
            }
          }
          
          ### Conduct likelihood-ratio test
          L.0 <- -2*(ll0-ll)
          pval.0 <- pchisq(L.0, df = 1, lower.tail = FALSE)
          
        } else if (type == "Wald" | type == "Wald/profile")
        { # Wald test
          L.0 <- est/se[1:n_bs]
          pval.0 <- 2*pnorm(abs(L.0), lower.tail = FALSE)
        }
        
        ##############################################################################
        
        ##### Test whether there is no (residual) between-study variance #####
        
        if (type == "profile" | type == "Wald/profile")
        { # Likelihood-ratio test
          
          par_fixed <- c(rep(NA, n_bs), 0)
          
          ### Remove the fixed parameters from par
          par_transf <- c(est, log(tau2))[is.na(par_fixed) == TRUE]
          
          if (length(par_transf) == 1)
          { # If only one parameter is estimated, use optimize() instead of optim()
            ### Multiplied by minus 1, because log-likelihood is minimized
            ll0 <- -1*optimize(ml_star, interval = c(-10,10), es = es, mods = mods, 
                               n_bs = n_bs, par_fixed = par_fixed, transf = FALSE, 
                               verbose = FALSE)$objective
          } else
          {
            ll0 <- -1*optim(par = par_transf, fn = ml_star, method = "Nelder-Mead", 
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
        
        ########################################################################
        
        ##### Compute 95% confidence intervals for fixed effects #####
        
        if (type == "profile")
        { # Compute profile likelihood confidence intervals for fixed effects
          
          par_fixed <- rep(NA, n_bs+1)
          
          message("Profile likelihood confidence intervals are computed")
          
          ### Compute lower bound of confidence interval for fixed effects
          ci.lb <- sapply(1:n_bs, FUN = function(ind)
          {
            tmp <- try(uniroot(f = get_profile_ci, 
                               interval = c(est[ind]-est.ci[1], est[ind]),
                               es = es, n_bs = n_bs, par_fixed = par_fixed, mods = mods, 
                               est = est, tau2 = tau2, ind = ind, 
                               chi_cv = qchisq(.95, df = 1), ll = ll,
                               model_type = "puni_star")$root, 
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
                               chi_cv = qchisq(.95, df = 1), ll = ll,
                               model_type = "puni_star")$root, 
                       silent = TRUE)
            
            if (inherits(tmp, what = "try-error"))
            {
              tmp <- NA
            }
            
            return(tmp)
          })
          
        } else if (type == "Wald" | type == "Wald/profile")
        { # Compute Wald confidence intervals for fixed effects
          
          if (all(is.na(se)) == FALSE)
          { # Only compute Wald confidence intervals if se could be computed
            ci.lb <- est - qnorm(.975)*se[1:n_bs]
            ci.ub <- est + qnorm(.975)*se[1:n_bs]
          } else
          {
            ci.lb <- ci.ub <- rep(NA, n_bs)
          }
        }
        
        if (type == "profile" | type == "Wald/profile")
        { # Compute profile likelihood confidence intervals for tau^2
          
          ### Check if lower bound of CI of tau2 is negative
          ll_at_zero <- get_profile_ci(x = log(0), es = es, n_bs = n_bs, 
                                       par_fixed = par_fixed, mods = mods, est = est, 
                                       tau2 = tau2, ind = n_bs+1, 
                                       chi_cv = qchisq(.95, df = 1), ll = ll,
                                       model_type = "puni_star")
          
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
                                   ll = ll, model_type = "puni_star")$root, silent = TRUE)
            
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
                                 ll = ll, model_type = "puni_star")$root, silent = TRUE)
          
          if (!inherits(tau2.ub, what = "try-error"))
          { # If upper bound could be computed transform to tau2 scale
            tau2.ub <- exp(tau2.ub)
          }
          
          if (inherits(tau2.ub, what = "try-error"))
          {
            tau2.ub <- NA
          }
        } else if (type == "Wald")
        { # Compute Wald confidence interval for tau^2
          
          if (all(is.na(se) == FALSE))
          { # Only compute Wald confidence intervals if se could be computed
            tau2.lb <- tau2 - qnorm(.975)*se[length(se)]
            tau2.ub <- tau2 + qnorm(.975)*se[length(se)]
            
            tau2.lb <- ifelse(tau2.lb < 0, 0, tau2.lb)
            tau2.ub <- ifelse(tau2.ub < 0, 0, tau2.ub)
          } else
          {
            tau2.lb <- tau2.ub <- NA
          }
        }
        
        ########################################################################
        
        ### This is the old implementation of getting profile likelihood confidence
        # intervals
        
        # ### Function to compute profile likelihood confidence intervals for the 
        # # average effect size
        # get_profile_ci_est <- function(d, tau, yi, vi, chi_cv, ll, con)
        # {
        #   
        #   ll0 <- optimize(f = ml_star_tau, interval = con$tau.int, d = d, yi = yi, vi = vi,
        #                   ycv = ycv, maximum = TRUE)$objective
        #   
        #   return(-2*(ll0-ll)-chi_cv)
        # }
        # 
        # if(con$proc.ml == "prof")
        # { # Get log-likelihood if optimization via the profile likelihoods was done
        #   
        #   ### Starting values
        #   par <- c(con$stval.d, con$stval.tau)
        #   
        #   ### Control arguments of optim(). -1 times con$fnscale to maximize the 
        #   # log-likelihood function
        #   control.optim <- list(fnscale = -1*con$fnscale, maxit = con$maxit)
        #   
        #   ### Optimize log likelihood function
        #   ll <- optim(par = par, fn = ml_star, yi = yi, vi = vi, ycv = ycv, 
        #               lower = c(-Inf, 0), method = "L-BFGS-B", 
        #               verbose = con$verbose, control = control.optim)$value
        # } else
        # {
        #   ll <- out$value
        # }
        # 
        # tmp.lb <- try(uniroot(f = get_profile_ci_est, interval = c(est-est.ci[1],est), 
        #                       tau = tau.est, yi = yi, vi = vi, 
        #                       chi_cv = qchisq(.95, df = 1), ll = ll, con = con)$root, 
        #               silent = TRUE)
        # 
        # ### Return NA if lower bound could not be estimated
        # lb <- ifelse(inherits(tmp.lb, what = "try-error"), NA, tmp.lb)
        # 
        # tmp.ub <- try(uniroot(f = get_profile_ci_est, interval = c(est,est+est.ci[2]), 
        #                       tau = tau.est, yi = yi, vi = vi, 
        #                       chi_cv = qchisq(.95, df = 1), ll = ll, con = con)$root, 
        #               silent = TRUE)
        # 
        # ### Return NA if lower bound could not be estimated
        # ub <- ifelse(inherits(tmp.ub, what = "try-error"), NA, tmp.ub)
        # 
        # ##########################################################################
        # 
        # ### Estimation of CI tau
        # 
        # ### Function to compute profile likelihood confidence intervals for the 
        # # average effect size
        # get_profile_ci_tau <- function(tau, d, yi, vi, chi_cv, ll, con)
        # {
        #   
        #   ll0 <- optimize(f = ml_star_est, interval = con$int, tau = tau, yi = yi, 
        #                   vi = vi, ycv = ycv, maximum = TRUE)$objective
        #   
        #   return(-2*(ll0-ll)-chi_cv)
        # }
        # 
        # if (get_profile_ci_tau(tau = 0, d = est, yi = yi, vi = vi, 
        #                        chi_cv = qchisq(.95, df = 1), ll = ll, con = con) < 0)
        # { # Set lower bound to zero if it is smaller than 0
        #   tau.lb <- 0
        # } else
        # {
        #   tmp.lb <- try(uniroot(f = get_profile_ci_tau, 
        #                         interval = c(max(c(0, tau.est-con$tau.ci[1])), tau.est), 
        #                         d = est, yi = yi, vi = vi, chi_cv = qchisq(.95, df = 1), 
        #                         ll = ll, con = con)$root, silent = TRUE)
        #   
        #   ### Return NA if lower bound could not be estimated
        #   tau.lb <- ifelse(inherits(tmp.lb, what = "try-error"), NA, tmp.lb)
        # }
        # 
        # tmp.ub <- try(uniroot(f = get_profile_ci_tau, 
        #                       interval = c(tau.est, tau.est+con$tau.ci[2]), 
        #                       d = est, yi = yi, vi = vi, chi_cv = qchisq(.95, df = 1), 
        #                       ll = ll, con = con)$root, silent = TRUE)
        # 
        # ### Return NA if lower bound could not be estimated
        # tau.ub <- ifelse(inherits(tmp.ub, what = "try-error"), NA, tmp.ub)
        
      }
    }
  } else if (method == "P" | method == "LNP")
  {
    
    if (mods != ~1)
    { # Return an error if ML estimation is not used 
      stop("Moderators cannot be included with estimation methods 'P' and 'LNP'.
        Please use the default method 'ML'.")
    }
    
    ### Starting values for root-finding
    tau.ci <- sqrt(tau2.ci)
    tau.est <- 0 # Use tau=0 for first step
    est <- 0 # Use est=0 for first step
    stay <- TRUE # In order to stay in while loop
    i <- 0 # Counter for number of iterations
    est.max <- taus.max <- numeric(10) # Empty objects to examine whether estimates oscilate
    
    ### While loop for estimating effect size and tau
    while(stay) 
    {
      
      i <- i+1 # Counter for number of iterations
      
      if (i > maxit-99)
      { # Check if root finding is not oscillating between two values
        ### Compute difference between new and old estimates
        tau.dif.new <- round(abs(tau.est-tau.old), 2)
        
        if (tau.dif.new == tau.dif)
        { ### If root finding is oscillating, add one half of the difference between old
          # and new values to the estimates
          tau.est <- tau.est+(abs(tau.old-tau.est))/2
        }
      }
      
      if (i > maxit-100)
      { # Compute difference between new and old estimates
        tau.dif <- round(abs(tau.old-tau.est), 2)
      }
      
      ### Compute the bounds for estimating effect size with P and LNP in such a way 
      # that the observed effect sizes are smaller (or larger) than the lower bound and 
      # upper bound for computing the conditional probabilities for estimating 
      # effect size and estimating tau
      int <- bounds_nsig(yi = yi, vi = vi, tau.est = tau.est, ycv = ycv, 
                         method = method, bounds.int = c(bounds.int[1], sort(yi), 
                                                         bounds.int[2]))
      
      ### For next iteration new estimate becomes old estimate
      old <- est
      tau.old <- tau.est
      
      ### Estimate effect size
      est <- try(uniroot(pdist_nsig, interval = c(int[1], int[2]), tau = tau.est, yi = yi,
                         vi = vi, param = "est", ycv = ycv, method = method, val = "es", 
                         cv_P = 0)$root, silent = TRUE)
      
      if (inherits(est, what = "try-error")) 
      { # If effect size could not be estimated return NAs and break out while loop
        est <- NA
        tau.est <- NA
        break
      }
      
      ### Check if tau is larger than zero and otherwise break out while loop
      if (i == 1) 
      {
        tau0 <- round(pdist_nsig(est = est, tau = 0, yi = yi, vi = vi, param = "tau", 
                                 ycv = ycv, method = method, val = "es", cv_P = 0), 3)
        if (method == "P" & tau0 <= 0 | method == "LNP" & tau0 >= 0) 
        {
          tau.est <- 0
          est <- try(uniroot(pdist_nsig, interval = c(int[1], int[2]), tau = tau.est, yi = yi,
                             vi = vi, param = "est", ycv = ycv, method = method, val = "es", 
                             cv_P = 0)$root, silent = TRUE)
          
          if (inherits(est, what = "try-error")) 
          { # If effect size could not be estimated return NAs
            est <- NA
            tau.est <- NA
          }
          
          break
        }
      }
      
      ### Estimate tau
      tau.est <- try(uniroot(pdist_nsig, interval = c(tau.int[1], tau.int[2]), est = est, yi = yi,
                             vi = vi, param = "tau", ycv = ycv, method = method, val = "es", 
                             cv_P = 0)$root, silent = TRUE)
      
      if (inherits(tau.est, what = "try-error")) 
      { # If effect size could not be estimated return NAs and break out while loop
        est <- NA
        tau.est <- NA
        break
      }
      
      ### Print intermediate steps if requested
      if (con$verbose == TRUE) 
      { 
        cat("est = ", est, "tau.est = ", tau.est, fill = TRUE) 
      }
      
      if (i > maxit-10)
      { # Store last ten estimates of effect size and tau before maximum number 
        # of iterations is reached
        est.max[i-maxit+10] <- est
        taus.max[i-maxit+10] <- tau.est
      }
      
      if (i == maxit)
      { # If maximum number of iterations is reached
        if(mean(abs(diff(est.max))) < 0.1)
        { # If algorithm is oscilating, compute mean of last ten estimates of tau and 
          # use this value for tau for estimating est
          tau.est <- mean(taus.max)
          
          ### Compute the bounds for estimating effect size with P and LNP
          int <- bounds_nsig(yi = yi, vi = vi, ycv = ycv, method = method, tau.est = tau.est, 
                             bounds.int = c(bounds.int[1], sort(yi), bounds.int[2]))
          
          ### Estimate effect size
          est <- try(uniroot(pdist_nsig, interval = c(int[1], int[2]), tau = tau.est, yi = yi,
                             vi = vi, param = "est", method = method, val = "es")$root, silent = TRUE)
          
          if (inherits(est, what = "try-error")) 
          { # If effect size could not be estimated return NAs and break out while loop
            est <- NA
            tau.est <- NA
            break
          }
          
        } else
        { # If maximum number of iterations is reached return NA
          est <- NA
          tau.est <- NA
        }
      }
      
      ### Stay in while loop till difference between previous and new estimates
      # is less than tol or maxit equals maximum number of iterations
      stay <- ifelse(abs(tau.est-tau.old) < tol & abs(est-old) < tol | i == maxit, 
                     FALSE, TRUE)
    }
    
    ### If estimates are equal to bounds of interval that was used for optimization
    # return NA except if estimate of tau is equal to zero
    if (any(round(est, 3) %in% int | round(tau.est, 3) %in% tau.int) & round(tau.est, 3) != 0) 
    {
      tau.est <- NA
      est <- NA
    }
    
    if (is.na(est) == TRUE & is.na(tau.est) == TRUE)
    { # If effect size and tau could not be estimated return NAs for CIs
      ci.lb <- NA
      ci.ub <- NA
      tau.lb <- NA
      tau.ub <- NA
    } else 
    {
      ### Estimate CI of est ###
      ci.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(est-con$est.ci[1], est),  
                                         tau = tau.est, yi = yi, vi = vi, param = "est", 
                                         ycv = ycv, method = method, val = "ci.lb", 
                                         get_cv_P(length(yi)))$root, silent = TRUE))
      
      if (inherits(ci.lb, what = "try-error")) 
      { # Check if lower bound could be estimated
        ci.lb <- NA
      } 
      
      ci.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(est, est+con$est.ci[2]),  
                                         tau = tau.est, yi = yi, vi = vi, param = "est", 
                                         ycv = ycv, method = method, val = "ci.ub", 
                                         get_cv_P(length(yi)))$root, silent = TRUE))
      
      if (inherits(ci.ub, what = "try-error")) 
      { # Check if upper bound could be estimated
        ci.ub <- NA
      } 
      
      ### Estimate CI of tau ###
      if (method == "P")
      {
        if (pdist_nsig(est = est, tau = 0, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                       method = method, val = "ci.ub", cv_P = get_cv_P(length(yi))) < 0)
        { # Return 0 (null set) if lower and upper bound are negative
          tau.lb <- tau.ub <- 0 
        } else if (pdist_nsig(est = est, tau = 0, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                              method = method, val = "ci.lb", cv_P = get_cv_P(length(yi))) < 0) 
        { # Truncate lower bound to zero if it is negative
          tau.lb <- 0 
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+tau.ci[1]), 
                                                 est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                 method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                                         silent = TRUE))
          
          if (inherits(tau.ub, what = "try-error")) 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else 
        { # Estimate lower and upper bound
          
          tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(max(0, tau.est-tau.ci[2]), tau.est), 
                                                 est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                 method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                                         silent = TRUE))  
          
          if (inherits(tau.lb, what = "try-error")) 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          } 
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+tau.ci[1]), 
                                                 est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                 method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                                         silent = TRUE))
          
          if (inherits(tau.ub, what = "try-error")) 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          }
        }
      } else if (method == "LNP")
      {
        if (pdist_nsig(est = est, tau = 0, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                       method = method, val = "ci.ub", cv_P = get_cv_P(length(yi))) > 0)
        { # Return 0 (null set) if lower and upper bound are negative
          tau.lb <- tau.ub <- 0
        } else if (pdist_nsig(est = est, tau = 0, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                              method = method, val = "ci.lb", cv_P = get_cv_P(length(yi))) > 0) 
        {
          tau.lb <- 0 # Truncate lower bound to zero if it is negative
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+tau.ci[1]), 
                                                 est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                 method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                                         silent = TRUE))
          
          if (inherits(tau.ub, what = "try-error")) 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else 
        { # Estimate lower and upper bound
          if (tau.ci[2] == 0)
          { # If user did not specify a value for search interval, search from 0 to tau.est
            tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est), est = est, 
                                                   yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                   method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                                           silent = TRUE))  
          } else 
          { # Estimate lower and upper bound
            tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(max(0, tau.est-tau.ci[2]), tau.est), 
                                                   est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                   method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                                           silent = TRUE))  
          }
          
          if (inherits(tau.lb, what = "try-error")) 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          } 
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+tau.ci[1]), 
                                                 est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                                 method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                                         silent = TRUE))
          
          if (inherits(tau.ub, what = "try-error")) 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          }
        }
      }
    }
    
    ############################################################################
    
    ##### Test of no between-study variance #####
    
    ### Estimate effect size with tau=0
    est0 <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(-4, 4), tau = 0, 
                                         yi = yi, vi = vi, param = "est", ycv = ycv, 
                                         method = method, val = "es", cv_P = 0)$root, 
                                 silent = TRUE))
    
    if (inherits(est0, what = "try-error")) 
    {
      est0 <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(-10, 10), tau = 0, 
                                           yi = yi, vi = vi, param = "est", ycv = ycv, 
                                           method = method, val = "es", cv_P = 0)$root, 
                                   silent = TRUE))
    }
    
    if (inherits(est0, what = "try-error"))
    { # If effect size cannot be estimated, return NA
      L.het <- NA
      pval.het <- NA
      pval.boot <- NA
    } else 
    {
      ### Compute conditional probabilities at est0
      tr.q <- trq(est = est0, tau = 0, yi = yi, vi = vi, ycv = ycv, param = "est")
      
      het.q <- 2*abs(tr.q-0.5) # Compute heterogeneity statistic
      L.het <- sum(-log(1-het.q))
      pval.het <- pgamma(L.het, length(yi), 1, lower.tail = FALSE)
      
      if (boot == TRUE)
      { # If boot == TRUE, bootstrapped p-value is computed
        
        ### Conduct bootstrapping
        L.het.boot <- replicate(reps, expr = boot_het(k = length(yi), est0 = est0, vi = vi, 
                                                      ycv = ycv, method = method,
                                                      con = con))
        
        ### Compute p-value with bootstrapping
        pval.het <- length(L.het.boot[L.het.boot > L.het & !is.na(L.het.boot)])/reps
        
      } else 
      {
        pval.het <- NA
      }
    }
    
    ############################################################################
    
    ### Return NA, because optimization information is only returned if both
    # parameters are estimated at the same time with method = "ML"
    out <- NA
    
    ### Estimates of tau^2 are returned and not of tau in this function
    tau2 <- tau.est^2
    tau2.lb <- tau.lb^2
    tau2.ub <- tau.ub^2
    
    ### Standard errors, L.0, and pval.0 are NA for "P" and "LNP"
    se <- NA
    L.0 <- NA
    pval.0 <- NA
    
  }
  
  return(list(est = est, tau2 = tau2, se = se, L.0 = L.0, pval.0 = pval.0,
              L.het = L.het, pval.het = pval.het, ci.lb = ci.lb, ci.ub = ci.ub, 
              tau2.lb = tau2.lb, tau2.ub = tau2.ub, optim.info = out))
}