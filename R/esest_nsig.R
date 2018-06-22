### Function for estimation with p-uniform*
esest_nsig <- function(yi, vi, int, tau.int, ycv, method, con) 
{
  
  if (method == "ML")
  {
    
    ### Starting values for optimization
    tau.est <- con$stval.tau
    int <- con$int
    tau.int <- con$tau.int
    est.ci <- con$est.ci
    tau.ci <- con$tau.ci
    tol <- con$tol
    max.iter <- con$max.iter
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
      # is less than tol or max.iter equals maximum number of iterations
      stay <- ifelse(abs(tau.est-tau.old) < tol & abs(est-old) < tol | 
                       i == max.iter, FALSE, TRUE)
      
    }
    
    if (i == max.iter | any(round(est, 3) %in% int | round(tau.est, 3) %in% 
                                tau.int) & round(tau.est, 3) != 0) 
    { # If maximum number of iterations is reached or if estimates are equal to 
      # bounds of interval that was used for optimization return NA except if 
      # estimate of tau is equal to zero
      est <- NA
      tau.est <- NA
      lb <- NA
      ub <- NA
      tau.lb <- NA
      tau.ub <- NA
    } else 
    {
      ### Estimation CI (suppressWarnings() in order to be able to specify wide search intervals)
      tmp.lb <- suppressWarnings(try(uniroot(get_LR_est, interval = c(est-est.ci[1], est), 
                                             yi = yi, vi = vi, est = est, tau.est = tau.est, 
                                             ycv = ycv)$root, silent = TRUE))
      
      lb <- ifelse(class(tmp.lb) == "try-error", NA, tmp.lb) # Return NA if lower bound could not be estimated
      
      tmp.ub <- suppressWarnings(try(uniroot(get_LR_est, interval = c(est, est+est.ci[2]), 
                                             yi = yi, vi = vi, est = est, tau.est = tau.est, 
                                             ycv = ycv)$root, silent = TRUE))
      
      ub <- ifelse(class(tmp.ub) == "try-error", NA, tmp.ub) # Return NA if upper bound could not be estimated
      
      ### Estimation of CI tau
      if (get_LR_tau(prof.tau = 0, yi = yi, vi = vi, est = est, tau.est = tau.est, ycv = ycv) < 0)
      { # If lower bound is smaller than zero, set estimate lower bound to zero
        tau.lb <- 0
      } else 
      {
        tmp.lb <- suppressWarnings(try(uniroot(get_LR_tau, interval = c(max(0, tau.est-tau.ci[1]), 
                                                                        tau.est), yi = yi, 
                                               vi = vi, est = est, tau.est = tau.est, ycv = ycv)$root, 
                                       silent = TRUE))
        
        tau.lb <- ifelse(class(tmp.lb) == "try-error", NA, tmp.lb) # Return NA if lower bound could not be estimated
      }
      
      tmp.ub <- suppressWarnings(try(uniroot(get_LR_tau, interval = c(tau.est, tau.est+tau.ci[2]), 
                                             yi = yi, vi = vi, est = est, tau.est = tau.est, ycv = ycv)$root, 
                                     silent = TRUE))
      
      tau.ub <- ifelse(class(tmp.ub) == "try-error", NA, tmp.ub) # Return NA if upper bound could not be estimated
    }
  } else if (method == "P" | method == "LNP")
  {
   
    ### Starting values for root-finding
    bounds.int <- con$bounds.int
    tau.int <- con$tau.int
    est.ci <- con$est.ci
    tau.ci <- con$tau.ci
    tol <- con$tol
    max.iter <- con$max.iter
    verbose <- con$verbose
    
    tau.est <- 0 # Use tau=0 for first step
    est <- 0 # Use est=0 for first step
    stay <- TRUE # In order to stay in while loop
    i <- 0 # Counter for number of iterations
    est.max <- taus.max <- numeric(10) # Empty objects to examine whether estimates oscilate
    
    ### While loop for estimating effect size and tau
    while(stay) 
    {
      
      i <- i+1 # Counter for number of iterations
      
      if (i > max.iter-99)
      { # Check if root finding is not oscillating between two values
        ### Compute difference between new and old estimates
        tau.dif.new <- round(abs(tau.est-tau.old), 2)
        
        if (tau.dif.new == tau.dif)
        { ### If root finding is oscillating, add one half of the difference between old
          # and new values to the estimates
          tau.est <- tau.est+(abs(tau.old-tau.est))/2
        }
      }
      
      if (i > max.iter-100)
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
      
      if (class(est) == "try-error")
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
          
          if (class(est) == "try-error")
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
      
      if (class(tau.est) == "try-error")
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
      
      if (i > max.iter-10)
      { # Store last ten estimates of effect size and tau before maximum number 
        # of iterations is reached
        est.max[i-max.iter+10] <- est
        taus.max[i-max.iter+10] <- tau.est
      }
      
      if (i == max.iter)
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
          
          if (class(est) == "try-error")
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
      # is less than tol or max.iter equals maximum number of iterations
      stay <- ifelse(abs(tau.est-tau.old) < tol & abs(est-old) < tol | i == max.iter, 
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
      lb <- NA
      ub <- NA
      tau.lb <- NA
      tau.ub <- NA
    } else 
    {
      ### Estimate CI of est ###
      lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(est-con$est.ci[1], est),  
                        tau = tau.est, yi = yi, vi = vi, param = "est", 
                        ycv = ycv, method = method, val = "ci.lb", 
                        get_cv_P(length(yi)))$root, silent = TRUE))
      
      if (class(lb) == "try-error") 
      { # Check if lower bound could be estimated
        lb <- NA
      } 
      
      ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(est, est+con$est.ci[2]),  
                        tau = tau.est, yi = yi, vi = vi, param = "est", 
                        ycv = ycv, method = method, val = "ci.ub", 
                        get_cv_P(length(yi)))$root, silent = TRUE))
      
      if (class(ub) == "try-error") 
      { # Check if upper bound could be estimated
        ub <- NA
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
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+con$tau.ci[1]), 
                                est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                        silent = TRUE))
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else 
        { # Estimate lower and upper bound
          
          tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(max(0, tau.est-con$tau.ci[2]), tau.est), 
                                est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                        silent = TRUE))  
          
          if (class(tau.lb) == "try-error") 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          } 
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+con$tau.ci[1]), 
                                est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                        silent = TRUE))
          
          if (class(tau.ub) == "try-error") 
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
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+con$tau.ci[1]), 
                                est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                        silent = TRUE))
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else 
        { # Estimate lower and upper bound
          if (con$tau.ci[2] == 0)
          { # If user did not specify a value for search interval, search from 0 to tau.est
            tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est), est = est, 
                                  yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                  method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                          silent = TRUE))  
          } else 
          { # Estimate lower and upper bound
            tau.lb <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(max(0, tau.est-con$tau.ci[2]), tau.est), 
                                  est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                  method = method, val = "ci.lb", cv_P = get_cv_P(length(yi)))$root, 
                          silent = TRUE))  
          }
          
          if (class(tau.lb) == "try-error") 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          } 
          
          tau.ub <- suppressWarnings(try(uniroot(pdist_nsig, interval = c(0, tau.est+con$tau.ci[1]), 
                                est = est, yi = yi, vi = vi, param = "tau", ycv = ycv, 
                                method = method, val = "ci.ub", cv_P = get_cv_P(length(yi)))$root, 
                        silent = TRUE))
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          }
        }
      }
    }
  }
  
  return(data.frame(est = est, tau.est = tau.est, lb = lb, ub = ub, tau.lb = tau.lb, 
                    tau.ub = tau.ub))
}
