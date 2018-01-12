### Function for estimation with p-uniform*
esest_nsig <- function(yi, vi, int, tau.int, ycv, method, con) 
{
  
  if (method == "ML")
  {
    
    stay <- TRUE # In order to stay in while loop
    tau.est <- con$stval.tau # Use starting value of tau for first step
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
      est <- suppressWarnings(optimize(ml_est, con$int, tau.est, yi, vi, ycv, 
                                       maximum = TRUE)$maximum)
      
      ### Optimize profile likelihood function of tau
      # (suppressWarnings() in order to be able to specify wide search intervals)
      tau.est <- suppressWarnings(optimize(ml_tau, con$tau.int, est, yi, vi, ycv, 
                                           maximum = TRUE)$maximum)
      
      ### Print intermediate steps if requested
      if (con$verbose == TRUE)
      {
        cat("est = ", est, "tau.est = ", tau.est, fill = TRUE)
      }
      
      ### Stay in while loop till difference between previous and new estimates
      # is less than tol or max.iter equals maximum number of iterations
      stay <- ifelse(abs(tau.est-tau.old) < con$tol & abs(est-old) < con$tol | 
                       i == con$max.iter, FALSE, TRUE)
      
    }
    
    if (i == con$max.iter | any(round(est, 3) %in% con$int | round(tau.est, 3) %in% 
                                con$tau.int) & round(tau.est, 3) != 0) 
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
      tmp.lb <- suppressWarnings(try(uniroot(get_LR_est, interval = c(est-con$est.ci[1], est), 
                                             yi = yi, vi = vi, est = est, tau.est = tau.est, 
                                             ycv = ycv)$root, silent = TRUE))
      
      lb <- ifelse(class(tmp.lb) == "try-error", NA, tmp.lb) # Return NA if lower bound could not be estimated
      
      tmp.ub <- suppressWarnings(try(uniroot(get_LR_est, interval = c(est, est+con$est.ci[2]), 
                                             yi = yi, vi = vi, est = est, tau.est = tau.est, 
                                             ycv = ycv)$root, silent = TRUE))
      
      ub <- ifelse(class(tmp.ub) == "try-error", NA, tmp.ub) # Return NA if upper bound could not be estimated
      
      ### Estimation of CI tau
      if (get_LR_tau(prof.tau = 0, yi = yi, vi = vi, est = est, tau.est = tau.est, ycv = ycv) < 0)
      { # If lower bound is smaller than zero, set estimate lower bound to zero
        tau.lb <- 0
      } else {
        tmp.lb <- suppressWarnings(try(uniroot(get_LR_tau, interval = c(max(0, tau.est-con$tau.ci[1]), 
                                                                        tau.est), yi = yi, 
                                               vi = vi, est = est, tau.est = tau.est, ycv = ycv)$root, 
                      silent = TRUE))
        
        tau.lb <- ifelse(class(tmp.lb) == "try-error", NA, tmp.lb) # Return NA if lower bound could not be estimated
      }
      
      tmp.ub <- suppressWarnings(try(uniroot(get_LR_tau, interval = c(tau.est, tau.est+con$tau.ci[2]), 
                            yi = yi, vi = vi, est = est, tau.est = tau.est, ycv = ycv)$root, 
                    silent = TRUE))
      
      tau.ub <- ifelse(class(tmp.ub) == "try-error", NA, tmp.ub) # Return NA if upper bound could not be estimated
    }
  } else if (method == "P" | method == "LNP")
  {
    
    stay <- TRUE # In order to stay in while loop
    tau.est <- con$stval.tau # Use tau=0 for first step
    est <- -999 # Unreasonable estimate to force procedure to use at least two iterations
    i <- 0 # Counter for number of iterations
    int <- con$int # First interval that is used for estimating effect size
    tau.int <- con$tau.int # First interval that is used for estimating tau
    
    ### While loop for optimizing profile likelihood functions
    while(stay) 
    {
      
      i <- i+1 # Counter for number of iterations
      
      if (i > 2) 
      { # Check if root finding is not oscillating between two values
        ### Compute difference between new and old estimates
        est.dif.new <- round(abs(est-old), 2)
        tau.dif.new <- round(abs(tau.est-tau.old), 2)
        
        if (est.dif.new == est.dif & tau.dif.new == tau.dif) {
          ### If root finding is oscillating, add one half of the difference between old
          # and new values to the estimates
          tau.est <- tau.est+(abs(tau.old-tau.est))/runif(1, min = 2, max = 5)
        }
      }
      
      if (i > 1) 
      {
        ### Compute difference between new and old estimates
        est.dif <- round(abs(old-est), 2)
        tau.dif <- round(abs(tau.old-tau.est), 2)
        
        ### Compute interval for root finding
        int <- c(est-5, est+5)
        tau.int <- c(max(0, tau.est-0.7), tau.est+0.7)
      }
      
      ### For next iteration new estimate becomes old estimate
      old <- est
      tau.old <- tau.est
      
      ### Estimate effect size
      est <- try(bisect_est(pdist_nsig, lo = int[1], hi = int[2], tol = con$tol, 
                            tau.est = tau.est, yi = yi, vi = vi, param = "est", 
                            ycv = ycv, method = method, val = "es"), silent = TRUE)
      
      if (class(est) == "try-error")
      { # Return NAs and break out while loop if effect size could not be estimated
        est <- NA
        tau.est <- NA
        break
      }
      
      ### Print intermediate steps if requested
      if (con$verbose == TRUE) 
      { 
        cat("est = ", est, "tau.est = ", tau.est, fill = TRUE) 
      }
      
      ### Check if tau is larger than zero and otherwise break out while loop
      if (i == 1) 
      {
        tau0 <- round(pdist_nsig(est, 0, yi, vi, "tau", ycv, method, "es", cv_P = 0), 3)
        if (method == "P" & tau0 >= 0 | method == "LNP" & tau0 <= 0) 
        {
          tau.est <- 0
          est <- try(bisect_est(pdist_nsig, lo = int[1], hi = int[2], tol = con$tol, 
                                tau.est = tau.est, yi = yi, vi = vi, param = "est", 
                                ycv = ycv, method = method, val = "es"), silent = TRUE)
          
          if (class(est) == "try-error")
          { # Return NAs if effect size could not be estimated
            est <- NA
            tau.est <- NA
          }
          break
        }
      }
      
      ### Estimate tau
      tau.est <- try(bisect_tau(pdist_nsig, lo = tau.int[1], hi = tau.int[2], 
                                tol = con$tol, est = est, yi = yi, vi = vi, 
                                param = "tau", ycv = ycv, method = method, val = "es"), 
                     silent = TRUE)
      
      if (class(tau.est) == "try-error")
      { # Return NAs and break out while loop if tau could not be estimated
        est <- NA
        tau.est <- NA
        break
      }
      
      ### Stay in while loop till difference between previous and new estimates
      # is less than tol or max.iter equals maximum number of iterations
      stay <- ifelse(abs(tau.est-tau.old) < con$tol & abs(est-old) < con$tol | 
                       i == con$max.iter, FALSE, TRUE)
    }
    
    if (i == con$max.iter | any(round(est, 3) %in% int | round(tau.est, 3) %in% 
                                tau.int) & round(tau.est, 3) != 0) 
    { # If maximum number of iterations is reached or if estimates are equal to 
      # bounds of interval that was used for optimization return NA except if 
      # estimate of tau is equal to zero
      est <- NA
      tau.est <- NA
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
      lb <- try(bisect_est(pdist_nsig, lo = est-con$est.ci[1], hi = est, tol = con$tol, 
                           tau.est = tau.est, yi = yi, vi = vi, param = "est", 
                           ycv = ycv, method = method, val = "ci.lb"), silent = TRUE)
      
      if (class(lb) == "try-error") 
      { # Check if lower bound could be estimated
        lb <- NA
      } 
      
      ub <- try(bisect_est(pdist_nsig, lo = est, hi = est+con$est.ci[2], tol = con$tol, 
                           tau.est = tau.est, yi = yi, vi = vi, param = "est", 
                           ycv = ycv, method = method, val = "ci.ub"), silent = TRUE)
      
      if (class(ub) == "try-error") 
      { # Check if upper bound could be estimated
        ub <- NA
      } 
      
      ### Estimate CI of tau ###
      if (method == "P")
      {
        
        if (pdist_nsig(est, 0, yi, vi, "tau", ycv, method, "ci.ub", get_cv_P(length(yi))) > 0)
        {
          tau.lb <- tau.ub <- 0 # Return 0 (null set) if lower and upper bound are negative
        } else if (pdist_nsig(est, 0, yi, vi, "tau", ycv, method, "ci.lb", get_cv_P(length(yi))) > 0) 
        {
          tau.lb <- 0 # Truncate lower bound to zero if it is negative
          
          tau.ub <- try(bisect_tau(pdist_nsig, lo = tau.est, hi = tau.est+con$tau.ci[2], 
                                   tol = con$tol, est = est, yi = yi, vi = vi, 
                                   param = "tau", ycv = ycv, method = method, 
                                   val = "ci.ub"), silent = TRUE)
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else { # Estimate lower and upper bound
          tau.lb <- try(bisect_tau(pdist_nsig, lo = max(0, tau.est-con$tau.ci[1]), 
                                   hi = tau.est, tol = con$tol, est = est, yi = yi, 
                                   vi = vi, param = "tau", ycv = ycv,
                                   method = method, val = "ci.lb"), silent = TRUE)
          
          if (class(tau.lb) == "try-error") 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          } 
          
          tau.ub <- try(bisect_tau(pdist_nsig, lo = tau.est, hi = tau.est+con$tau.ci[2], 
                                   tol = con$tol, est = est, yi = yi, vi = vi, 
                                   param = "tau", ycv = ycv, method = method, 
                                   val = "ci.ub"), silent = TRUE)
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          }
        }
      } else if (method == "LNP")
      {
        ### Check if CI is a null set or estimates are too extreme to compute conditional probabilities at tau=0
        null.ub <- pdist_nsig(est, 0, yi, vi, "tau", ycv, method, "ci.ub", 0)
        null.lb <- pdist_nsig(est, 0, yi, vi, "tau", ycv, method, "ci.lb", 0)
        
        if (is.nan(null.ub) == FALSE & null.ub < 0) 
        {
          tau.lb <- tau.ub <- 0 # Return 0 (null set) if lower and upper bound are negative
        } else if (is.nan(null.lb) == FALSE & null.lb < 0) 
        {
          tau.lb <- 0 # Truncate lower bound to zero if it is negative
          
          tau.ub <- try(bisect_tau(pdist_nsig, lo = tau.est, hi = tau.est+con$tau.ci[2], 
                                   tol = con$tol, est = est, yi = yi, vi = vi, 
                                   param = "tau", ycv = ycv, method = method, 
                                   val = "ci.ub"), silent = TRUE)
          
          if (class(tau.ub) == "try-error") 
          { # Check if upper bound could be estimated
            tau.ub <- NA
          } 
          
        } else 
        { # Estimate lower and upper bound
          tau.lb <- try(bisect_tau(pdist_nsig, lo = max(0, tau.est-con$tau.ci[1]), 
                                   hi = tau.est, tol = con$tol, est = est, yi = yi, 
                                   vi = vi, param = "tau", ycv = ycv,
                                   method = method, val = "ci.lb"), silent = TRUE)
          
          if (class(tau.lb) == "try-error") 
          { # Check if lower bound could be estimated
            tau.lb <- NA
          }
          
          tau.ub <- try(bisect_tau(pdist_nsig, lo = tau.est, hi = tau.est+con$tau.ci[2], 
                                   tol = con$tol, est = est, yi = yi, vi = vi, param = "tau", 
                                   ycv = ycv, method = method, val = "ci.ub"), silent = TRUE)
          
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
