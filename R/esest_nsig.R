esest_nsig <- function(stval_tau, yi, vi, d.int, tau.int, ycv, method) {

  if (method == "ML")
  {

    stay <- TRUE # In order to stay in while loop
    tau.new <- stval_tau # Use starting value of tau for first step
    d.new <- -999 # Unreasonable estimate to force procedure to use at least two iterations
    i <- 0 # Counter for number of iterations

    ### While loop for optimizing profile likelihood functions
    while(stay) {

      i <- i+1 # Counter for number of iterations

      ### For next iteration new estimate becomes old estimate
      d.old <- d.new
      tau.old <- tau.new

      ### Optimize profile likelihood function of delta
      d.new <- optimize(ml_d_C, d.int, tau.new, yi, vi, ycv, maximum = TRUE)$maximum

      ### Optimize profile likelihood function of tau
      tau.new <- optimize(ml_tau_C, tau.int, d.new, yi, vi, ycv, maximum = TRUE)$maximum

      ### Print intermediate steps if requested
      if (verbose == TRUE)
      {
        cat("d.new = ", d.new, "tau.new = ", tau.new, fill = TRUE)
      }

      ### Stay in while loop till difference between previous and new estimates
      # is less than tol or max.iter equals maximum number of iterations
      stay <- ifelse(abs(tau.new-tau.old) < tol & abs(d.new-d.old) < tol | i == max.iter, FALSE, TRUE)

    }

    if (i == max.iter) { # If maximum number of iterations is reached return NA
      d.new <- NA
      tau.new <- NA
    }

    ### If estimates are equal to bounds of interval that was used for optimization
    # return NA except if estimate of tau is equal to zero
    if (any(round(d.new, 3) %in% d.int | round(tau.new, 3) %in% tau.int) & round(tau.new, 3) != 0)
    {
      tau.new <- NA
      d.new <- NA
    }






  }

  return(data.frame(d.est = d.new, tau.est = tau.new))

  # return(list(est = est, ci.lb = ci.lb, ci.ub = ci.ub, approx.est = max(approx.est),
  #             approx.ci.lb = max(approx.ci.lb), ext.lb = ext.lb, tr.q = tr.q))
}
