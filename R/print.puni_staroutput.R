#' @export

print.puni_staroutput <- function(x, ...) {
  cat("\n")
  cat("Method: ", x$method, " (k = ", x$k, "; ksig = ", x$ksig, ")", sep = "")
  cat("\n")
  cat("\n")
  cat("Model results p-uniform*:")
  cat("\n")
  cat("\n")
  if (x$con$type == "Wald" | x$con$type == "Wald/profile") 
  {
    print(format(data.frame(est = report(x$est), 
                            se = report(x$se[1:(length(x$se)-1)]), 
                            ci.lb = report(x$ci.lb),
                            ci.ub = report(x$ci.ub), 
                            zval = report(x$L.0),
                            pval = report(x$pval.0, type = "p"), 
                            row.names = x$var_names), width = 9))
  } else if (x$con$type == "profile")
  {
    print(format(data.frame(est = report(x$est), 
                            se = report(x$se[1:(length(x$se)-1)]), 
                            ci.lb = report(x$ci.lb),
                            ci.ub = report(x$ci.ub), 
                            LR = report(x$L.0),
                            pval = report(x$pval.0, type = "p"), 
                            row.names = x$var_names), width = 9))
  }
  if (x$method == "P" | x$method == "LNP")
  {
    cat("\n")
    cat("Note:")
    cat("\n")
    cat("- Test of no effect is not available for method", x$method)
    cat("\n")
  }
  cat("\n")
  cat("===")
  cat("\n")
  cat("\n")
  if (length(x$est) == 1)
  {
    cat("Estimating between-study variance:")
  } else if (length(x$est) > 1)
  {
    cat("Estimating residual between-study variance:")
  }
  cat("\n")
  cat("\n")
  if (x$con$type == "Wald")
  {
    print(format(data.frame(tau2 = report(x$tau2, type = "tau2"),
                            se = report(x$se[length(x$se)]), 
                            tau2.lb = report(x$tau2.lb, type = "tau2"),
                            tau2.ub = report(x$tau2.ub, type = "tau2"), 
                            zval = report(x$L.het),
                            pval = report(x$pval.het, type = "p"), 
                            row.names = ""), width = 9))
  } else if (x$con$type == "profile" | x$con$type == "Wald/profile")
  {
    print(format(data.frame(tau2 = report(x$tau2, type = "tau2"),
                            se = report(x$se[length(x$se)]), 
                            tau2.lb = report(x$tau2.lb, type = "tau2"),
                            tau2.ub = report(x$tau2.ub, type = "tau2"), 
                            LR = report(x$L.het),
                            pval = report(x$pval.het, type = "p"), 
                            row.names = ""), width = 9))
  }
  # cat("\n")
  # cat("===")
  # cat("\n")
  # cat("\n")
  # cat("Publication bias test p-uniform*")
  # cat("\n")
  # cat("\n")
  # x$pval.pb <- ifelse(x$pval.pb < 0.001, "  <.001", round(x$pval.pb, 4))
  # print(format(data.frame(L.pb = round(x$L.pb, 4), pval = x$pval.pb,
  #                         row.names = ""), width = 9))
  # if (x$method == "P" | x$method == "LNP")
  # {
  #   cat("\n")
  #   cat("Note:")
  #   cat("\n")
  #   cat("- Publication bias test is not yet implemented for method", x$method)
  #   cat("\n")
  # }
}