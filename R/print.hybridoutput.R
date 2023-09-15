#' @export

print.hybridoutput <- function(x, ...) {
  
  if (x$con$implementation == "two")
  {
    cat("\n")
    cat("Results Hybrid method (k = ", x$k, "; Number conventional studies = ", x$k.conventional, ")", sep = "")
    cat("\n")
    cat("\n")
    x$pval.0 <- ifelse(x$pval.0 < 0.001, "  <.001", round(x$pval.0, 4))
    print(format(data.frame(est = round(x$est, 4), x = round(x$L.0, 4), pval = x$pval.0, 
                            ci.lb = round(x$ci.lb, 4), ci.ub = round(x$ci.ub, 4),
                            row.names = ""), width = 9))
    cat("\n")
    cat("===")
    cat("\n")
    cat("\n")
    cat("Results Hybrid^0 method")
    cat("\n")
    cat("\n")
    x$pval.hy0 <- ifelse(x$pval.hy0 < 0.001, "  <.001", round(x$pval.hy0, 4))
    print(format(data.frame(est = round(x$est.hy0, 4), x = round(x$L.0.hy0, 4), 
                            pval = round(x$pval.0.hy0, 4), ci.lb = round(x$ci.lb.hy0, 4), 
                            ci.ub = round(x$ci.ub.hy0, 4), row.names = ""), width = 9))
    cat("\n")
    cat("===")
    cat("\n")
    cat("\n")
    cat("Results Hybrid^R method")
    cat("\n")
    cat("\n")
    x$pval.hyr <- ifelse(x$pval.hyr < 0.001, "  <.001", round(x$pval.hyr, 4))
    x$pval.o <- ifelse(x$pval.o < 0.001, "<.001", round(x$pval.o, 4))
    print(format(data.frame(est = round(x$est.hyr, 4), x = round(x$L.0.hyr, 4), 
                            pval = round(x$pval.0.hyr, 4), ci.lb = round(x$ci.lb.hyr, 4), 
                            ci.ub = round(x$ci.ub.hyr, 4), row.names = ""), width = 9))
    cat("\n")
    cat("- Two-tailed p-value FE meta-analysis conventional studies:", x$pval.o)
    cat("\n")
    cat("\n")
    cat("===")
    cat("\n")
    cat("\n")
    cat("Fixed-effect meta-analysis")
    cat("\n")
    cat("\n")
    x$pval.fe <- ifelse(x$pval.fe < 0.001, "  <.001", round(x$pval.fe, 4))
    print(format(data.frame(est = round(x$est.fe, 4), se = round(x$se.fe, 4), 
                            zval = round(x$zval.fe, 4), pval = x$pval.fe, ci.lb = round(x$ci.lb.fe, 4),
                            ci.ub = round(x$ci.ub.fe, 4), row.names = ""), width = 9))
    cat("\n")
    cat("===")
    cat("\n")
    cat("\n")
    cat("Replications")
    cat("\n")
    cat("\n")
    x$pval.repl <- ifelse(x$pval.repl < 0.001, "  <.001", round(x$pval.repl, 4))
    print(format(data.frame(est = round(x$est.repl, 4), se = round(x$se.repl, 4),
                            stat = round(x$stat.repl, 4), pval = x$pval.repl, 
                            ci.lb = round(x$ci.lb.repl, 4), ci.ub = round(x$ci.ub.repl, 4), 
                            row.names = ""), width = 9))
    cat("\n")
  } else if (x$con$implementation == "multiple")
  {
    cat("\n")
    cat("Results Hybrid method (k = ", x$k, "; number conventional studies = ", x$k.conventional, ")", sep = "")
    cat("\n")
    cat("\n")
    cat("Model results Hybrid method:")
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
  }
  
} 