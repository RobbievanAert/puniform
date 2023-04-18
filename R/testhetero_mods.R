### Function to test null hypothesis of no effect in case of moderators
testhetero_mods <- function(es, mods, n_bs, ll, tau2.est, se, type, con)
{
  
  if (type == "profile")
  { # Likelihood-ratio test
    
    if (con$optimizer == "L-BFGS-B")
    {
      ll0 <- esest_mods(es = es, mods = mods, n_bs = n_bs, 
                        par_fixed = c(rep(NA, n_bs), 0), type = type, con = con)$ll  
    } else if (con$optimizer != "L-BFGS-B")
    {
      ll0 <- esest_mods(es = es, mods = mods, n_bs = n_bs, 
                        par_fixed = c(rep(NA, n_bs), log(0)), type = type, con = con)$ll  
    }
    
    ### Conduct likelihood-ratio test
    L.het <- -2*(ll0-ll)
    pval.het <- pchisq(L.het, df = 1, lower.tail = FALSE)
    
  } else if (type == "Wald")
  { # Wald test
    L.het <- tau2.est/se[length(se)]
    pval.het <- 2*pnorm(abs(L.het), lower.tail = FALSE)
  }
  
  return(data.frame(L.het = L.het, pval.het = pval.het))
}