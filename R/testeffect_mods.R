### Function to test null hypothesis of no effect in case of moderators
testeffect_mods <- function(es, mods, n_bs, ll, est, se, type, con)
{
  
  if (type == "profile")
  { # Likelihood-ratio test
    ll0 <- numeric(n_bs)
    
    for (b in 1:n_bs)
    {
      par_fixed <- rep(NA, n_bs+1)
      par_fixed[b] <- 0
      
      ll0[b] <- esest_mods(es = es, mods = mods, n_bs = n_bs, par_fixed = par_fixed, 
                           type = type, con = con)$ll
    }
    
    ### Conduct likelihood-ratio test
    L.0 <- -2*(ll0-ll)
    pval.0 <- pchisq(L.0, df = 1, lower.tail = FALSE)
    
  } else if (type == "Wald")
  { # Wald test
    L.0 <- est/se[1:n_bs]
    pval.0 <- 2*pnorm(abs(L.0), lower.tail = FALSE)
  }
  
  return(data.frame(L.0 = L.0, pval.0 = pval.0))
}