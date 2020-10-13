#' ml_mods
#' 
#' @export

ml_mods <- function(par, yi, vi, xi, ycv)
{
  b0 <- par[1]
  b1 <- par[2]
  tau <- par[3]
  
  q.sig <- q.nsig <- numeric() # Empty object for storing conditional probabilities
  
  for (i in 1:length(yi)) 
  {
    
    if (yi[i] > ycv[i]) 
    { # Compute conditional probabilities for significant effect sizes
      
      q.sig[i] <- dnorm(yi[i], mean = b0+b1*xi[i], sd = sqrt(vi[i]+tau^2))/
        pnorm(ycv[i], mean = b0+b1*xi[i], sd = sqrt(vi[i]+tau^2), lower.tail = FALSE)
      
    } else if (yi[i] < ycv[i]) 
    { # Compute conditional probabilities for nonsignificant effect sizes
      
      q.nsig[i] <- dnorm(yi[i], mean = b0+b1*xi[i], sd = sqrt(vi[i]+tau^2))/
        pnorm(ycv[i], mean = b0+b1*xi[i], sd = sqrt(vi[i]+tau^2))
      
    }
    
  }
  
  ### Remove NAs of vectors with conditional probabilities
  q.sig <- q.sig[!is.na(q.sig)]
  q.nsig <- q.nsig[!is.na(q.nsig)]
  
  return(-sum(log(c(q.sig, q.nsig)))) # Compute log likelihood
  
}