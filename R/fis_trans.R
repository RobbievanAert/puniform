#' fis_trans
#'
#' Function for transforming raw correlation coefficients to Fisher-z transformed
#' correlation coefficients and vice versa.
#'
#' @param r A numerical value being a raw correlation coefficient
#' @param fis A numerical value being a Fisher-z transformed correlation coefficient
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @export

fis_trans <- function(r, fis) 
{
  if (!missing("r")) 
  {
    out <- 0.5*log((1+r)/(1-r))
  }
  
  if (!missing("fis")) 
  {
    out <- (exp(2*fis) - 1)/(exp(2*fis) + 1)
  }
  return(out)
}
