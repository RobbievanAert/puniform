### Function for reporting results that rounds by default to four digits, shows 
# trailing zeros when desired, and uses <.001 for p-values if necessary
report <- function(x, digits = 4, format = "f", type = "number")
{
  out <- sapply(x, FUN = function(x) {
    
  if (is.na(x) == TRUE)
  {
    NA
  } else if (type == "number")
  {
    formatC(x, format = format, digits = digits)
  } else if (type == "p")
  {
    ifelse(x < 0.001, "<.001", 
                  sub("^(-?)0.", "\\1.", sprintf(paste0("%.", digits, "f"), x)))
  } else if (type == "tau2")
  {
    tmp <- formatC(x, format = format, digits = digits)
    
    ### Do not report 0.0000 but 0
    tmp <- ifelse(x == 0, "0", tmp)
  }
    
  })
  
  return(out)
}