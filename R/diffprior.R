#' diffprior
#'
#' Function for computing posterior probabilities based on the Snapshot Bayesian
#' Hybrid Meta-Analysis Method or Snapshot Bayesian Meta-Analysis Method (uncorrected)
#' for another than a uniform prior.
#'
#' @param prior A vector of length four containing the prior probabilities for no,
#' small, medium, and large true effect size.
#' @param res.snap A data frame with posterior probabilities obtained with the
#' \code{snapshot} or \code{uncor.snapshot} function.
#'
#' @details The function computes posterior probabilities for four true effect
#' sizes (no, small, medium, and large) based on the \code{snapshot} or
#' \code{uncor.snapshot} function for another than a uniform prior. For more
#' information see van Aert and van Assen (2016).
#'
#' @return The \code{diffprior} function returns a data frame with posterior probabilities
#' for no (\code{p.0}), small (\code{p.sm}), medium (\code{p.me}), and large (\code{p.la})
#' true effect size.
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R.C.M. & van Assen, M.A.L.M. (2017). Bayesian evaluation
#' of effect size after replicating an original study. PLoS ONE, 12(4), e0175302. 
#' doi:10.1371/journal.pone.0175302
#'
#' @examples ### Example as presented on page 491 in Maxwell, Lau, and Howard (2015)
#' res.snap <- snapshot(ri=c(0.243, 0.114), ni=c(80, 172), alpha=.05)
#'
#' ### Prior probabilities with probablity for no effect twice as large as for the other true effects
#' prior <- c(0.4, 0.2, 0.2, 0.2)
#'
#' ### Compute posterior probabilities based on new prior
#' diffprior(prior = prior, res.snap = res.snap)
#'
#' @export
#'
diffprior <- function(prior, res.snap)
{

  p.0 <- (prior[1]*res.snap$p.0)/(prior[1]*res.snap$p.0+prior[2]*res.snap$p.sm+
                                    prior[3]*res.snap$p.me+prior[4]*res.snap$p.la)

  p.sm <- (prior[2]*res.snap$p.sm)/(prior[1]*res.snap$p.0+prior[2]*res.snap$p.sm+
                                      prior[3]*res.snap$p.me+prior[4]*res.snap$p.la)

  p.me <- (prior[3]*res.snap$p.me)/(prior[1]*res.snap$p.0+prior[2]*res.snap$p.sm+
                                      prior[3]*res.snap$p.me+prior[4]*res.snap$p.la)

  p.la <- (prior[4]*res.snap$p.la)/(prior[1]*res.snap$p.0+prior[2]*res.snap$p.sm+
                                      prior[3]*res.snap$p.me+prior[4]*res.snap$p.la)

  return(data.frame(p.0=p.0, p.sm=p.sm, p.me=p.me, p.la=p.la))
}
