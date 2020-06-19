#' hybrid
#'
#' Function to statistically combine original studies and replications by means
#' of the hybrid methods and fixed-effect meta-analysis as described in van Aert
#' and van Assen (2018).
#'
#' @param moi A vector of group means for one-sample means for the original studies
#' @param mri A vector of group means for one-sample means for the replications
#' @param roi A vector of raw correlations for the original studies
#' @param rri A vector of raw correlations for the replications
#' @param noi A vector of sample sizes for one-sample means and correlations for 
#' original studies
#' @param nri A vector of sample sizes for one-sample means and correlations for 
#' replications
#' @param sdoi A vector of standard deviations for one-sample means for original 
#' studies
#' @param sdri A vector of standard deviations for one-sample means for replications 
#' @param mo1i A vector of means in group 1 for two-independent means for the 
#' original studies
#' @param mr1i A vector of means in group 1 for two-independent means for the 
#' replications
#' @param mo2i A vector of means in group 2 for two-independent means for the 
#' original studies
#' @param mr2i A vector of means in group 2 for two-independent means for the 
#' replications
#' @param no1i A vector of sample sizes in group 1 for two-independent means for 
#' the original studies
#' @param nr1i A vector of sample sizes in group 1 for two-independent means for 
#' the replications
#' @param no2i A vector of sample sizes in group 2 for two-independent means for 
#' the original studies
#' @param nr2i A vector of sample sizes in group 2 for two-independent means for 
#' the replications
#' @param sdo1i A vector of standard deviations in group 1 for two-independent
#' means for the original studies
#' @param sdr1i A vector of standard deviations in group 1 for two-independent
#' means for the replications
#' @param sdo2i A vector of standard deviations in group 2 for two-independent
#' means for the original studies
#' @param sdr2i A vector of standard deviations in group 2 for two-independent
#' means for the replications
#' @param toobs A vector of t-values for one-sample means and two-independent means 
#' for the original studies
#' @param trobs A vector of t-values for one-sample means and two-independent means 
#' for the replications
#' @param yoi A vector of standardized effect sizes of the original studies 
#' (see Details)
#' @param voi A vector of sampling variances belonging to the standardized effect
#' sizes of the original studies (see Details)
#' @param yri A vector of standardized effect sizes of the replications 
#' (see Details)
#' @param vri A vector of sampling variances belonging to the standardized effect
#' sizes of the replications (see Details)
#' @param alpha A numerical value specifying the alpha level as used in the original
#' study (default is 0.05, see Details).
#' @param side A character indicating whether the observed effect size of the
#' original study is in the right-tail of the distribution (i.e., positive) or
#' in the left-tail of the distribution (i.e., negative) (either \code{"right"}
#' or \code{"left"})
#' @param mi A vector of group means for one-sample means (deprecated, see Details)
#' @param ri A vector of raw correlations (deprecated, see Details)
#' @param ni A vector of sample sizes for one-sample means and correlations 
#' (deprecated, see Details)
#' @param sdi A vector of standard deviations for one-sample means (deprecated, 
#' see Details)
#' @param m1i A vector of means in group 1 for two-independent means (deprecated, 
#' see Details)
#' @param m2i A vector of means in group 2 for two-independent means (deprecated, 
#' see Details)
#' @param n1i A vector of sample sizes in group 1 for two-independent means 
#' (deprecated, see Details)
#' @param n2i A vector of sample sizes in group 2 for two-independent means 
#' (deprecated, see Details)
#' @param sd1i A vector of standard deviations in group 1 for two-independent
#' means (deprecated, see Details)
#' @param sd2i A vector of standard deviations in group 2 for two-independent
#' means (deprecated, see Details)
#' @param tobs A vector of t-values (deprecated, see Details)
#'
#' @details Three different effect sizes can be used as input for the
#' \code{hybrid} function: one-sample means, two-independent means, and raw
#' correlation coefficients. For each effect size measure, data of the original 
#' studies and replications have to be provided separately. For analyzing 
#' one-sample means, either the group means (\code{moi} and \code{mri}), standard 
#' deviations (\code{sdoi} and \code{sdri}), and sample sizes (\code{noi} and 
#' \code{nri}) for the original studies and replications or t-values (\code{toobs} 
#' and \code{trobs}) and sample sizes (\code{noi} and \code{nri}) have to be 
#' provided. For analyzing two-independent means, either the group means of group 
#' 1 (\code{mo1i} and \code{mr1i}) and group 2 (\code{mo1i} and \code{mr1i}), 
#' standard deviations of group 1 (\code{sdo1i} and \code{sdr1i}) and group 2 
#' (\code{sdo2i} and \code{sdr2i}), and sample sizes of group 1 (\code{no1i} and 
#' \code{nr1i}) and group 2 (\code{no2i} and \code{nr2i}) for the original studies 
#' and replications have to be provided. It is also possible to analyze 
#' two-independent means by providing t-values (\code{toobs} and \code{trobs}) in 
#' combination with sample sizes of group 1 (\code{no1i} and \code{nr1i}) and group
#' 2 (\code{no2i} and \code{nr2i}) for the original studies and replications. 
#' Correlation coefficients can also be analyzed by supplying the function with 
#' raw correlation coefficients (\code{roi} and \code{rri}) and sample sizes 
#' (\code{noi} and \code{nri}) of the original studies and replications. The 
#' \code{side} argument to specify whether the observed effect size of the
#' original study is in the right-tail of the distribution (i.e., positive) or
#' in the left-tail of the distribution should also be specified for every effect size 
#' measure.
#' 
#' It is also possible to specify the standardized effect sizes and its sampling
#' variances directly via the \code{yoi}, \code{yri}, \code{voi}, and \code{vri} 
#' arguments. However, extensive knowledge about computing standardized effect 
#' sizes and its sampling variances is required and specifying standardized effect 
#' sizes and sampling variances is not recommended to be used if the p-values in 
#' the primary studies are not computed with a z-test. In case the p-values in the 
#' original studies were computed with, for instance, a t-test, the p-values of a 
#' z-test and t-test do not exactly coincide and studies may be not statistically 
#' significant according to a z-test.
#' 
#' The hybrid methods assume that the original studies are statistically
#' significant, so original studies that are not statistically signifcant are 
#' discarded from the analysis. Furthermore, it is assumed that two-tailed 
#' hypothesis tests were conducted in the original studies. In case one-tailed 
#' hypothesis tests were conducted in the original studies, the alpha level has 
#' to be multiplied by two. For example, if one-tailed hypothesis tests were 
#' conducted with an alpha level of .05, an alpha of 0.1 has to be entered into 
#' the \code{hybrid} function.
#' 
#' \strong{Previous version}
#' 
#' The usage of a previous version of the \code{hybrid} function was more restricted.
#' Users could only apply the method to a single original study and replication. 
#' Before the addition of the extra functionality to also analyze multiple original 
#' studies and replications, data of the original study and replication were 
#' specified in vectors containing two elements with the first element being 
#' the data of the original study and the second one the data of the replication. 
#' In order to maintain backwards compatibility, it is still possible to analyze 
#' data like this by using the arguments \code{m1i, m2i, mi, ri, sd1i, sd2i, sdi, 
#' n1i, n2i, ni, tobs}. However, using the \code{hybrid} function in this way is 
#' now deprecated.
#'
#' @return
#' \item{k}{total number of effect sizes}
#' \item{krep}{number of effect sizes of replications}
#' \item{est.hy}{effect size estimate of hybrid method}
#' \item{ci.lb.hy}{lower bound of hybrid method's confidence interval}
#' \item{ci.ub.hy}{upper bound of hybrid method's confidence interval}
#' \item{x.hy}{test statistic of hybrid method's test of null-hypothesis of no
#' effect}
#' \item{pval.hy}{two-tailed p-value of hybrid method's test of null-hypothesis
#'  of no effect}
#' \item{measure}{effect size measure}
#' \item{est.hyr}{effect size estimate of hybridR method}
#' \item{ci.lb.hyr}{lower bound of hybridR method's confidence interval}
#' \item{ci.ub.hyr}{upper bound of hybridR method's confidence interval}
#' \item{x.hyr}{test statistic of hybridR method's test of null-hypothesis of
#'  no effect}
#' \item{pval.hyr}{two-tailed p-value of hybridR method's test of
#' null-hypothesis of no effect}
#' \item{pval.o}{two-tailed p-value of original study}
#' \item{est.hy0}{effect size estimate of hybrid0 method}
#' \item{ci.lb.hy0}{lower bound of hybrid0 method's confidence interval}
#' \item{ci.ub.hy0}{upper bound of hybrid0 method's confidence interval}
#' \item{x.hy0}{test statistic of hybrid0 method's test of null-hypothesis of no
#' effect}
#' \item{pval.hy0}{two-tailed p-value of hybrid0 method's test of
#' null-hypothesis of no effect}
#' \item{est.fe}{effect size estimate based on traditional fixed-effect
#' meta-analysis}
#' \item{se.fe}{standard error of effect size estimate based on traditional
#' fixed-effect meta-analysis}
#' \item{zval.fe}{test statistic of the null-hypothesis of no effect based on
#' traditional fixed-effect meta-analysis}
#' \item{pval.fe}{two-tailed p-value of the null-hypothesis of no effect based
#' on traditional fixed-effect meta-analysis}
#' \item{ci.lb.fe}{lower bound of confidence interval based on traditional
#' fixed-effect meta-analysis}
#' \item{ci.ub.fe}{upper bound of confidence interval based on
#' traditional fixed-effect meta-analysis}
#' \item{est.repl}{effect size estimate of replication}
#' \item{se.repl}{standard error of replication's effect size estimate}
#' \item{ci.lb.repl}{lower bound of replication's confidence interval}
#' \item{ci.ub.repl}{upper bound of replication's confidence interval}
#' \item{stat.repl}{test statistic of replication for testing null-hypothesis of
#'  no effect}
#' \item{pval.repl}{two-tailed p-value of replication for testing
#' null-hypothesis of no effect}
#'
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references van Aert, R. C. M., & van Assen, M. A. L. M. (2018). Examining
#' reproducibility in psychology: A hybrid method for statistically combining a
#' biased original study and replication. Behavior Research Methods, 50(4): 1515-1539.
#' doi:10.3758/s13428-017-0967-6
#'
#' @examples
#' ### Apply hybrid method to example on page 5 of van Aert and van Assen (2018).
#'
#' pval.o <- 0.03 # p-value original study
#' pval.r <- 0.3 # p-value replication
#'
#' no1i <- no2i <- 40 # Sample size per group in original study
#' nr1i <- nr2i <- 80 # Sample size per group in replication
#' toobs <- qt(pval.o/2, df = no1i+no2i-2, lower.tail = FALSE) # Observed t-values original study
#' trobs <- qt(pval.r/2, df = nr1i+nr2i-2, lower.tail = FALSE) # Observed t-values replication
#'
#' ### Apply hybrid method
#' hybrid(toobs = toobs, trobs = trobs, no1i = no1i, no2i = no2i, nr1i = nr1i, 
#' nr2i = nr2i, side = "right")
#' 
#' ### Apply hybrid method to two original studies and two replications
#' 
#' noi <- nri <- 50 # Sample size original studies and replicaitons
#' sdoi <- sdri <- 1 
#' sei <- sdoi/sqrt(50) # Standard error
#' 
#' ### Generate data
#' pso <- c(0.025/3, 0.025/3*2)
#' psr <- c(1/3, 1/3*2)
#' moi <- qnorm(pso, mean = 0, sd = sei, lower.tail = FALSE)
#' mri <- qnorm(psr, mean = 0, sd = sei, lower.tail = FALSE)
#' 
#' ### Apply hybrid method
#' hybrid(moi = moi, noi = noi, sdoi = sdoi, mri = mri, nri = nri, sdri = sdri, side = "right")
#'
#' @export

hybrid <- function(mo1i, mo2i, moi, roi, sdo1i, sdo2i, sdoi, no1i, no2i, noi, 
                   toobs, mr1i, mr2i, mri, rri, sdr1i, sdr2i, sdri, nr1i, nr2i, 
                   nri, trobs, m1i, m2i, mi, ri, sd1i, sd2i, sdi, n1i, n2i, ni, 
                   tobs, yoi, yri, voi, vri, alpha = 0.05, side) 
{
  
  if (!missing("mi") | !missing("tobs") | !missing("m1i") | !missing("ri"))
  { # If only two studies are supplied as input
    
    if (!missing("mi") & !missing("ni") & !missing("sdi")) 
    {
      measure <- "M"
      es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$ori <- c(1, 0) # Create variable to indicate that first study is an original study
      
      res.repl <- repl(es = es, mi = mi, sdi = sdi, ni = ni, measure = measure,
                       side = side)
      
    } else if (!missing("ni") & !missing("tobs")) 
    {
      measure <- "MT"
      es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$ori <- c(1, 0) # Create variable to indicate that first study is an original study
      
      res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
      
    } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
               !missing("sd1i") & !missing("sd2i")) 
    {
      measure <- "MD"
      es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i,
                      sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
      
      es$ori <- c(1, 0) # Create variable to indicate that first study is an original study
      
      res.repl <- repl(es = es, m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i,
                       sd1i = sd1i, sd2i = sd2i, measure = measure, side = side)
      
    } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) 
    {
      measure <- "MDT"
      es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2,
                      side = side, measure = measure)
      
      es$ori <- c(1, 0) # Create variable to indicate that first study is an original study
      
      res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
      
    } else if (!missing("ri") & !missing("ni")) 
    {
      measure <- "COR"
      es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side,
                      measure = measure)
      
      es$ori <- c(1, 0) # Create variable to indicate that first study is an original study
      
      res.repl <- repl(es = es, measure = measure, side = side)
    }
    
    if (es$pval[1] > alpha/2) 
    {
      stop("Original study is not statistically significant")
    }
    
    ### If data of more than one original study and replication are provided
  } else
  {
    if (!missing("moi") & !missing("noi") & !missing("sdoi") & !missing("mri") & 
        !missing("nri") & !missing("sdri")) 
    {
      measure <- "M"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(mi = moi, ni = noi, sdi = sdoi, alpha = alpha/2, 
                          side = side, measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(mi = mri, ni = nri, sdi = sdri, alpha = alpha/2, 
                          side = side, measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames
      
    } else if (!missing("noi") & !missing("toobs") & !missing("nri") & !missing("trobs")) 
    {
      measure <- "MT"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(ni = noi, tobs = toobs, alpha = alpha/2, 
                          side = side, measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(ni = nri, tobs = trobs, alpha = alpha/2, 
                          side = side, measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames
      
    } else if (!missing("mo1i") & !missing("mo2i") & !missing("no1i") & 
               !missing("no2i") & !missing("sdo1i") & !missing("sdo2i") & 
               !missing("mr1i") & !missing("mr2i") & !missing("nr1i") & 
               !missing("nr2i") & !missing("sdr1i") & !missing("sdr2i")) 
    {
      measure <- "MD"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(m1i = mo1i, m2i = mo2i, n1i = no1i, n2i = no2i, sd1i = sdo1i,
                          sd2i = sdo2i, alpha = alpha/2, side = side, measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(m1i = mr1i, m2i = mr2i, n1i = nr1i, n2i = nr2i, sd1i = sdr1i,
                          sd2i = sdr2i, alpha = alpha/2, side = side, measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames

    } else if (!missing("no1i") & !missing("no2i") & !missing("toobs") & 
               !missing("nr1i") & !missing("nr2i") & !missing("trobs")) 
    {
      measure <- "MDT"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(n1i = no1i, n2i = no2i, tobs = toobs, alpha = alpha/2,
                          side = side, measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(n1i = nr1i, n2i = nr2i, tobs = trobs, alpha = alpha/2,
                          side = side, measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames

    } else if (!missing("roi") & !missing("noi") & !missing("rri") & !missing("nri")) 
    {
      measure <- "COR"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(ri = roi, ni = noi, alpha = alpha/2, side = side,
                          measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(ri = rri, ni = nri, alpha = alpha/2, side = side,
                          measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames

    } else if (!missing("yoi") & !missing("voi") & !missing("yri") & 
               !missing("vri"))
    {
      measure <- "SPE"
      
      ### Compute effect size for original studies
      es.yoi <- escompute(yi = yoi, vi = voi, alpha = alpha/2, side = side, 
                          measure = measure)
      es.yoi$ori <- 1 # Create variable to indicate that it is an original study
      
      ### Compute effect size for replications
      es.yri <- escompute(yi = yri, vi = vri, alpha = alpha/2, side = side, 
                          measure = measure)
      es.yri$ori <- 0 # Create variable to indicate that it is a replication
      
      es <- rbind(es.yoi, es.yri) # Bind data frames
      
    }
    
    ### Select only statistically significant original studies
    es <- subset(es, es$ori == 0 | (es$ori == 1 & es$pval < alpha/2))
    
    ### Conduct fixed-effect meta-analysis
    res.repl <- repl(es = es, measure = measure, side = side)
  }
  
  k <- nrow(es) # Total number of effect sizes
  krep <- sum(es$ori == 0) # Number of replications
  
  ##############################################################################

  ### Apply hybrid method
  res1 <- hy(es = es, measure = measure, side = side)
  
  ### Hybrid^R method
  if (res.repl$pval.o < alpha/2) 
  { # Use results of hybrid if two-tailed p-value of original studies < alpha/2
    res2 <- data.frame(est.hyr = res1$est.hy, ci.lb.hyr = res1$ci.lb.hy,
                       ci.ub.hyr = res1$ci.ub.hy, x.hyr = res1$x.hy, 
                       pval.hyr = res1$pval.hy, pval.o = res.repl$pval.o)
  } else 
  { # Use results of only replications if two-tailed p-value of orginal studies 
    # > alpha/2
    res2 <- data.frame(est.hyr = res.repl$est.repl, ci.lb.hyr = res.repl$ci.lb.repl,
                       ci.ub.hyr = res.repl$ci.ub.repl, x.hyr = res.repl$stat.repl,
                       pval.hyr = res.repl$pval.repl, pval.o = res.repl$pval.o)
  }
  
  ### Hybrid^0 method
  res3 <- hy0(es = es, res1 = res1, alpha = alpha/2)
  
  ### Apply fixed-effect meta-analysis
  res4 <- fe_ma(yi = es$yi, vi = es$vi)
  
  ### Transform results of fixed-effect meta-analysis
  if (measure == "COR") 
  { # Back transform Fisher z to correlation
    est.fe <- (exp(2*res4$est.fe) - 1)/(exp(2*res4$est.fe) + 1)
    se.fe <- (exp(2*res4$se.fe) - 1)/(exp(2*res4$se.fe) + 1)
    ci.lb.fe <- (exp(2*res4$ci.lb.fe) - 1)/(exp(2*res4$ci.lb.fe) + 1)
    ci.ub.fe <- (exp(2*res4$ci.ub.fe) - 1)/(exp(2*res4$ci.ub.fe) + 1)
    zval.fe <- res4$zval.fe
    
    if (side == "left") 
    { # Re-mirror estimates
      est.fe <- est.fe*-1
      tmp <- ci.ub.fe
      ci.ub.fe <- ci.lb.fe*-1
      ci.lb.fe <- tmp*-1
      zval.fe <- zval.fe*-1
    }
    
  } else if (side == "left" & measure != "COR") 
  { # Re-mirror estimates
    est.fe <- res4$est.fe*-1
    tmp <- res4$ci.ub.fe
    ci.ub.fe <- res4$ci.lb.fe*-1
    ci.lb.fe <- tmp*-1
    zval.fe <- res4$zval.fe*-1
  } else 
  {
    est.fe <- res4$est.fe
    ci.ub.fe <- res4$ci.ub.fe
    ci.lb.fe <- res4$ci.lb.fe
    zval.fe <- res4$zval.fe
  }
  
  x <- list(est.hy = res1$est, ci.lb.hy = res1$ci.lb, ci.ub.hy = res1$ci.ub,
            x.hy = res1$x, pval.hy = res1$pval, measure = measure, est.hyr = res2$est.hyr,
            ci.lb.hyr = res2$ci.lb.hyr, ci.ub.hyr = res2$ci.ub.hyr, x.hyr = res2$x.hyr,
            pval.hyr = res2$pval.hyr, pval.o = res2$pval.o, est.hy0 = res3$est.hy0,
            ci.lb.hy0 = res3$ci.lb.hy0, ci.ub.hy0 = res3$ci.ub.hy0, x.hy0 = res3$x.hy0,
            pval.hy0 = res3$pval.hy0, est.fe = est.fe, se.fe = res4$se.fe,
            ci.lb.fe = ci.lb.fe, ci.ub.fe = ci.ub.fe, zval.fe = zval.fe,
            pval.fe = res4$pval.fe, est.repl = res.repl$est.repl, se.repl = res.repl$se.repl,
            ci.lb.repl = res.repl$ci.lb.repl, ci.ub.repl = res.repl$ci.ub.repl,
            stat.repl = res.repl$stat.repl, pval.repl = res.repl$pval.repl, 
            k = k, krep = krep)
  class(x) <- "hybridoutput"
  return(x)
}