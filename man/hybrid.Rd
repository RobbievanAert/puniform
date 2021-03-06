% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hybrid.R
\name{hybrid}
\alias{hybrid}
\title{hybrid}
\usage{
hybrid(
  mo1i,
  mo2i,
  moi,
  roi,
  sdo1i,
  sdo2i,
  sdoi,
  no1i,
  no2i,
  noi,
  toobs,
  mr1i,
  mr2i,
  mri,
  rri,
  sdr1i,
  sdr2i,
  sdri,
  nr1i,
  nr2i,
  nri,
  trobs,
  m1i,
  m2i,
  mi,
  ri,
  sd1i,
  sd2i,
  sdi,
  n1i,
  n2i,
  ni,
  tobs,
  yoi,
  yri,
  voi,
  vri,
  alpha = 0.05,
  side
)
}
\arguments{
\item{mo1i}{A vector of means in group 1 for two-independent means for the 
original studies}

\item{mo2i}{A vector of means in group 2 for two-independent means for the 
original studies}

\item{moi}{A vector of group means for one-sample means for the original studies}

\item{roi}{A vector of raw correlations for the original studies}

\item{sdo1i}{A vector of standard deviations in group 1 for two-independent
means for the original studies}

\item{sdo2i}{A vector of standard deviations in group 2 for two-independent
means for the original studies}

\item{sdoi}{A vector of standard deviations for one-sample means for original 
studies}

\item{no1i}{A vector of sample sizes in group 1 for two-independent means for 
the original studies}

\item{no2i}{A vector of sample sizes in group 2 for two-independent means for 
the original studies}

\item{noi}{A vector of sample sizes for one-sample means and correlations for 
original studies}

\item{toobs}{A vector of t-values for one-sample means and two-independent means 
for the original studies}

\item{mr1i}{A vector of means in group 1 for two-independent means for the 
replications}

\item{mr2i}{A vector of means in group 2 for two-independent means for the 
replications}

\item{mri}{A vector of group means for one-sample means for the replications}

\item{rri}{A vector of raw correlations for the replications}

\item{sdr1i}{A vector of standard deviations in group 1 for two-independent
means for the replications}

\item{sdr2i}{A vector of standard deviations in group 2 for two-independent
means for the replications}

\item{sdri}{A vector of standard deviations for one-sample means for replications}

\item{nr1i}{A vector of sample sizes in group 1 for two-independent means for 
the replications}

\item{nr2i}{A vector of sample sizes in group 2 for two-independent means for 
the replications}

\item{nri}{A vector of sample sizes for one-sample means and correlations for 
replications}

\item{trobs}{A vector of t-values for one-sample means and two-independent means 
for the replications}

\item{m1i}{A vector of means in group 1 for two-independent means (deprecated, 
see Details)}

\item{m2i}{A vector of means in group 2 for two-independent means (deprecated, 
see Details)}

\item{mi}{A vector of group means for one-sample means (deprecated, see Details)}

\item{ri}{A vector of raw correlations (deprecated, see Details)}

\item{sd1i}{A vector of standard deviations in group 1 for two-independent
means (deprecated, see Details)}

\item{sd2i}{A vector of standard deviations in group 2 for two-independent
means (deprecated, see Details)}

\item{sdi}{A vector of standard deviations for one-sample means (deprecated, 
see Details)}

\item{n1i}{A vector of sample sizes in group 1 for two-independent means 
(deprecated, see Details)}

\item{n2i}{A vector of sample sizes in group 2 for two-independent means 
(deprecated, see Details)}

\item{ni}{A vector of sample sizes for one-sample means and correlations 
(deprecated, see Details)}

\item{tobs}{A vector of t-values (deprecated, see Details)}

\item{yoi}{A vector of standardized effect sizes of the original studies 
(see Details)}

\item{yri}{A vector of standardized effect sizes of the replications 
(see Details)}

\item{voi}{A vector of sampling variances belonging to the standardized effect
sizes of the original studies (see Details)}

\item{vri}{A vector of sampling variances belonging to the standardized effect
sizes of the replications (see Details)}

\item{alpha}{A numerical value specifying the alpha level as used in the original
study (default is 0.05, see Details).}

\item{side}{A character indicating whether the observed effect size of the
original study is in the right-tail of the distribution (i.e., positive) or
in the left-tail of the distribution (i.e., negative) (either \code{"right"}
or \code{"left"})}
}
\value{
\item{k}{total number of effect sizes}
\item{krep}{number of effect sizes of replications}
\item{est.hy}{effect size estimate of hybrid method}
\item{ci.lb.hy}{lower bound of hybrid method's confidence interval}
\item{ci.ub.hy}{upper bound of hybrid method's confidence interval}
\item{x.hy}{test statistic of hybrid method's test of null-hypothesis of no
effect}
\item{pval.hy}{two-tailed p-value of hybrid method's test of null-hypothesis
 of no effect}
\item{measure}{effect size measure}
\item{est.hyr}{effect size estimate of hybridR method}
\item{ci.lb.hyr}{lower bound of hybridR method's confidence interval}
\item{ci.ub.hyr}{upper bound of hybridR method's confidence interval}
\item{x.hyr}{test statistic of hybridR method's test of null-hypothesis of
 no effect}
\item{pval.hyr}{two-tailed p-value of hybridR method's test of
null-hypothesis of no effect}
\item{pval.o}{two-tailed p-value of original study}
\item{est.hy0}{effect size estimate of hybrid0 method}
\item{ci.lb.hy0}{lower bound of hybrid0 method's confidence interval}
\item{ci.ub.hy0}{upper bound of hybrid0 method's confidence interval}
\item{x.hy0}{test statistic of hybrid0 method's test of null-hypothesis of no
effect}
\item{pval.hy0}{two-tailed p-value of hybrid0 method's test of
null-hypothesis of no effect}
\item{est.fe}{effect size estimate based on traditional fixed-effect
meta-analysis}
\item{se.fe}{standard error of effect size estimate based on traditional
fixed-effect meta-analysis}
\item{zval.fe}{test statistic of the null-hypothesis of no effect based on
traditional fixed-effect meta-analysis}
\item{pval.fe}{two-tailed p-value of the null-hypothesis of no effect based
on traditional fixed-effect meta-analysis}
\item{ci.lb.fe}{lower bound of confidence interval based on traditional
fixed-effect meta-analysis}
\item{ci.ub.fe}{upper bound of confidence interval based on
traditional fixed-effect meta-analysis}
\item{est.repl}{effect size estimate of replication}
\item{se.repl}{standard error of replication's effect size estimate}
\item{ci.lb.repl}{lower bound of replication's confidence interval}
\item{ci.ub.repl}{upper bound of replication's confidence interval}
\item{stat.repl}{test statistic of replication for testing null-hypothesis of
 no effect}
\item{pval.repl}{two-tailed p-value of replication for testing
null-hypothesis of no effect}
}
\description{
Function to statistically combine original studies and replications by means
of the hybrid methods and fixed-effect meta-analysis as described in van Aert
and van Assen (2018).
}
\details{
Three different effect sizes can be used as input for the
\code{hybrid} function: one-sample means, two-independent means, and raw
correlation coefficients. For each effect size measure, data of the original 
studies and replications have to be provided separately. For analyzing 
one-sample means, either the group means (\code{moi} and \code{mri}), standard 
deviations (\code{sdoi} and \code{sdri}), and sample sizes (\code{noi} and 
\code{nri}) for the original studies and replications or t-values (\code{toobs} 
and \code{trobs}) and sample sizes (\code{noi} and \code{nri}) have to be 
provided. For analyzing two-independent means, either the group means of group 
1 (\code{mo1i} and \code{mr1i}) and group 2 (\code{mo1i} and \code{mr1i}), 
standard deviations of group 1 (\code{sdo1i} and \code{sdr1i}) and group 2 
(\code{sdo2i} and \code{sdr2i}), and sample sizes of group 1 (\code{no1i} and 
\code{nr1i}) and group 2 (\code{no2i} and \code{nr2i}) for the original studies 
and replications have to be provided. It is also possible to analyze 
two-independent means by providing t-values (\code{toobs} and \code{trobs}) in 
combination with sample sizes of group 1 (\code{no1i} and \code{nr1i}) and group
2 (\code{no2i} and \code{nr2i}) for the original studies and replications. 
Correlation coefficients can also be analyzed by supplying the function with 
raw correlation coefficients (\code{roi} and \code{rri}) and sample sizes 
(\code{noi} and \code{nri}) of the original studies and replications. The 
\code{side} argument to specify whether the observed effect size of the
original study is in the right-tail of the distribution (i.e., positive) or
in the left-tail of the distribution should also be specified for every effect size 
measure.

It is also possible to specify the standardized effect sizes and its sampling
variances directly via the \code{yoi}, \code{yri}, \code{voi}, and \code{vri} 
arguments. However, extensive knowledge about computing standardized effect 
sizes and its sampling variances is required and specifying standardized effect 
sizes and sampling variances is not recommended to be used if the p-values in 
the primary studies are not computed with a z-test. In case the p-values in the 
original studies were computed with, for instance, a t-test, the p-values of a 
z-test and t-test do not exactly coincide and studies may be not statistically 
significant according to a z-test.

The hybrid methods assume that the original studies are statistically
significant, so original studies that are not statistically signifcant are 
discarded from the analysis. Furthermore, it is assumed that two-tailed 
hypothesis tests were conducted in the original studies. In case one-tailed 
hypothesis tests were conducted in the original studies, the alpha level has 
to be multiplied by two. For example, if one-tailed hypothesis tests were 
conducted with an alpha level of .05, an alpha of 0.1 has to be entered into 
the \code{hybrid} function.

\strong{Previous version}

The usage of a previous version of the \code{hybrid} function was more restricted.
Users could only apply the method to a single original study and replication. 
Before the addition of the extra functionality to also analyze multiple original 
studies and replications, data of the original study and replication were 
specified in vectors containing two elements with the first element being 
the data of the original study and the second one the data of the replication. 
In order to maintain backwards compatibility, it is still possible to analyze 
data like this by using the arguments \code{m1i, m2i, mi, ri, sd1i, sd2i, sdi, 
n1i, n2i, ni, tobs}. However, using the \code{hybrid} function in this way is 
now deprecated.
}
\examples{
### Apply hybrid method to example on page 5 of van Aert and van Assen (2018).

pval.o <- 0.03 # p-value original study
pval.r <- 0.3 # p-value replication

no1i <- no2i <- 40 # Sample size per group in original study
nr1i <- nr2i <- 80 # Sample size per group in replication
toobs <- qt(pval.o/2, df = no1i+no2i-2, lower.tail = FALSE) # Observed t-values original study
trobs <- qt(pval.r/2, df = nr1i+nr2i-2, lower.tail = FALSE) # Observed t-values replication

### Apply hybrid method
hybrid(toobs = toobs, trobs = trobs, no1i = no1i, no2i = no2i, nr1i = nr1i, 
nr2i = nr2i, side = "right")

### Apply hybrid method to two original studies and two replications

noi <- nri <- 50 # Sample size original studies and replicaitons
sdoi <- sdri <- 1 
sei <- sdoi/sqrt(50) # Standard error

### Generate data
pso <- c(0.025/3, 0.025/3*2)
psr <- c(1/3, 1/3*2)
moi <- qnorm(pso, mean = 0, sd = sei, lower.tail = FALSE)
mri <- qnorm(psr, mean = 0, sd = sei, lower.tail = FALSE)

### Apply hybrid method
hybrid(moi = moi, noi = noi, sdoi = sdoi, mri = mri, nri = nri, sdri = sdri, side = "right")

}
\references{
van Aert, R. C. M., & van Assen, M. A. L. M. (2018). Examining
reproducibility in psychology: A hybrid method for statistically combining a
biased original study and replication. Behavior Research Methods, 50(4): 1515-1539.
doi:10.3758/s13428-017-0967-6
}
\author{
Robbie C.M. van Aert \email{R.C.M.vanAert@tilburguniversity.edu}
}
