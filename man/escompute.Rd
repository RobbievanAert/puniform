% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/escompute.R
\name{escompute}
\alias{escompute}
\title{escompute}
\usage{
escompute(
  mi,
  ri,
  ni,
  sdi,
  m1i,
  m2i,
  n1i,
  n2i,
  sd1i,
  sd2i,
  tobs,
  yi,
  vi,
  alpha,
  side,
  measure
)
}
\arguments{
\item{mi}{A vector of group means for one-sample mean}

\item{ri}{A vector of raw correlation coefficients}

\item{ni}{A vector of sample sizes for one-sample mean}

\item{sdi}{A vector of standard deviations for one-sample mean}

\item{m1i}{A vector of means in group 1 for two-independent means}

\item{m2i}{A vector of means in group 2 for two-independent means}

\item{n1i}{A vector of sample sizes in group 1 for two-independent means}

\item{n2i}{A vector of sample sizes in group 2 for two-independent means}

\item{sd1i}{A vector of standard deviations in group 1 for two-independent means}

\item{sd2i}{A vector of standard deviations in group 2 for two-independent means}

\item{tobs}{A vector of t-values}

\item{yi}{A vector of standardized effect sizes}

\item{vi}{A vector of sampling variances belonging to the standardized effect sizes (\code{yi})}

\item{alpha}{A numerical value specifying the alpha level as used in primary studies}

\item{side}{A character indicating the direction of the tested hypothesis in the primary studies (either "\code{right}" or "\code{left}")}

\item{measure}{A character indicating what kind of effect size should be computed (Hedges' g or Fisher's r-to-z transformed correlation coefficients) and which arguments are used as input ("\code{M}", "\code{MT}", "\code{MD}", "\code{MDT}", or "\code{COR}"). See Details below.}
}
\value{
Function returns a data frame with standardized effect sizes (yi), variances of these standardized effect sizes (vi), z-values (zval), p-values as computed in primary studies (pval), and critical z-values (zcv).
}
\description{
Function that computes Hedges' g and its sampling variance for an one-sample mean and a two-independent means, Fisher's r-to-z transformed correlation coefficient and its sampling variance for a raw correlation coefficient and computes a p-value as in the primary studies was done.
}
\details{
The \code{measure} argument has to be used to specify the desired effect size and what input parameters are used. There are six options:
\itemize{
\item{\code{"M"}}{ for one-sample mean with \code{mi}, \code{ni}, \code{sdi}, \code{alpha}, and \code{side} as input parameters}
\item{\code{"MT"}}{ for one-sample mean with \code{tobs}, \code{ni}, \code{alpha}, and \code{side} as input parameters}
\item{\code{"MD"}}{ for two-sample mean with \code{m1i}, \code{m2i}, \code{n1i}, \code{n2i}, \code{sd1i}, \code{sd2i}, \code{alpha}, and \code{side} as input parameters}
\item{\code{"MDT"}}{ for two-sample mean with \code{tobs}, \code{n1i}, \code{n2i}, \code{alpha}, and \code{side} as input parameters}
\item{\code{"COR"}}{ for raw correlation coefficients with \code{ri}, \code{ni}, \code{alpha}, and \code{side} as input parameters}
\item{\code{"SPE"}}{ for user-specified standardized effect sizes and sampling variances with \code{yi}, \code{vi}, \code{alpha}, and \code{side} as input parameters}
}
}
\author{
Robbie C.M. van Aert \email{R.C.M.vanAert@tilburguniversity.edu}
}
