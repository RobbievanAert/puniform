README for "puniform" package

The puniform package provides methods for estimating effect size and a
confidence interval, testing the null-hypothesis of no effect, 
and testing for publication bias while taking into account that the primary studies are statistically significant. The methods in the puniform package are  based on the statistical theory that the distribution of p-values is uniform conditional on the population effect size.

The puniform package contains two different methods: the p-uniform method and 
the hybrid method. The p-uniform method is a meta-analysis method that estimates
 the effect size by in a set of studies by searching for the effect size where 
the distribution of conditional p-values is uniform. The p-uniform method can 
also be used for creating a confidence interval around its effect size estimate,
testing the null-hypothesis of no effect, and testing for publication bias. The
hybrid method is a meta-analysis method for combining an original study and 
replication and while taking into account statistical significance of the  original study. 

Please note that this package is still under development and that this is a beta version. If you suspect a bug, please send me an email (R.C.M.vanAert@tilburguniversity.edu).

Papers about the p-uniform method:

van Assen, M. A. L. M., Van Aert, R. C. M., & Wicherts, J. M. (2015). Meta-analysis using effect size distributions of only statistically significant studies. Psychological Methods, 20(3), 293-309. doi: http://dx.doi.org/10.1037/met0000025

Van Aert, R. C. M., Wicherts, J. M., & Van Assen, M. A. L. M. (2015). Conducting meta-analyses on p-values: Reservations and recommendations for applying p-uniform and p-curve. Manuscript submitted for publication.

Paper about the hybrid method:

Van Aert, R. C. M., & Van Assen, M. A. L. M. (2016). Examining reproducibility 
in psychology: A hybrid method for statistically combining a biased original study and replication. Manuscript submitted for publication.
