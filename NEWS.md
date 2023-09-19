# Changes in puniform 0.2.7 (2023-09-18)

- Updated and improved documentation

- Include a new implementation of the hybrid() function such that it also works
for multiple conventional and preregistered studies.

- The R package "metadat" is now included as suggested package, because the data
from Lehmann et al. (2018) are now used as example in the documentation of the
hybrid() function.


# Changes in puniform 0.2.6 (2023-07-14)

- Fix a bug in hybrid() where the standard error was also transformed from Fisher's
z-value to Pearson correlation

- Improve pdist_hy() such that logarithm of the probabilities are computed. The 
approximation of extreme tail probabilities is now no longer needed.

- hybrid() now has a control argument that provide the user more control about the
estimation. Note that this may cause minor differences in estimates compared to 
the previous version due to different bounds that are used for estimation.

- The internal function bounds_hy() is no longer used for determining the bounds
for estimation and is removed from the package

- The default optimization procedure of puni_star() and method = "ML" is now to 
optimize both parameters at the same time. The previous version where the profile
log-likelihood functions of both parameters were iteratively optimized is still 
available by setting the control argument proc.ml = "prof".

- Fix a bug in the computation of the profile likelihood confidence intervals of 
esest_nsig() and in the computation of the likelihood-ratio test in testeffect_nsig()
and testhetero(). Note that the functions get_LR_est() and get_LR_tau() are not
necessary anymore and are deleted from the package.

- Updated and improved documentation


# Changes in puniform 0.2.5 (2022-03-18)

- Modifying meta_plot() function to avoid a warning by the cumul() function of the 
metafor package

- Updated and improved documentation


# Changes in puniform 0.2.4 (2020-12-20)

- pub_bias argument in meta_plot() function is now more flexible

- Publication bias test of puni_star() removed for now

- Updated and improved documentation


# Changes in puniform 0.2.3 (2020-10-19)

- Added extra information to the plot created with the meta_plot() function

- Bug fix for conducting likelihood ratio tests with the puni_star() function


# Changes in puniform 0.2.2 (2020-06-19)

- Updated and improved documentation

- var_dif_fis(), var_boot_fis(), var_dif_rmd(), var_boot_rmd(), and var_pop() 
were added to estimate the variability of the outcomes' effect size. These 
functions can be used in a meta-regression model to correct for outcome reporting 
bias with the CORB method.


# Changes in puniform 0.2.1 (2019-08-23)

- meta_plot() function is added to package


# Changes in puniform 0.1.1 (2019-03-06)

- Updated and improved documentation

- Setting default parameters in puniform(), hybrid(), snapshot(), and puni_star() without
losing backwards compatibility


# Changes in puniform 0.1.0

- First release on CRAN