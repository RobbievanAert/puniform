# Changes in puniform 0.2.6 (XXXX-XX-XX)

- Fix a bug in hybrid() where the standard error was also transformed from Fisher's
z-value to Pearson correlation

- Improve pdist_hy() such that logarithm of the probabilities are computed. The 
approximation of extreme tail probabilities is now no longer needed.

- hybrid() now has a control argument that provide the user more control about the
estimation. Note that this may cause minor differences in estimates compared to 
the previous version due to different bounds that are used for estimation.

- The internal function bounds_hy() is no longer used for determining the bounds
for estimation and is removed from the package

- Improved the output of print.puni_staroutput() using the newly created report
function

- Added moderator analyses to puni_star()

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