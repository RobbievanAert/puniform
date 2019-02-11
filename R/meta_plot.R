#' Meta-plot
#'
#' Function to create meta-plots for two-independent means, raw correlations, and 
#' odds ratios. See van Assen et al. (2019) for more information.
#'
#' @param m1i A vector of means in group 1 for two-independent means
#' @param m2i A vector of means in group 2 for two-independent means
#' @param n1i A vector of sample sizes in group 1 for two-independent means
#' @param n2i A vector of sample sizes in group 2 for two-independent means
#' @param sd1i A vector of standard deviations in group 1 for two-independent means
#' @param sd2i A vector of standard deviations in group 2 for two-independent means
#' @param gi A vector of Hedges' g values for two-independent means if group means 
#' and standard deviations are not available
#' @param vgi A vector of Hedges' g sampling variances for two-independent means 
#' if group means and standard deviations are not available
#' @param ri A vector of raw correlations
#' @param ni A vector of sample sizes if raw correlations are the effect size measure
#' @param ai A vector of frequencies in upper left cell of 2x2 frequency table
#' @param bi A vector of frequencies in upper right cell of 2x2 frequency table
#' @param ci A vector of frequencies in lower left cell of 2x2 frequency table
#' @param di A vector of frequencies in lower right cell of 2x2 frequency table
#' @param alpha A integer specifying the alpha level as used in primary studies
#' (default is 0.05 but see Details).
#' @param method_tau2 A character indicating the estimation method for the 
#' between-study variance in true effect size in the meta-analysis 
#' (default is \code{"PM"}, but see Details).
#' @param nr_lines A character indicating whether all primary study's effect sizes
#' (\code{"all"}, default) or a selection of primary study's effect sizes 
#' (\code{"summary"}) are plotted (see Details)
#' @param main A character indicating the title of the plot (default is no title)
#' @param cex.pch An integer to control the size of the points in the plot 
#' 
#' @details The \code{meta_plot} function assumes that two-tailed hypothesis tests 
#' were conducted in the primary studies. In case one-tailed hypothesis tests were 
#' conducted in the primary studies, the submitted \code{alpha} argument to the 
#' \code{meta_plot} function has to be multiplied by two. For example, if one-tailed 
#' hypothesis tests were  conducted with an alpha level of .05, an alpha of 0.1 
#' has to be submitted to the \code{meta_plot} function.
#' 
#' Different estimators can be used for estimating the between-study variance in 
#' true effect size. The default estimator is the Paule-Mandel estimator 
#' (Paule & Mandel, 1982), because this estimator was recommended in Veroniki 
#' et al. (2016) and Langan, Higgins, and Simmonds (2016). However, all estimators 
#' that are included in the \code{rma.uni} function of the \code{metafor} package 
#' can be used, because this function is called in the \code{meta_plot} function.
#' 
#' When \code{nr_lines = "summary"} is specified, the estimates of meta-analyses 
#' based on primary studies with sufficient statistical power are displayed.
#' Next to the estimate and 95\% confidence interval of the meta-analysis including 
#' all studies (leftmost), it shows these results for studies with sufficient 
#' statistical power (80\%) to detect a large true effect size (left vertical line), 
#' medium true effect size (middle), and small true effect size (right). Note 
#' that the summary meta-plot is just the meta-plot with many meta-analyses and 
#' confidence intervals left out, and keeping the leftmost meta-analysis and 
#' those just immediately to the right of the vertical lines. 
#' 
#' @return An invisibly returned data frame consisting of the submitted data and
#' \item{yi}{Standardized effect sizes used in the analyses}
#' \item{vi}{Sampling variances of the standardized effect sizes used in the analyses}
#' \item{est_cum}{Estimates of the cumulative meta-analyses}
#' \item{lb_cum}{Lower bounds of the 95\% confidence intervals of the cumulative meta-analyses}
#' \item{ub_cum}{Upper bounds of the 95\% confidence intervals of the cumulative meta-analyses}
#' \item{pub_est}{Estimates of cumulative meta-analyses based on Mill's ratios}
#' \item{info}{Information of a primary study (only for two-independent means)}
#' \item{stand_info}{Standardized information of a primary study (only for 
#' two-independent means)}
#' \item{preci}{Precision of a primary study (only for odds ratios)}
#' 
#' @author Robbie C.M. van Aert \email{R.C.M.vanAert@@tilburguniversity.edu}
#'
#' @references Langan, D., Higgins, J. P. T., & Simmonds, M. (2016). Comparative 
#' performance of heterogeneity variance estimators in meta-analysis: A review of 
#' simulation studies. Research Synthesis Methods, 8(2), 181-198. doi:10.1002/jrsm.1198
#' @references Veroniki, A. A., Jackson, D., Viechtbauer, W., Bender, R., Bowden, J., 
#' Knapp, G., . . . Salanti, G. (2016). Methods to estimate the between-study variance 
#' and its uncertainty in meta-analysis. Research Synthesis Methods, 7(1), 55-79. 
#' doi:10.1002/jrsm.1164
#' @references van Assen, ..., & van Aert (2019). The meta-plot. Manuscript in 
#' preparation.
#'
#' @examples ### Load data from meta-analysis by McCall and Carriger (1993)
#' data(data.mccall93)
#'
#' ### Create meta-plot
#' meta_plot(ri = data.mccall93$ri, ni = data.mccall93$ni)
#'
#' ### Create summary meta-plot
#' meta_plot(ri = data.mccall93$ri, ni = data.mccall93$ni, nr_lines = "summary")
#'
#' @export

meta_plot <- function(m1i, m2i, sd1i, sd2i, n1i, n2i, gi, vgi, ri, ni, ai, bi, 
                      ci, di, alpha = .05, method_tau2 = "PM", nr_lines = "all", 
                      main = "", cex.pch = 1)
{
  ### Necessary to pass R CMD check otherwise the following warning will be present:
  # no visible binding for global variable yi vi
  yi <- vi <- NULL 
  
  ### Assume that two-tailed tests are conducted but results in predicted 
  # direction are only reported
  alpha <- alpha/2 
  
  if ((!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") &
       !missing("sd1i") & !missing("sd2i")) | (!missing("gi") & !missing("vgi"))) 
  { # Check whether standardized mean differences were entered
    
    if (!missing("gi") & !missing("vgi"))
    {
      dat <- data.frame(yi = gi, vi = vgi, n1i = n1i, n2i = n2i)
      pval <- pnorm(dat$yi/sqrt(dat$vi), lower.tail = FALSE)
    } else 
    {
      ### Create data frame of data
      dat <- data.frame(m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, 
                        n2i = n2i)
      
      ### Compute standardized effect sizes
      dat <- metafor::escalc(m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, 
                             n2i = n2i, measure = "SMD", vtype = "UB", data = dat)
      
      ### Compute one-tailed p-values
      spool <- with(dat, sqrt(((n1i-1) * sd1i^2 + (n2i-1) * sd2i^2)/(n1i+n2i-2)))
      tval <- with(dat, (m1i-m2i)/sqrt(spool^2 * (1/n1i + 1/n2i)))
      pval <- pt(tval, df = n1i+n2i-2, lower.tail = FALSE)
    }
    
    ### Information per study
    dat$info <- with(dat, n1i*n2i/(n1i+n2i))
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pval < alpha)
    
    ### Square root of standardized information per study divided by 400 because 
    # 650*650/(650+650) = 325 and total information that can be obtained if maximum 
    # total sample size is 1300
    dat$stand_info <- sqrt(dat$info/325)
    
    ### Order data based on stand_info from large to small
    dat <- dat[order(dat$stand_info, decreasing = TRUE), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$stand_info, decreasing = TRUE))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis based 
    # on the most precise effect size and end with a meta-analysis based on all
    # included effect sizes
    dat$est_cum <- cum_dat$estimate
    dat$lb_cum <- cum_dat$ci.lb
    dat$ub_cum <- cum_dat$ci.ub
    
    ### Select the rows of dat such that there are no equal values of stand_info in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$stand_info)) == FALSE), ]
    
    ### Compute Mill's ratio (expected value of truncated normal distribution) based 
    # on sampling variance of each effect size 
    m <- dat_uniq$n1i + dat_uniq$n2i - 2 # Degrees of freedom
    J <- exp(lgamma(m/2) - log(sqrt(m/2)) - lgamma((m-1)/2)) # Hedges' g correction factor
    ex_vi <- (m*J^2)/((m-2)*dat_uniq$info) # Exact variance of g (See (23) in Viechtbauer, 2007)
    ev <- sqrt(ex_vi)*(1/alpha)*dnorm(qnorm(1-alpha))
    
    ### Cumulative meta-analysis based on Mill's ratios starting with all Mill's ratios  
    # and then repeatedly dropping the least precise Mill's ratio
    ev_cum <- function(i) 
    {
      sum(ev[1:(length(ev)+1-i)]*(1/ex_vi[1:(length(ev)+1-i)]))/(sum(1/ex_vi[1:(length(ev)+1-i)]))
    }
    dat_uniq$pub_est <- rev(sapply(1:length(ev), FUN = ev_cum, simplify = TRUE))
    
    if (prop_sig > 0.8)
    { # If cumulative meta-analyses based on Mill's ratio will be ploted, determine 
      # y-axis limits based on cumulative meta-analyses of data and Mill's ratio
      ylim <- round(c(min(dat_uniq$lb_cum, dat_uniq$pub_est)-0.1, 
                      max(dat_uniq$ub_cum, dat_uniq$pub_est)+0.1), 1)
    } else
    {
      ylim <- round(c(min(dat_uniq$lb_cum)-0.1, max(dat_uniq$ub_cum)+0.1), 1)
    }
    
    ### Make sure that y-axis always start at least at zero
    ylim[1] <- ifelse(ylim[1] > 0, 0, ylim[1])
    
    ### Group results of studies together if information is larger than 650*650/(650+650)
    if (any(dat_uniq$stand_info > sqrt((650*650/(650+650))/325)) == TRUE)
    {
      ind <- which(dat_uniq$stand_info > sqrt((650*650/(650+650))/325))
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[1,"stand_info"] <- 1
    }
    
    if (nr_lines == "all")
    {
      draw_plot_g(dat = dat_uniq, ylim = ylim, alpha = alpha, prop_sig = prop_sig, 
                  main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$stand_info >= sqrt((651*651/(651+651))/325)), n = 1), 
               tail(which(dat_uniq$stand_info > sqrt((64*64/(64+64))/325) & 
                            dat_uniq$stand_info <= sqrt((392*392/(392+392))/325)), n = 1), 
               tail(which(dat_uniq$stand_info > sqrt((26*26/(26+26))/325) & 
                            dat_uniq$stand_info <= sqrt((64*64/(64+64))/325)), n = 1),
               tail(which(dat_uniq$stand_info <= sqrt((26*26/(26+26))/325)), n = 1), 
               nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_g(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, 
                  prop_sig = prop_sig, main = main, cex.pch = cex.pch)
    }
    
    ### Text presenting percentage of studies particular statistical power
    txt1 <- round(sum(dat$stand_info <= sqrt((26*26/(26+26))/325))/length(dat$n1i)*100) # Less than 80% for detecting large effect
    txt2 <- round(sum(dat$stand_info > sqrt((26*26/(26+26))/325))/length(dat$n1i)*100) # More than 80% for detecting large effect
    txt3 <- round(sum(dat$stand_info >= sqrt((64*64/(64+64))/325))/length(dat$n1i)*100) # More than 80% for detecting medium effect
    txt4 <- round(sum(dat$stand_info >= sqrt((392*392/(392+392))/325))/length(dat$n1i)*100) # More than 80% for detecting small effect
    mtext(paste(txt1, sep = ""), side = 3, at = sqrt((26*26/(26+26))/325)/2, 
          line = -1, cex = par()$cex.lab)
    mtext(paste(txt2, sep = ""), side = 3, at = (sqrt((64*64/(64+64))/325)+
                                                   sqrt((26*26/(26+26))/325))/2, 
          line = -1, cex = par()$cex.lab)
    mtext(paste(txt3, sep = ""), side = 3, at = (sqrt((392*392/(392+392))/325)+
                                                   sqrt((64*64/(64+64))/325))/2, 
          line = -1, cex = par()$cex.lab)
    mtext(paste(txt4, sep = ""), side = 3, at = (sqrt((650*650/(650+650))/325)+
                                                   sqrt((392*392/(392+392))/325))/2, 
          line = -1, cex = par()$cex.lab)
    
    ### Text presenting percentage of statistically significant effect sizes
    mtext(paste("Sig. = ", round(prop_sig*100,1), sep = ""), side = 3, 
          at = sqrt((190*190/(190+190))/325), line = -3, cex = par()$cex.lab)
    
  } 
  else if (!missing("ri") & !missing("ni"))
  { # Check whether correlation coefficients were entered
    
    ### Create data frame of data
    dat <- data.frame(ri = ri, ni = ni)
    
    ### Compute Fisher-z transformed correlations and their sampling variances
    dat <- metafor::escalc(ri = ri, ni = ni, data = dat, measure = "ZCOR")
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pnorm(dat$yi/sqrt(dat$vi), lower.tail = FALSE) < alpha)
    
    ### Order data based on vi from large to small
    dat <- dat[order(dat$vi), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$vi))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis based 
    # on the most precise effect size and end with a meta-analysis based on all
    # included effect sizes
    dat$est_cum <- metafor::transf.ztor(cum_dat$estimate)
    dat$lb_cum <- metafor::transf.ztor(cum_dat$ci.lb)
    dat$ub_cum <- metafor::transf.ztor(cum_dat$ci.ub)
    
    ### Select the rows of dat such that there are no equal values of vi in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$vi)) == FALSE), ]
    
    ### Compute Mill's ratio (expected value of truncated normal distribution) based 
    # on sampling variance of each effect size 
    ev <- sqrt(dat_uniq$vi)*(1/alpha)*dnorm(qnorm(1-alpha))
    
    ### Cumulative meta-analysis based on Mill's ratios starting with all Mill's ratios  
    # and then repeatedly dropping the least precise Mill's ratio
    ev_cum <- function(i) 
    {
      sum(ev[1:(length(ev)+1-i)]*(1/dat_uniq$vi[1:(length(ev)+1-i)]))/(sum(1/dat_uniq$vi[1:(length(ev)+1-i)]))
    }
    dat_uniq$pub_est <- metafor::transf.ztor(rev(sapply(1:length(ev), FUN = ev_cum, simplify = TRUE)))
    
    if (prop_sig > 0.8)
    { # If cumulative meta-analyses based on Mill's ratio will be ploted, determine 
      # y-axis limits based on cumulative meta-analyses of data and Mill's ratio
      ylim <- round(c(min(dat_uniq$lb_cum, dat_uniq$pub_est)-0.1, 
                      max(dat_uniq$ub_cum, dat_uniq$pub_est)+0.1), 1)
    } else
    {
      ylim <- round(c(min(dat_uniq$lb_cum)-0.1, max(dat_uniq$ub_cum)+0.1), 1)
    }
    
    ### Make sure that y-axis always start at least at zero
    ylim[1] <- ifelse(ylim[1] > 0, 0, ylim[1])
    
    ### Group results of studies together if study's sample size is larger than 
    # 1293 (95% power for detecting a small effect)
    if (any(dat_uniq$ni > 1293) == TRUE) 
    {
      ind <- which(dat_uniq$ni > 1293)
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[1,"ni"] <- 1293
    }
    
    if (nr_lines == "all")
    {
      draw_plot_r(dat = dat_uniq, ylim = ylim, alpha = alpha, prop_sig = prop_sig, 
                  main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$ni > 782), n = 1), 
               tail(which(dat_uniq$ni > 84 & dat_uniq$ni <= 782), n = 1), 
               tail(which(dat_uniq$ni > 29 & dat_uniq$ni <= 84), n = 1), 
               tail(which(dat_uniq$ni <= 29), n = 1), nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_r(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, 
                  prop_sig = prop_sig, main = main, cex.pch = cex.pch)
    }
    
    ### Text presenting percentage of studies particular statistical power
    txt1 <- round(sum(dat$ni <= 29)/length(dat$ni)*100) # Less than 80% for detecting large effect
    txt2 <- round(sum(dat$ni > 29)/length(dat$ni)*100) # More than 80% for detecting large effect
    txt3 <- round(sum(dat$ni >= 84)/length(dat$ni)*100) # More than 80% for detecting medium effect
    txt4 <- round(sum(dat$ni >= 782)/length(dat$ni)*100) # More than 80% for detecting small effect
    mtext(paste(txt1, sep = ""), side = 3, at = sqrt(29)/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt2, sep = ""), side = 3, at = (sqrt(84)+sqrt(29))/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt3, sep = ""), side = 3, at = (sqrt(782)+sqrt(84))/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt4, sep = ""), side = 3, at = (sqrt(1300)+sqrt(782))/2, line = -1, cex = par()$cex.lab)
    
    ### Text presenting percentage of statistically significant effect sizes
    mtext(paste("Sig. = ", round(prop_sig*100,1), sep = ""), side = 3, 
          at = sqrt(335), line = -3, cex = par()$cex.lab)
    
  }
  
  ##############################################################################
  
  else if (!missing("ai") & !missing("bi") & !missing("ci") & !missing("di"))
  { # Check whether cell frequencies of 2x2 table were entered
    
    ### Create data frame of data
    dat <- data.frame(ai = ai, bi = bi, ci = ci, di = di)
    
    ### Compute log odds ratios and their sampling variances
    dat <- metafor::escalc(ai = ai, bi = bi, ci = ci, di = di, data = dat, 
                           measure = "OR")
    
    ### Compute precision (1/se)
    dat$preci <- 1/sqrt(dat$vi)
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pnorm(dat$yi/sqrt(dat$vi), lower.tail = FALSE) < alpha)
    
    ### Order data based on vi from large to small
    dat <- dat[order(dat$vi), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$vi))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis based 
    # on the most precise effect size and end with a meta-analysis based on all
    # included effect sizes
    dat$est_cum <- exp(cum_dat$estimate)
    dat$lb_cum <- exp(cum_dat$ci.lb)
    dat$ub_cum <- exp(cum_dat$ci.ub)
    
    ### Select the rows of dat such that there are no equal values of vi in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$vi)) == FALSE), ]
    
    ### Compute Mill's ratio (expected value of truncated normal distribution) based 
    # on sampling variance of each effect size 
    ev <- sqrt(dat_uniq$vi)*(1/alpha)*dnorm(qnorm(1-alpha))
    
    ### Cumulative meta-analysis based on Mill's ratios starting with all Mill's ratios  
    # and then repeatedly dropping the least precise Mill's ratio
    ev_cum <- function(i) 
    {
      sum(ev[1:(length(ev)+1-i)]*(1/dat_uniq$vi[1:(length(ev)+1-i)]))/
        (sum(1/dat_uniq$vi[1:(length(ev)+1-i)]))
    }
    dat_uniq$pub_est <- exp(rev(sapply(1:length(ev), FUN = ev_cum, simplify = TRUE)))
    
    if (prop_sig > 0.8)
    { # If cumulative meta-analyses based on Mill's ratio will be ploted, determine 
      # y-axis limits based on cumulative meta-analyses of data and Mill's ratio
      ylim <- round(c(min(dat_uniq$lb_cum, dat_uniq$pub_est)-0.1, 
                      max(dat_uniq$ub_cum, dat_uniq$pub_est)+0.1), 1)
    } else
    {
      ylim <- round(c(min(dat_uniq$lb_cum)-0.1, max(dat_uniq$ub_cum)+0.1), 1)
    }
    
    ### Make sure that y-axis always start at least at no effect
    ylim[1] <- ifelse(ylim[1] > 1, 1, ylim[1])
    
    ### Group results of studies together if study's precision is larger than 
    # 6.949 (95% power for detecting a small effect)
    if (any(dat_uniq$preci > 6.949) == TRUE) 
    {
      ind <- which(dat_uniq$preci > 6.949)
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[1,"preci"] <- 6.949
    }
    
    if (nr_lines == "all")
    {
      draw_plot_or(dat = dat_uniq, ylim = ylim, alpha = alpha, prop_sig = prop_sig, 
                   main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$preci > 5.4), n = 1), 
               tail(which(dat_uniq$preci > 2.252 & dat_uniq$preci <= 5.4), n = 1), 
               tail(which(dat_uniq$preci > 1.472 & dat_uniq$preci <= 5.4), n = 1), 
               tail(which(dat_uniq$preci <= 1.472), n = 1), nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_or(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, 
                   prop_sig = prop_sig, main = main, cex.pch = cex.pch)
    }
    
    ### Text presenting percentage of studies particular statistical power
    txt1 <- round(sum(dat$preci <= 1.472)/nrow(dat)*100) # Less than 80% for detecting large effect
    txt2 <- round(sum(dat$preci > 1.472)/nrow(dat)*100) # More than 80% for detecting large effect
    txt3 <- round(sum(dat$preci >= 2.252)/nrow(dat)*100) # More than 80% for detecting medium effect
    txt4 <- round(sum(dat$preci >= 5.4)/nrow(dat)*100) # More than 80% for detecting small effect
    mtext(paste(txt1, sep = ""), side = 3, at = (0.17+1.472)/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt2, sep = ""), side = 3, at = (1.472+2.252)/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt3, sep = ""), side = 3, at = (2.252+5.4)/2, line = -1, cex = par()$cex.lab)
    mtext(paste(txt4, sep = ""), side = 3, at = (5.4+6.949)/2, line = -1, cex = par()$cex.lab)
    
    ### Text presenting percentage of statistically significant effect sizes
    mtext(paste("Sig. = ", round(prop_sig*100,1), sep = ""), side = 3, 
          at = (2.252+5.4)/2, line = -3, cex = par()$cex.lab)
    
  }
  
  invisible(dat_uniq)
  
}