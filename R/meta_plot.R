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
#' (see Details)
#' @param bi A vector of frequencies in upper right cell of 2x2 frequency table
#' (see Details)
#' @param ci A vector of frequencies in lower left cell of 2x2 frequency table
#' (see Details)
#' @param di A vector of frequencies in lower right cell of 2x2 frequency table
#' (see Details)
#' @param alpha A numerical value specifying the alpha level as used in primary studies
#' (default is 0.05 but see Details)
#' @param method_tau2 A character indicating the estimation method for the 
#' between-study variance in true effect size in the meta-analysis 
#' (default is \code{"PM"}, but see Details)
#' @param nr_lines A character indicating whether all primary study's effect sizes
#' (\code{"all"}, default) or a selection of primary study's effect sizes 
#' (\code{"summary"}) are plotted (see Details)
#' @param pub_bias A logical indicating whether the expected results of the cumulative 
#' meta-analysis based on a zero true effect in combination with extreme publication 
#' bias should be plotted (default is TRUE). Note that these results are only 
#' included if at least 80\% of the primary studies is statistically significant 
#' regardless of the \code{pub_bias} parameter  
#' @param main A character indicating the title of the plot (default is no title)
#' @param cex.pch A numerical value to control the size of the points in the plot 
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
#' For creating a meta-plot based on odds ratios as effect size measure, the 2x2 
#' frequency table should follow a specific format. The reason for this is that 
#' the probability for the outcome of interest in the control conditions has 
#' to be estimated. Hence, the 2x2 frequency table should look like this:
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
                      pub_bias = TRUE, main = "", cex.pch = 1)
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
    
    ### Get relative positions on the x-axis for small, medium, and large effect 
    # (80% power) and small effect (95% power)
    pos_95 <- (sqrt(650*650/(650+650))-sqrt(2))/(sqrt(325)-sqrt(2))
    pos_sm <- (sqrt(392*392/(392+392))-sqrt(2))/(sqrt(325)-sqrt(2))
    pos_me <- (sqrt(64*64/(64+64))-sqrt(2))/(sqrt(325)-sqrt(2))
    pos_la <- (sqrt(26*26/(26+26))-sqrt(2))/(sqrt(325)-sqrt(2))
    
    ### Get relative position on the x-axis of observerd effect sizes
    dat$posi <- (sqrt(dat$n1i*dat$n2i/(dat$n1i+dat$n2i))-sqrt(2))/(sqrt(325)-sqrt(2))
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pval < alpha)
    
    ### Order data based on stand_info from large to small
    dat <- dat[order(dat$posi, decreasing = TRUE), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$posi, decreasing = TRUE))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis 
    # based on the most precise effect size and end with a meta-analysis based on 
    # all included effect sizes
    dat$est_cum <- cum_dat$estimate
    dat$lb_cum <- cum_dat$ci.lb
    dat$ub_cum <- cum_dat$ci.ub
    
    ### Select the rows of dat such that there are no equal values of posi in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$posi)) == FALSE), ]
    
    ### Compute Mill's ratio (expected value of truncated normal distribution) based 
    # on sampling variance of each effect size 
    m <- dat_uniq$n1i + dat_uniq$n2i - 2 # Degrees of freedom
    J <- exp(lgamma(m/2) - log(sqrt(m/2)) - lgamma((m-1)/2)) # Hedges' g correction factor
    
    ### Exact variance of g (See (22) in Viechtbauer, 2007)
    ex_vi <- (m*J^2)/((m-2)*dat_uniq$n1i*dat_uniq$n2i/(dat_uniq$n1i+dat_uniq$n2i)) 
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
    if (any(dat_uniq$posi > 1) == TRUE)
    {
      ind <- which(dat_uniq$posi > 1)
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[1,"posi"] <- 1
    }
    
    if (nr_lines == "all")
    {
      draw_plot_g(dat = dat_uniq, ylim = ylim, alpha = alpha, pub_bias = pub_bias, 
                  prop_sig = prop_sig, pos_sm = pos_sm, pos_me = pos_me, 
                  pos_la = pos_la, main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$posi > pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_me & dat_uniq$posi <= pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_la & dat_uniq$posi <= pos_me), n = 1), 
               tail(which(dat_uniq$posi <= pos_sm), n = 1), nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_g(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, pub_bias = pub_bias, 
                  prop_sig = prop_sig, pos_sm = pos_sm, pos_me = pos_me, 
                  pos_la = pos_la, main = main, cex.pch = cex.pch)
    }
    
  } 
  else if (!missing("ri") & !missing("ni"))
  { # Check whether correlation coefficients were entered
    
    ### Create data frame of data
    dat <- data.frame(ri = ri, ni = ni)
    
    ### Compute Fisher-z transformed correlations and their sampling variances
    dat <- metafor::escalc(ri = ri, ni = ni, data = dat, measure = "ZCOR")
    
    ### Get relative positions on the x-axis for small, medium, and large effect 
    # (80% power) and small effect (95% power)
    pos_95 <- (sqrt(1293)-2)/(sqrt(1300)-2)
    pos_sm <- (sqrt(782)-2)/(sqrt(1300)-2)
    pos_me <- (sqrt(84)-2)/(sqrt(1300)-2)
    pos_la <- (sqrt(29)-2)/(sqrt(1300)-2)
    
    ### Get relative position on the x-axis of observerd effect sizes
    dat$posi <- (sqrt(dat$ni)-2)/(sqrt(1300)-2)
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pnorm(dat$yi/sqrt(dat$vi), lower.tail = FALSE) < alpha)
    
    ### Order data based on vi from large to small
    dat <- dat[order(dat$posi, decreasing = TRUE), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$posi, decreasing = TRUE))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis based 
    # on the most precise effect size and end with a meta-analysis based on all
    # included effect sizes
    dat$est_cum <- metafor::transf.ztor(cum_dat$estimate)
    dat$lb_cum <- metafor::transf.ztor(cum_dat$ci.lb)
    dat$ub_cum <- metafor::transf.ztor(cum_dat$ci.ub)
    
    ### Select the rows of dat such that there are no equal values of vi in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$posi)) == FALSE), ]
    
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
    if (any(dat_uniq$posi > 1) == TRUE) 
    {
      ind <- which(dat_uniq$posi > 1)
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[1,"posi"] <- 1
    }
    
    if (nr_lines == "all")
    {
      draw_plot_r(dat = dat_uniq, ylim = ylim, alpha = alpha, pub_bias = pub_bias, 
                  prop_sig = prop_sig, pos_sm = pos_sm, pos_me = pos_me, 
                  pos_la = pos_la, main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$posi > pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_me & dat_uniq$posi <= pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_la & dat_uniq$posi <= pos_me), n = 1), 
               tail(which(dat_uniq$posi <= pos_sm), n = 1), nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_r(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, pub_bias = pub_bias,
                  prop_sig = prop_sig, pos_sm = pos_sm, pos_me = pos_me, pos_la = pos_la, 
                  main = main, cex.pch = cex.pch)
    }
    
  }
  
  ##############################################################################
  
  else if (!missing("ai") & !missing("bi") & !missing("ci") & !missing("di"))
  { # Check whether cell frequencies of 2x2 table were entered
    
    ### Create data frame of data
    dat <- data.frame(ai = ai, bi = bi, ci = ci, di = di)
    
    ### Compute log odds ratios and their sampling variances
    dat <- metafor::escalc(ai = ai, bi = bi, ci = ci, di = di, data = dat, 
                           measure = "OR", to = "all")
    
    ### Conduct meta-analysis based on proportions in order to estimate pi_C
    ma_prop <- metafor::rma(xi = ai, ni = ai+bi, method = method_tau2, data = dat, 
                            measure = "PLO", to = "all")
    
    ### Back-transform to get median probability of outcome in control group
    pi_c <- metafor::transf.ilogit(ma_prop$b[1]) 
    
    ### Compute probability and odds ratio of outcome in control group if effect 
    # size is small
    pi_e <- pnorm(qnorm(pi_c) + 0.2)
    sm <- 1/((pi_c/(1-pi_c))/(pi_e/(1-pi_e)))
    
    ### Compute probability and odds ratio of outcome in control group if effect 
    # size is medium
    pi_e <- pnorm(qnorm(pi_c) + 0.5)
    me <- 1/((pi_c/(1-pi_c))/(pi_e/(1-pi_e)))
    
    ### Compute probability and odds ratio of outcome in control group if effect 
    # size is large
    pi_e <- pnorm(qnorm(pi_c) + 0.8)
    la <- 1/((pi_c/(1-pi_c))/(pi_e/(1-pi_e)))
    
    ### Function to get standard error given a particular effect size and power
    get_se <- function(se, es, pow)
    {
      pnorm(qnorm(.975, sd = se), mean = log(es), sd = se, lower.tail = FALSE)+
        pnorm(qnorm(.025, sd = se), mean = log(es), sd = se, lower.tail = TRUE)-pow
    }
    
    ### Get precision for 95% power and small effect and 80% power for small, 
    # medium, and large effect
    prec_95 <- 1/uniroot(get_se, interval = c(0, 100), es = sm, pow = 0.95)$root
    prec_sm <- 1/uniroot(get_se, interval = c(0, 100), es = sm, pow = 0.8)$root
    prec_me <- 1/uniroot(get_se, interval = c(0, 100), es = me, pow = 0.8)$root
    prec_la <- 1/uniroot(get_se, interval = c(0, 100), es = la, pow = 0.8)$root
    
    ### Get relative position on the x-axis
    pos_95 <- (prec_95 - 1/sqrt(4/2.5))/(prec_95 - 1/sqrt(4/2.5))
    pos_sm <- (prec_sm - 1/sqrt(4/2.5))/(prec_95 - 1/sqrt(4/2.5))
    pos_me <- (prec_me - 1/sqrt(4/2.5))/(prec_95 - 1/sqrt(4/2.5))
    pos_la <- (prec_la - 1/sqrt(4/2.5))/(prec_95 - 1/sqrt(4/2.5))
    
    ### Get relative position on the x-axis of observerd effect sizes
    dat$posi <- (1/sqrt(dat$vi) - 1/sqrt(4/2.5))/(prec_95 - 1/sqrt(4/2.5))
    
    ### Proportion of statistically significant effect sizes
    prop_sig <- mean(pnorm(dat$yi/sqrt(dat$vi), lower.tail = FALSE) < alpha)
    
    ### Order data based on vi from large to small
    dat <- dat[order(dat$posi, decreasing = TRUE), ]
    
    ### Random-effects meta-analysis
    res <- metafor::rma(yi = yi, vi = vi, method = method_tau2, data = dat)
    
    ### Cumulative meta-analysis starting with the most precise effect size
    cum_dat <- metafor::cumul(res, order(dat$posi, decreasing = TRUE))
    
    ### Store estimates of cumulative meta-analysis start with a meta-analysis based 
    # on the most precise effect size and end with a meta-analysis based on all
    # included effect sizes
    dat$est_cum <- exp(cum_dat$estimate)
    dat$lb_cum <- exp(cum_dat$ci.lb)
    dat$ub_cum <- exp(cum_dat$ci.ub)
    
    ### Select the rows of dat such that there are no equal values of posi in 
    # the data. Taking into account that the rows that are kept are the ones where the 
    # most effect sizes are included in the cumulative meta-analysis.
    dat_uniq <- dat[rev(duplicated(rev(dat$posi)) == FALSE), ]
    
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
    
    ### Group results of studies together if study's relative position is larger 
    # than 1 (95% power for detecting a small effect)
    if (any(dat_uniq$posi > 1) == TRUE) 
    {
      ind <- which(dat_uniq$posi > 1)
      dat_uniq <- dat_uniq[tail(ind, n = 1):nrow(dat_uniq), ]
      dat_uniq[ind,"posi"] <- 1
    }
    
    if (nr_lines == "all")
    {
      draw_plot_or(dat = dat_uniq, ylim = ylim, alpha = alpha, pub_bias = pub_bias, 
                   prop_sig = prop_sig, pos_sm = pos_sm, pos_me = pos_me, 
                   pos_la = pos_la, main = main, cex.pch = cex.pch)
    } else if (nr_lines == "summary")
    {
      ##### Plot summary results #####
      ### Indexes of cumulative meta-analyses with all effect sizes, the most precise 
      # effect size and the largest number of effect sizes in the regions of the plot 
      # created by the vertical lines
      ind <- c(tail(which(dat_uniq$posi > pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_me & dat_uniq$posi <= pos_sm), n = 1), 
               tail(which(dat_uniq$posi > pos_la & dat_uniq$posi <= pos_me), n = 1), 
               tail(which(dat_uniq$posi <= pos_sm), n = 1), nrow(dat_uniq))
      
      ind <- unique(ind) # Select only unique values in the vector
      
      draw_plot_or(dat = dat_uniq[ind, ], ylim = ylim, alpha = alpha, 
                   pub_bias = pub_bias, prop_sig = prop_sig, pos_sm = pos_sm, 
                   pos_me = pos_me, pos_la = pos_la, main = main, cex.pch = cex.pch)
    }
    
  }
  
  ### Text presenting percentage of studies particular statistical power
  txt1 <- round(sum(dat$posi <= pos_la)/nrow(dat)*100) # Less than 80% for detecting large effect
  txt2 <- round(sum(dat$posi > pos_la)/nrow(dat)*100) # More than 80% for detecting large effect
  txt3 <- round(sum(dat$posi >= pos_me)/nrow(dat)*100) # More than 80% for detecting medium effect
  txt4 <- round(sum(dat$posi >= pos_sm)/nrow(dat)*100) # More than 80% for detecting small effect
  
  text(x = pos_la/2, y = par("usr")[4],
       labels = paste(txt1, "%", sep = ""), cex = par()$cex.lab*1.2, xpd = TRUE)
  text(x = (pos_me+pos_la)/2, y = par("usr")[4],
       labels = paste(txt2, "%", sep = ""), cex = par()$cex.lab*1.2, xpd = TRUE)
  text(x = (pos_sm+pos_me)/2, y = par("usr")[4],
       labels = paste(txt3, "%", sep = ""), cex = par()$cex.lab*1.2, xpd = TRUE)
  text(x = (pos_95+pos_sm)/2, y = par("usr")[4],
       labels = paste(txt4, "%", sep = ""), cex = par()$cex.lab*1.2, xpd = TRUE)
  
  ### Arrows drawn right above the percentages
  arrows(x0 = pos_la/2+strwidth(paste(txt1, "%", sep = ""),
                                cex = par()$cex.lab*1.2)/2, 
         x1 = pos_la/2-strwidth(paste(txt1, "%", sep = ""),
                                cex = par()$cex.lab*1.2)/2, 
         y0 = par("usr")[4] + strheight(paste(txt1, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         y1 = par("usr")[4] + strheight(paste(txt1, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         xpd = TRUE, length = 0.5*strheight(paste("%", sep = ""), units = "inches"))
  
  arrows(x0 = (pos_me+pos_la)/2-strwidth(paste(txt2, "%", sep = ""),
                                         cex = par()$cex.lab*1.2)/2, 
         x1 = (pos_me+pos_la)/2+strwidth(paste(txt2, "%", sep = ""),
                                         cex = par()$cex.lab*1.2)/2, 
         y0 = par("usr")[4] + strheight(paste(txt2, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         y1 = par("usr")[4] + strheight(paste(txt2, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         xpd = TRUE, length = 0.5*strheight(paste("%", sep = ""), units = "inches"))
  
  arrows(x0 = (pos_sm+pos_me)/2-strwidth(paste(txt3, "%", sep = ""),
                                         cex = par()$cex.lab*1.2)/2, 
         x1 = (pos_sm+pos_me)/2+strwidth(paste(txt3, "%", sep = ""),
                                         cex = par()$cex.lab*1.2)/2, 
         y0 = par("usr")[4] + strheight(paste(txt3, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         y1 = par("usr")[4] + strheight(paste(txt3, "%", sep = ""),
                                        cex = par()$cex.lab*1.2), 
         xpd = TRUE, length = 0.5*strheight(paste("%", sep = ""), units = "inches"))
  
  arrows(x0 = (pos_95+pos_sm)/2-strwidth(paste(txt4, "%", sep = ""),
                                         cex = par()$cex.lab*1.2)/2, 
         x1 = (pos_95+pos_sm)/2+strwidth(paste(txt4, "%", sep = ""), 
                                         cex = par()$cex.lab*1.2)/2, 
         y0 = par("usr")[4] + strheight(paste(txt4, "%", sep = ""), 
                                        cex = par()$cex.lab*1.2), 
         y1 = par("usr")[4] + strheight(paste(txt4, "%", sep = ""), 
                                        cex = par()$cex.lab*1.2), 
         xpd = TRUE, length = 0.5*strheight(paste("%", sep = ""), units = "inches"))
  
  ### Text presenting percentage of statistically significant effect sizes
  text(x = 0.5, 
       y = par("usr")[4]-2*strheight(paste(txt1, "%", sep = ""), cex = par()$cex.lab*1.2),
       labels = paste("Sig. = ", round(prop_sig*100,1), "%", sep = ""), 
       cex = par()$cex.lab*1.2)
  
  invisible(dat_uniq)
  
}