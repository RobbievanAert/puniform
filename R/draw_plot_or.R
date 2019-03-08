### Function for creating plot when odds ratio is the effect size measure
draw_plot_or <- function(dat, ylim, alpha = alpha, pub_bias = pub_bias, prop_sig, 
                         main, cex.pch = cex.pch) 
{
  
  ### Create plot
  with(dat, plot(x = preci, y = est_cum, type = "p", pch = 16, 
                 ylim = ylim, xlim = c(0.17, 7.2), xaxt = "n", yaxt = "n", 
                 bty = "n", xlab = "", cex = cex.pch, cex.lab = par()$cex.lab,
                 ylab = ""))
  
  ### Add title to plot
  title(main, line = 2.5)
  
  ### Add label x-axis
  mtext(expression(italic(N)), side = 1, cex = par()$cex.lab, line = 3.8)
  
  ### Add label y-axis
  mtext("Effect size (OR)", side = 2, cex = par()$cex.lab, line = par()$mgp[1])     
  
  ### Draw confidence intervals
  with(dat[nrow(dat), ], arrows(x0 = preci, y0 = lb_cum, y1 = ub_cum, 
                                code = 3, angle = 90, length = 0.1)) # First CI is black
  
  with(dat[1:(nrow(dat)-1), ], arrows(x0 = preci, y0 = lb_cum, y1 = ub_cum, # Other CIs are gray 
                                      code = 3, angle = 90, length = 0.1, col = "gray"))
  
  ### Create x-axis
  labels <- c(0, 20, 40, 80, 160, 320, 640, 800)
  at <- 1/sqrt(1/(labels/4)+1/(labels/4)+1/(labels/4)+1/(labels/4))
  axis(1, at = at, labels = rep("", 8))
  text(x = at, y = ylim[1]-sum(abs(ylim))/10, srt = 45, adj = 1, xpd = TRUE, 
       labels = labels, cex = par()$cex.axis)
  
  ### Create y-axis
  axis(2, at = round(seq(ylim[1], ylim[2], length.out = 8), 2), las = 1)
  
  if (prop_sig > 0.8 & pub_bias == TRUE)
  { ### Add points for cumulative meta-analysis based on Mill's ratios if
    # proportion statistically significant effect sizes is larger than 0.8
    with(dat, points(x = preci, y = pub_est, cex = cex.pch, pch = 8))  
  }
  
  ### Vertical line reflecting required sample size for particular statistical power
  abline(v = 1.472, col = "gray", lty = 2) # 80% for large effect
  abline(v = 2.252, col = "gray", lty = 2) # 80% for medium effect
  abline(v = 5.4, col = "gray", lty = 2) # 80% for small effect
  
  ### Horizontal line reflecting estimate of meta-analysis
  abline(h = dat$est_cum[nrow(dat)], lty = 3)
  
  ### Horizontal line at no effect
  abline(h = 1)
  
  ### Letters indicating to what effect size each vertical line belongs
  mtext("L", side = 3, at = 1.472, cex = par()$cex.lab)
  mtext("M", side = 3, at = 2.252, cex = par()$cex.lab)
  mtext("S", side = 3, at = 5.4, cex = par()$cex.lab)
  
}