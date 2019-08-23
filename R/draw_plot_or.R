### Function for creating plot when odds ratio is the effect size measure
draw_plot_or <- function(dat, ylim, alpha = alpha, pub_bias = pub_bias, prop_sig, 
                         pos_sm, pos_me, pos_la, main, cex.pch = cex.pch) 
{
  
  ### Create plot
  with(dat, plot(x = posi, y = est_cum, type = "p", pch = 16, 
                 ylim = ylim, xlim = c(0, 1.05), xaxt = "n", yaxt = "n", 
                 bty = "n", xlab = "", cex = cex.pch, cex.lab = par()$cex.lab,
                 ylab = ""))
  
  ### Add title to plot
  title(main, line = 2.5)
  
  ### Add label x-axis
  mtext("Precision", side = 1, cex = par()$cex.lab, line = 3.8)
  
  ### Add label y-axis
  mtext("Effect size (OR)", side = 2, cex = par()$cex.lab, line = par()$mgp[1])     
  
  ### Draw confidence intervals
  with(dat[nrow(dat), ], arrows(x0 = posi, y0 = lb_cum, y1 = ub_cum, 
                                code = 3, angle = 90, length = 0.1)) # First CI is black
  
  with(dat[1:(nrow(dat)-1), ], arrows(x0 = posi, y0 = lb_cum, y1 = ub_cum, # Other CIs are gray 
                                      code = 3, angle = 90, length = 0.1, col = "gray"))
  
  ### Create x-axis
  axis(1, at = seq(0, 1, 0.1), cex = par()$cex.axis)
  
  ### Create y-axis
  axis(2, at = round(seq(ylim[1], ylim[2], length.out = 8), 2), las = 1)
  
  if (prop_sig > 0.8 & pub_bias == TRUE)
  { ### Add points for cumulative meta-analysis based on Mill's ratios if
    # proportion statistically significant effect sizes is larger than 0.8
    with(dat, points(x = posi, y = pub_est, cex = cex.pch, pch = 8))  
  }
  
  ### Vertical line reflecting required sample size for particular statistical power
  abline(v = pos_la, col = "gray", lty = 2) # 80% for large effect
  abline(v = pos_me, col = "gray", lty = 2) # 80% for medium effect
  abline(v = pos_sm, col = "gray", lty = 2) # 80% for small effect
  
  ### Horizontal line reflecting estimate of meta-analysis
  abline(h = dat$est_cum[nrow(dat)], lty = 3)
  
  ### Horizontal line at no effect
  abline(h = 1)
  
  ### Letters indicating to what effect size each vertical line belongs
  mtext("L", side = 3, at = pos_la, cex = par()$cex.lab)
  mtext("M", side = 3, at = pos_me, cex = par()$cex.lab)
  mtext("S", side = 3, at = pos_sm, cex = par()$cex.lab)
  
}