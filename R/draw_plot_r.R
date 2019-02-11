### Function for creating plot when correlation coefficient is the effect 
# size measure
draw_plot_r <- function(dat, ylim, alpha = alpha, prop_sig, main, cex.pch = cex.pch) 
{
  
  ### Create plot
  with(dat, plot(x = sqrt(ni), y = est_cum, type = "p", pch = 16, ylim = ylim,
                 xlim = c(0, sqrt(1300)), xaxt = "n", yaxt = "n", bty = "n",
                 xlab = "", cex = cex.pch, cex.lab = par()$cex.lab,
                 ylab = ""))
  
  ### Add title to plot
  title(main, line = 2.5)
  
  ### Add label x-axis
  mtext(expression(italic(N)), side = 1, cex = par()$cex.lab, line = par()$mgp[1])
  
  ### Add label y-axis
  mtext(expression(paste("Effect size (", italic(r), ")"), sep = ""), side = 2, 
        cex = par()$cex.lab, line = par()$mgp[1])     
  
  ### Draw confidence intervals
  with(dat[nrow(dat), ], arrows(x0 = sqrt(ni), y0 = lb_cum, y1 = ub_cum, code = 3, 
                                angle = 90, length = 0.1)) # First CI is black
  
  with(dat[1:(nrow(dat)-1), ], arrows(x0 = sqrt(ni), y0 = lb_cum, y1 = ub_cum, # Other CIs are gray 
                                      code = 3, angle = 90, length = 0.1, col = "gray"))
  
  ### Create x-axis
  axis(1, at = sqrt(c(0, 25, 50, 100, 200, 400, 800, 1300)), labels = rep("", 8))
  text(x = sqrt(c(0, 25, 50, 100, 200, 400, 800, 1300)), y = ylim[1]-sum(abs(ylim))/10, 
       srt = 45, adj = 1, xpd = TRUE, labels = c(0, 25, 50, 100, 200, 400, 800, 1300), 
       cex = par()$cex.axis)
  
  ### Create y-axis
  axis(2, at = round(seq(ylim[1], ylim[2], length.out = 8), 2), las = 1)
  
  if (prop_sig > 0.8)
  { ### Add points for cumulative meta-analysis based on Mill's ratios if
    # proportion statistically significant effect sizes is larger than 0.8
    with(dat, points(x = sqrt(ni), y = pub_est, cex = cex.pch, pch = 8))  
  }
  
  ### Vertical line reflecting required sample size for particular statistical power
  abline(v = sqrt(29), col = "gray", lty = 2) # 80% for large effect
  abline(v = sqrt(84), col = "gray", lty = 2) # 80% for medium effect
  abline(v = sqrt(782), col = "gray", lty = 2) # 80% for small effect
  
  ### Horizontal line reflecting estimate of meta-analysis
  abline(h = dat$est_cum[nrow(dat)], lty = 3)
  
  ### Horizontal line at no effect
  abline(h = 0)
  
  ### Letters indicating to what effect size each vertical line belongs
  mtext("L", side = 3, at = sqrt(29), cex = par()$cex.lab)
  mtext("M", side = 3, at = sqrt(84), cex = par()$cex.lab)
  mtext("S", side = 3, at = sqrt(782), cex = par()$cex.lab)
  
}
