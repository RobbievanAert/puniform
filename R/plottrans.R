### Function to generate plot of transformed p-values
plottrans <- function(tr.q, ksig) {
  plot(x = seq(0,1,1/(ksig+1))[2:(ksig+1)], y = sort(tr.q), ylim = c(0,1), 
       xlim = c(0,1), xlab = expression(paste("Expected conditional ", italic(p), "-values")), 
                                              ylab = expression(paste("Observed conditional ", italic(p), "-values")),
       pch = 16)
  lines(seq(0,1,1/(ksig+1)), seq(0,1,1/(ksig+1)), lty = 2, lwd = 1.5)
}