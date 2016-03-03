#' hybrid
#'
#' @export

hybrid <- function(m1i, m2i, mi, ri, sd1i, sd2i, sdi, n1i, n2i, ni, tobs, alpha, side) {

  if (!missing("mi") & !missing("ni") & !missing("sdi")) {
    measure <- "M"
    es <- escompute(mi = mi, ni = ni, sdi = sdi, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, mi = mi, sdi = sdi, ni = ni, measure = measure, side = side)
  } else if (!missing("ni") & !missing("tobs")) {
    measure <- "MT"
    es <- escompute(ni = ni, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
  } else if (!missing("m1i") & !missing("m2i") & !missing("n1i") & !missing("n2i") & !missing("sd1i") & !missing("sd2i")) {
    measure <- "MD"
    es <- escompute(m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, m1i = m1i, m2i = m2i, n1i = n1i, n2i = n2i, sd1i = sd1i, sd2i = sd2i, measure = measure, side = side)
  } else if (!missing("n1i") & !missing("n2i") & !missing("tobs")) {
    measure <- "MDT"
    es <- escompute(n1i = n1i, n2i = n2i, tobs = tobs, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, tobs = tobs, measure = measure, side = side)
  } else if (!missing("ri") & !missing("ni")) {
    measure <- "COR"
    es <- escompute(ri = ri, ni = ni, alpha = alpha/2, side = side, measure = measure)
    res.repl <- repl(es = es, measure = measure, side = side)
  }

  if (es$pval[1] > alpha/2) { stop("Original study is not statistically significant") }

  res1 <- hy(es = es, measure = measure, side = side, alpha = alpha/2)

  if (es$pval[1] < alpha/4) {
    res2 <- data.frame(est.hyr = res1$est.hy, ci.lb.hyr = res1$ci.lb.hy, ci.ub.hyr = res1$ci.ub.hy, stat.hyr = res1$x.hy, pval.hyr = res1$pval.hy, pval.o = res.repl$pval.o)
  } else { res2 <- data.frame(est.hyr = res.repl$est.repl, ci.lb.hyr = res.repl$ci.lb.repl, ci.ub.hyr = res.repl$ci.ub.repl, stat.hyr = res.repl$stat.repl,
                              pval.hyr = res.repl$pval.repl, pval.o = res.repl$pval.o)
  }

  res3 <- hy0(es = es, res1 = res1, alpha = alpha/2)

  res4 <- fix(yi = es$yi, vi = es$vi, measure = measure, side = side)

  w <- list(est.hy = res1$est, ci.lb.hy = res1$ci.lb, ci.ub.hy = res1$ci.ub, x.hy = res1$x, pval.hy = res1$pval, measure = measure, est.hyr = res2$est.hyr,
            ci.lb.hyr = res2$ci.lb.hyr, ci.ub.hyr = res2$ci.ub.hyr, stat.hyr = res2$stat.hyr, pval.hyr = res2$pval.hyr, pval.o = res2$pval.o, est.hy0 = res3$est.hy0,
            ci.lb.hy0 = res3$ci.lb.hy0, ci.ub.hy0 = res3$ci.ub.hy0, x.hy0 = res3$x.hy0, pval.hy0 = res3$pval.hy0, est.fe = res4$est.fe, se.fe = res4$se.fe,
            ci.lb.fe = res4$ci.lb.fe, ci.ub.fe = res4$ci.ub.fe, zval.fe = res4$zval.fe, pval.fe = res4$pval.fe, est.repl = res.repl$est.repl, se.repl = res.repl$se.repl,
            ci.lb.repl = res.repl$ci.lb.repl, ci.ub.repl = res.repl$ci.ub.repl, stat.repl = res.repl$stat.repl, pval.repl = res.repl$pval.repl)
  class(w) <- "hybridoutput"
  return(w)
}
