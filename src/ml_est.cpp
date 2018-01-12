#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ml_est(double est, double tau, NumericVector yi, NumericVector vi, NumericVector ycv)
{
  int n = yi.size();
  NumericVector q(n);
  NumericVector tot_var(n);
  double prod_q = 1;
  
  for(int i = 0; i < n; ++i){
    tot_var[i] = vi[i] + pow(tau, 2);
    
    if (yi[i] > ycv[i]) {
      q[i] = Rf_dnorm4(yi[i], est, sqrt(tot_var[i]), 0)/Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0);
    } else {
      q[i] = Rf_dnorm4(yi[i], est, sqrt(tot_var[i]), 0)/Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0);
    }
    
    prod_q *= q[i];
    
  }
  
  return log(prod_q);
}