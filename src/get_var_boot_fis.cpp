#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double get_var_boot_fis(arma::mat Sigma, int n, int dv, int reps)
{
  
  arma::vec mu(dv+1);
  arma::mat sams;
  arma::rowvec x1;
  arma::rowvec y;
  arma::vec ri_boot(dv);
  arma::vec yi_boot(dv);
  arma::vec boot(reps);
  double boot_var;
  
  for(int i = 0; i < reps; ++i){
    sams = arma::mvnrnd(mu, Sigma, n);
    
    x1 = sams.row(0);
    
    for(int q = 0; q < dv; ++q){
      
      y = sams.row(q+1);
      ri_boot[q] = arma::cor(x1,y).eval()(0,0);
      yi_boot[q] = 0.5 * log((1+ri_boot[q])/(1-ri_boot[q]));
    }
    
    boot[i] = var(yi_boot);
    
  }
  
  boot_var = arma::mean(boot);
  
  return boot_var;
}