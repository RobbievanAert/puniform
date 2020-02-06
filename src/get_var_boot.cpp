#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat mvrnorm_cpp(int n, arma::mat Sigma, int dv) 
{
  arma::mat x = arma::randn(n, dv);
  arma::mat y = x * arma::chol(Sigma);
  return y.t();
}

// [[Rcpp::export]]
double get_var_boot_rmd(arma::mat Sigma, int dv, int reps)
{
  arma::vec sams;
  arma::vec boot(reps);
  double boot_var;
  
  for(int i = 0; i < reps; ++i){
    sams = mvrnorm_cpp(1, Sigma, dv);
    boot[i] = var(sams);
  }
  
  boot_var = arma::mean(boot);
  
  return boot_var;
}

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
    sams = mvrnorm_cpp(n, Sigma, dv+1);
    
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