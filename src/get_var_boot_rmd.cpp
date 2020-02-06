#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec mvrnorm_cpp(arma::mat Sigma, int dv) 
{
  arma::mat x = arma::randn(1, dv);
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
    sams = mvrnorm_cpp(Sigma, dv);
    boot[i] = var(sams);
  }
  
  boot_var = arma::mean(boot);
  
  return boot_var;
}