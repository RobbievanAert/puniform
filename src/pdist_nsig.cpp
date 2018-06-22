#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double approx_C(double yi, double tot_var, double ycv, double est)
{
  
  double zval;
  double zest;
  double zcv;
  double a;
  double oz;
  double oz_a;
  double out;
  
  zval = yi/sqrt(tot_var);
  zest = est/sqrt(tot_var);
  zcv = ycv/sqrt(tot_var);
  a = zval-zcv;
  
  oz = 1-1/(pow(zcv-zest,2)+2)+1/((pow(zcv-zest,2)+2)*(pow(zcv-zest,2)+4))-
    5/((pow(zcv-zest,2)+2)*(pow(zcv-zest,2)+4)*(pow(zcv-zest,2)+6));
  oz_a = 1-1/(pow(zval-zest,2)+2)+1/((pow(zval-zest,2)+2)*(pow(zval-zest,2)+4))-
    5/((pow(zval-zest,2)+2)*(pow(zval-zest,2)+4)*(pow(zval-zest,2)+6));
  out = exp(-(zcv-zest)*a-0.5*pow(a,2))*(zcv-zest)/(zval-zest)*oz_a/oz;
  
  return out;
}

// [[Rcpp::export]]
double pdist_nsig(double est, double tau, NumericVector yi, NumericVector vi, String param,
                NumericVector ycv, String method, String val, double cv_P)
{
  int n = yi.size();
  NumericVector q(n);
  NumericVector tot_var(n);
  double sum_q = 0;
  double out = 0;
  
  for(int i = 0; i < n; ++i){
    tot_var[i] = vi[i] + pow(tau, 2);
    
    if (yi[i] > ycv[i]) 
    {
      if (param == "est") 
        {
        if (ycv[i]/sqrt(tot_var[i]) - est/sqrt(tot_var[i]) <= 38)
        {
          q[i] = Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0)/Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0);
        } else {
          q[i] = approx_C(yi[i], tot_var[i], ycv[i], est);
        }
      } else if (param == "tau")
      {
        if (est < ycv[i] && est < yi[i])
        {
          q[i] = (Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0)-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0);
        } else if (est > ycv[i] && est < yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/0.5;
        } else if (est > ycv[i] && est > yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/
            (0.5-Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0));
        }
      }
      
    } else if (yi[i] < ycv[i]) {
      
      if (param == "est") 
      {
        if (yi[i]/sqrt(tot_var[i]) - est/sqrt(tot_var[i]) >= -30)
        {
          q[i] = 1-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0)/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0);
        } else {
          q[i] = 1-approx_C(-yi[i], tot_var[i], -ycv[i], -est);
        }
      } else if (param == "tau") 
      { 
        if (est < ycv[i] && est < yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/
            (0.5-Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0));
        } else if (est < ycv[i] && est > yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/0.5;
        } else if (est > ycv[i] && est > yi[i])
        {
          q[i] = (Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0)-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0);
        }
      }
      
    }
    
    if (method == "P")
    {
      sum_q += q[i];
    } else if (method == "LNP")
    {
      sum_q += -log(q)[i];
    }
    
  }
  
  if (val == "es")
  {
    if (method == "P")
    {
      out = sum_q-0.5*n;
    } else if (method == "LNP")
    {
      out = sum_q-n;
    }
  } else if (val == "ci.lb")
  {
    if (method == "P")
    {
      if (param == "est")
      {
        out = sum_q-cv_P;
      } else if (param == "tau")
      {
        out = sum_q-(n-cv_P);
      }
    } else if (method == "LNP")
    {
      if (param == "est")
      {
        out = sum_q-R::qgamma(0.975, n, 1, 1, 0);
      } else if (param == "tau")
      {
        out = sum_q-R::qgamma(0.025, n, 1, 1, 0);
      }
    }
  } else if (val == "ci.ub")
  {
    if (method == "P")
    {
      if (param == "est")
      {
        out = sum_q-(n-cv_P);
      } else if (param == "tau")
      {
        out = sum_q-cv_P;
      }
    } else if (method == "LNP")
    {
      if (param == "est")
      {
        out = sum_q-R::qgamma(0.025, n, 1, 1, 0);
      } else if (param == "tau")
      {
        out = sum_q-R::qgamma(0.975, n, 1, 1, 0);
      }
    }
  }
  
  return out;
  
}

// [[Rcpp::export]]
NumericVector trq(double est, double tau, NumericVector yi, NumericVector vi, 
                  NumericVector ycv, String param)
{
  int n = yi.size();
  NumericVector q(n);
  NumericVector tot_var(n);
  
  for(int i = 0; i < n; ++i){
    tot_var[i] = vi[i] + pow(tau, 2);
    
    if (yi[i] > ycv[i]) 
    {
      if (param == "est") 
        {
        if (ycv[i]/sqrt(tot_var[i]) - est/sqrt(tot_var[i]) <= 38)
        {
          q[i] = Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0)/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0);
        } else {
          q[i] = approx_C(yi[i], tot_var[i], ycv[i], est);
        }
      } else if (param == "tau")
      {
        if (est < ycv[i] && est < yi[i])
        {
          q[i] = (Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0)-
            Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0);
        } else if (est > ycv[i] && est < yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/0.5;
        } else if (est > ycv[i] && est > yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/
            (0.5-Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0));
        }
      }
      
    } else if (yi[i] < ycv[i]) 
      {
      if (param == "est") 
      {
        if (yi[i]/sqrt(tot_var[i]) - est/sqrt(tot_var[i]) >= -30)
        {
          q[i] = 1-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0)/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0);
        } else {
          q[i] = 1-approx_C(-yi[i], tot_var[i], -ycv[i], -est);
        }
      } else if (param == "tau") 
      { 
        if (est < ycv[i] && est < yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 0, 0))/
            (0.5-Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 0, 0));
        } else if (est < ycv[i] && est > yi[i])
        {
          q[i] = (0.5-Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/0.5;
        } else if (est > ycv[i] && est > yi[i])
        {
          q[i] = (Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0)-
            Rf_pnorm5(yi[i], est, sqrt(tot_var[i]), 1, 0))/
            Rf_pnorm5(ycv[i], est, sqrt(tot_var[i]), 1, 0);
        }
      }
      
    }
  }
  
  return q;
}