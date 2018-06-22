// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ml_est
double ml_est(double est, double tau, NumericVector yi, NumericVector vi, NumericVector ycv);
RcppExport SEXP _puniform_ml_est(SEXP estSEXP, SEXP tauSEXP, SEXP yiSEXP, SEXP viSEXP, SEXP ycvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type est(estSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vi(viSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ycv(ycvSEXP);
    rcpp_result_gen = Rcpp::wrap(ml_est(est, tau, yi, vi, ycv));
    return rcpp_result_gen;
END_RCPP
}
// ml_tau
double ml_tau(double tau, double d, NumericVector yi, NumericVector vi, NumericVector ycv);
RcppExport SEXP _puniform_ml_tau(SEXP tauSEXP, SEXP dSEXP, SEXP yiSEXP, SEXP viSEXP, SEXP ycvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vi(viSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ycv(ycvSEXP);
    rcpp_result_gen = Rcpp::wrap(ml_tau(tau, d, yi, vi, ycv));
    return rcpp_result_gen;
END_RCPP
}
// approx_C
double approx_C(double yi, double tot_var, double ycv, double est);
RcppExport SEXP _puniform_approx_C(SEXP yiSEXP, SEXP tot_varSEXP, SEXP ycvSEXP, SEXP estSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< double >::type tot_var(tot_varSEXP);
    Rcpp::traits::input_parameter< double >::type ycv(ycvSEXP);
    Rcpp::traits::input_parameter< double >::type est(estSEXP);
    rcpp_result_gen = Rcpp::wrap(approx_C(yi, tot_var, ycv, est));
    return rcpp_result_gen;
END_RCPP
}
// pdist_nsig
double pdist_nsig(double est, double tau, NumericVector yi, NumericVector vi, String param, NumericVector ycv, String method, String val, double cv_P);
RcppExport SEXP _puniform_pdist_nsig(SEXP estSEXP, SEXP tauSEXP, SEXP yiSEXP, SEXP viSEXP, SEXP paramSEXP, SEXP ycvSEXP, SEXP methodSEXP, SEXP valSEXP, SEXP cv_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type est(estSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vi(viSEXP);
    Rcpp::traits::input_parameter< String >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ycv(ycvSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< String >::type val(valSEXP);
    Rcpp::traits::input_parameter< double >::type cv_P(cv_PSEXP);
    rcpp_result_gen = Rcpp::wrap(pdist_nsig(est, tau, yi, vi, param, ycv, method, val, cv_P));
    return rcpp_result_gen;
END_RCPP
}
// trq
NumericVector trq(double est, double tau, NumericVector yi, NumericVector vi, NumericVector ycv, String param);
RcppExport SEXP _puniform_trq(SEXP estSEXP, SEXP tauSEXP, SEXP yiSEXP, SEXP viSEXP, SEXP ycvSEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type est(estSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vi(viSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ycv(ycvSEXP);
    Rcpp::traits::input_parameter< String >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(trq(est, tau, yi, vi, ycv, param));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_puniform_ml_est", (DL_FUNC) &_puniform_ml_est, 5},
    {"_puniform_ml_tau", (DL_FUNC) &_puniform_ml_tau, 5},
    {"_puniform_approx_C", (DL_FUNC) &_puniform_approx_C, 4},
    {"_puniform_pdist_nsig", (DL_FUNC) &_puniform_pdist_nsig, 9},
    {"_puniform_trq", (DL_FUNC) &_puniform_trq, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_puniform(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
