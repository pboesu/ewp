// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dewp3_cpp
double dewp3_cpp(int x, double lambda, double beta1, double beta2);
RcppExport SEXP _ewp_dewp3_cpp(SEXP xSEXP, SEXP lambdaSEXP, SEXP beta1SEXP, SEXP beta2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    rcpp_result_gen = Rcpp::wrap(dewp3_cpp(x, lambda, beta1, beta2));
    return rcpp_result_gen;
END_RCPP
}
// pllik3_part_cpp
double pllik3_part_cpp(Rcpp::NumericVector X, Rcpp::NumericVector lambda, double beta1, double beta2);
RcppExport SEXP _ewp_pllik3_part_cpp(SEXP XSEXP, SEXP lambdaSEXP, SEXP beta1SEXP, SEXP beta2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    rcpp_result_gen = Rcpp::wrap(pllik3_part_cpp(X, lambda, beta1, beta2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ewp_dewp3_cpp", (DL_FUNC) &_ewp_dewp3_cpp, 4},
    {"_ewp_pllik3_part_cpp", (DL_FUNC) &_ewp_pllik3_part_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ewp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
