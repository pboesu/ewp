// We can now use the BH package
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/factorials.hpp>

using namespace Rcpp;

double w_k3(double beta1, double beta2, double k, double lambda){
  if(k <= lambda){
    return exp(-1 * beta1 * (lambda - k));
  } else {
    return exp(-1 * beta2 * (k - lambda));
  }
}//EWP_3 weights


double W_inner3(double beta1, double beta2, double k, double lambda){
    return exp(-1 * lambda) * pow(lambda, k) * w_k3(beta1, beta2, k, lambda)/boost::math::factorial<double>(k);
  }


double W3(double beta1, double beta2, double lambda, int sum_limit = 30){
  double out = 0;
  for (int i = 0; i < sum_limit+1; i++){
    out += W_inner3(beta1, beta2, i,lambda);
  }
  return out;
}


//' Probability mass function of the three-parameter EWP
//'
//' @param x vector of (positive integer) quantiles.
//' @param lambda centrality parameter
//' @param beta1 lower-tail dispersion parameter
//' @param beta2 upper tail dispersion parameter
//'
//' @return a probability mass
//' @export
// [[Rcpp::export]]
double dewp3_cpp(int x, double lambda, double beta1, double beta2){
  return exp(-1 * lambda) * pow(lambda, x) * w_k3(beta1, beta2, x, lambda)/(W3(beta1, beta2, lambda)*boost::math::factorial<double>(x));
}


// [[Rcpp::export]]
double pllik3_part_cpp(Rcpp::NumericVector X, Rcpp::NumericVector lambda, double beta1, double beta2){
  double ll = 0;
  for (int i = 0; i < X.size(); i++){
    ll += log(dewp3_cpp(X[i],lambda[i],beta1,beta2));
  }
  return(-1*ll);
}
