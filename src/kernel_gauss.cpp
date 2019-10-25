#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

//' Gaussian kernel
//'
//' @param ksize Side length of the squared kernel.
//' @param sigma Standard deviation of Gaussian weighting function.
//' @return A squared matrix sum to 1.
//' @export
//'
// [[Rcpp::export]]
arma::mat kernel_gauss(int ksize, double sigma){
  int rad = floor(0.5 * ksize);
  vec weight(ksize);
  for(int i = 0; i < ksize; ++i){
    weight(i) = exp((i - rad) * (i - rad) / (-2 * sigma * sigma));
  }
  weight = weight / accu(weight);
  arma::mat kernel = weight * weight.t();
  return kernel;
}
