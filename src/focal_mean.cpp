#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::plugins(openmp)]]

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


//' Focal mean with a given kernel
//'
//' @param x A matrix.
//' @param kernel A kernel which can be produced by [kernel_gauss]
//' @param globe Are data at globe scale? If `TRUE`, both vertical
//' borders will be padded.
//' @return A matrix of focal mean.
//' @export
// [[Rcpp::export]]
arma::mat focal_mean(arma::mat x, arma::mat kernel, bool globe){
  arma::mat out = arma::mat(size(x)).fill(datum::nan);
  int rad = floor(0.5*kernel.n_cols);
  // pad cols
  arma::mat headcols = x.head_cols(rad);
  x.insert_cols(0, x.tail_cols(rad));
  x.insert_cols(x.n_cols, headcols);
  if (!globe) {
    x.head_cols(rad).fill(datum::nan);
    x.tail_cols(rad).fill(datum::nan);
  }
  // pad rows
  arma::mat nanrows = arma::mat(rad, x.n_cols).fill(datum::nan);
  x.insert_rows(0, nanrows);
  x.insert_rows(x.n_rows, nanrows);
  uvec locs = find_finite(x);
  locs = locs.elem(find(locs >= rad * x.n_rows && locs <= x.n_elem - 1 - rad * x.n_rows));
  int nlocs = locs.n_elem;
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int loc = 0; loc < nlocs; ++loc){
    // declare inside as private variables for parallel
    // "%" here is modulo operator, not elewise multiplication
    int col = locs(loc) / x.n_rows;
    int row = locs(loc) % x.n_rows;
    arma::mat xsubmat = x.submat(row-rad, col-rad, row+rad, col+rad);
    uvec inx = find_finite(xsubmat);
    arma::vec xsub = vectorise(xsubmat.elem(inx));
    arma::vec w = kernel.elem(inx);
    w = w / accu(w);
    out(row-rad, col-rad) = accu(w % xsub);
  }
  return out;
}
