#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::plugins(openmp)]]

//' Focal statistics on spatial window
//'
//' @param x A matrix.
//' @param y The other matrix.
//' @param rescale Rescale or just clamp before computation.
//' @param xmin,xmax,ymin,ymax Rescale or clamp parameters.
//' If `NA_real_`, the default, will be set as extrema of data.
//' @param ksize Side length of the spatial window.
//' @param globe Are data at the global scale? If `TRUE`, both vertical
//' borders will be padded.
//' @param stat Statistic to return. One of 'xmn', 'ymn', 'xsd', or 'ysd'.
//' @return A matirx.
//' @export
// [[Rcpp::export]]
arma::mat focal_stat_sw(arma::mat x, arma::mat y, bool rescale = false,
                        double xmin = NA_REAL, double xmax = NA_REAL,
                        double ymin = NA_REAL, double ymax = NA_REAL,
                        double ksize = 9, bool globe = false,
                        std::string stat = "xmn"){
  if (x.has_nan() || y.has_nan()){
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
    if (find_finite(x).is_empty()) Rcpp::stop("x and y have no valid overlap!");
  }
  if (!std::isfinite(xmin)) xmin = x.min();
  if (!std::isfinite(xmax)) xmax = x.max();
  if (!std::isfinite(ymin)) ymin = y.min();
  if (!std::isfinite(ymax)) ymax = y.max();
  if (xmin > xmax) Rcpp::stop("xmin > xmax, please reset them!");
  if (ymin > ymax) Rcpp::stop("ymin > ymax, please reset them!");
  if (xmax < x.min() || xmin > x.max()) Rcpp::stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < y.min() || ymin > y.max()) Rcpp::stop("[ymin, ymax] is beyond the range of x!");

  x = clamp(x, xmin, xmax);
  y = clamp(y, ymin, ymax);
  if (rescale) {
    // avoid zero-denominator
    if (xmin == xmax) x.elem(find_finite(x)).ones(); else x = (x-xmin) / (xmax-xmin);
    if (ymin == ymax) y.elem(find_finite(y)).ones(); else y = (y-ymin) / (ymax-ymin);
  }

  arma::mat xmn = arma::mat(size(x)).fill(datum::nan);
  arma::mat ymn = arma::mat(size(x)).fill(datum::nan);
  arma::mat xsd = arma::mat(size(x)).fill(datum::nan);
  arma::mat ysd = arma::mat(size(x)).fill(datum::nan);
  int rad = floor(0.5 * ksize);

  // pad cols
  arma::mat xheadcols = x.head_cols(rad);
  x.insert_cols(0, x.tail_cols(rad));
  x.insert_cols(x.n_cols, xheadcols);
  arma::mat yheadcols = y.head_cols(rad);
  y.insert_cols(0, y.tail_cols(rad));
  y.insert_cols(y.n_cols, yheadcols);
  if (!globe) {
    x.head_cols(rad).fill(datum::nan);
    x.tail_cols(rad).fill(datum::nan);
    y.head_cols(rad).fill(datum::nan);
    y.tail_cols(rad).fill(datum::nan);
  }
  // pad rows
  arma::mat nanrows = arma::mat(rad, x.n_cols).fill(datum::nan);
  x.insert_rows(0, nanrows);
  x.insert_rows(x.n_rows, nanrows);
  y.insert_rows(0, nanrows);
  y.insert_rows(y.n_rows, nanrows);
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
    arma::mat ysubmat = y.submat(row-rad, col-rad, row+rad, col+rad);
    uvec inx = find_finite(xsubmat);
    if (inx.is_empty()) continue;
    arma::vec xsub = vectorise(xsubmat.elem(inx));
    arma::vec ysub = vectorise(ysubmat.elem(inx));
    xmn(row-rad, col-rad) = mean(xsub);
    ymn(row-rad, col-rad) = mean(ysub);
    xsd(row-rad, col-rad) = stddev(xsub);
    ysd(row-rad, col-rad) = stddev(ysub);
  }
  if (stat == "xmn") return xmn;
  if (stat == "ymn") return ymn;
  if (stat == "xsd") return xsd;
  if (stat == "ysd") return ysd;
  Rcpp::stop("stat should be 'xmn' or 'ymn', 'xsd', 'ysd'!");
}


//' Focal statistics on temporal window
//'
//' This function computes the GCSM at each cell treating the 3rd dimension
//' as a temporal window.
//' @inheritParams focal_stat_sw
//' @param xmin,xmax,ymin,ymax Parameters that determine fuzzy sets.
//' If `NA_real_`, the default, will be set as extrema of data.
//' @param xxx A 3d array. The 3rd one is time.
//' @param yyy The other array.
//' @return A matrix.
//' @export
// [[Rcpp::export]]
arma::mat focal_stat_tw(arma::cube xxx, arma::cube yyy, bool rescale = false,
                        double xmin = NA_REAL, double xmax = NA_REAL,
                        double ymin = NA_REAL, double ymax = NA_REAL,
                        std::string stat = "xmn"){
  // make sure x and y are not altered, use const make nonsense since they are SEXP, i.e. pointer
  // see https://github.com/RcppCore/RcppArmadillo/issues/253#issue-425224004
  // and discussion in https://github.com/RcppCore/RcppArmadillo/pull/109#issuecomment-255484817
  cube x = xxx;
  cube y = yyy;
  if (x.has_nan() || y.has_nan()){
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
    if (find_finite(x).is_empty()) Rcpp::stop("x and y have no valid overlap!");
  }
  if (!std::isfinite(xmin)) xmin = x.min();
  if (!std::isfinite(xmax)) xmax = x.max();
  if (!std::isfinite(ymin)) ymin = y.min();
  if (!std::isfinite(ymax)) ymax = y.max();
  if (xmin > xmax) Rcpp::stop("xmin > xmax, please reset them!");
  if (ymin > ymax) Rcpp::stop("ymin > ymax, please reset them!");
  if (xmax < x.min() || xmin > x.max()) Rcpp::stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < y.min() || ymin > y.max()) Rcpp::stop("[ymin, ymax] is beyond the range of x!");
  x = clamp(x, xmin, xmax);
  y = clamp(y, ymin, ymax);
  if (rescale) {
    // avoid zero-denominator
    if (xmin == xmax) x.elem(find_finite(x)).ones(); else x = (x-xmin) / (xmax-xmin);
    if (ymin == ymax) y.elem(find_finite(y)).ones(); else y = (y-ymin) / (ymax-ymin);
  }

  arma::mat xmn = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat ymn = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat xsd = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat ysd = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  // uvec locs = find_finite(x.slice(0));
  // int nlocs = locs.n_elem;
  int nlocs = x.slice(0).n_elem;
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int loc = 0; loc < nlocs; ++loc){
    // declare inside as private variables for parallel
    // "%" here is modulo operator, not elewise multiplication
    int col = loc / x.n_rows;
    int row = loc % x.n_rows;
    arma::vec xtube = x.tube(row, col);
    arma::vec ytube = y.tube(row, col);
    uvec inx = find_finite(xtube);
    if (inx.is_empty()) continue;
    arma::vec xsub = vectorise(xtube.elem(inx));
    arma::vec ysub = vectorise(ytube.elem(inx));
    xmn(row, col) = mean(xsub);
    ymn(row, col) = mean(ysub);
    xsd(row, col) = stddev(xsub);
    ysd(row, col) = stddev(ysub);
  }
  if (stat == "xmn") return xmn;
  if (stat == "ymn") return ymn;
  if (stat == "xsd") return xsd;
  if (stat == "ysd") return ysd;
  Rcpp::stop("stat should be 'xmn' or 'ymn', 'xsd', 'ysd'!");
}

