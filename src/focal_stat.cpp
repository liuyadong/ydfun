#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

//' Focal statistics on spatial window
//'
//' @param x A matrix.
//' @param y The other matrix.
//' @param fuzzy Fuzzify or not prior to calculation.
//' @param xmin Lower limit of \code{x}. Defaults to \code{min(x)}.
//' @param xmax Upper limit of \code{x}. Defaults to \code{max(x)}.
//' @param ymin Lower limit of \code{y}. Defaults to \code{min(y)}.
//' @param ymax Upper limit of \code{y}. Defaults to \code{max(y)}.
//' @param ksize Side length of the spatial window.
//' @param global Are data at global scale? If \code{TRUE}, both vertical
//' borders will be padded.
//' @param stat Statistic to return. One of \code{c('xmn', 'ymn', 'xsd', 'ysd')}.
//' @return A matirx of \code{stat}.
//' @export
// [[Rcpp::export]]
arma::mat focal_stat_2d(arma::mat x, arma::mat y, bool fuzzy,
                        double xmin = NA_REAL, double xmax = NA_REAL,
                        double ymin = NA_REAL, double ymax = NA_REAL,
                        double ksize = 9, bool global = false,
                        Rcpp::String stat = "xmn"){
  if (x.has_nan() || y.has_nan()){
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
  }
  if (!is_finite(xmin)) xmin = x.min();
  if (!is_finite(xmax)) xmax = x.max();
  if (!is_finite(ymin)) ymin = y.min();
  if (!is_finite(ymax)) ymax = y.max();
  // x or y may be a constant or NA
  if (xmin < xmax) x = clamp(x, xmin, xmax);
  if (ymin < ymax) y = clamp(y, ymin, ymax);
  if (fuzzy){
    if (xmin < xmax) x = (x-xmin) / (xmax-xmin); else x.elem(arma::find_finite(x)).ones();
    if (ymin < ymax) y = (y-ymin) / (ymax-ymin); else y.elem(arma::find_finite(y)).ones();
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
  if (!global) {
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
#pragma omp parallel for
  for(int loc = 0; loc < nlocs; ++loc){
    // decalre inside as private variables for parallel
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
  else Rcpp::stop("stat should be 'xmn' or 'ymn', 'xsd', 'ysd'!");
}


//' Focal statistics on temporal window
//'
//' This function computes the GCSM at each cell treating the 3rd dimension
//' as a temporal window.
//' @inheritParams focal_stat_2d
//' @param x A 3d array. The 3rd one is time.
//' @param y The other array.
//' @return A matrix of \code{stat}.
//' @export
// [[Rcpp::export]]
arma::mat focal_stat_3d(arma::cube x, arma::cube y, bool fuzzy,
                        double xmin = NA_REAL, double xmax = NA_REAL,
                        double ymin = NA_REAL, double ymax = NA_REAL,
                        Rcpp::String stat = "xmn"){
  if (x.has_nan() || y.has_nan()){
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
  }
  if (!is_finite(xmin)) xmin = x.min();
  if (!is_finite(xmax)) xmax = x.max();
  if (!is_finite(ymin)) ymin = y.min();
  if (!is_finite(ymax)) ymax = y.max();
  // x or y may be a constant or NA
  if (xmin < xmax) x = clamp(x, xmin, xmax);
  if (ymin < ymax) y = clamp(y, ymin, ymax);
  if (fuzzy) {
    if (xmin < xmax) x = (x-xmin) / (xmax-xmin); else x.elem(arma::find_finite(x)).ones();
    if (ymin < ymax) y = (y-ymin) / (ymax-ymin); else y.elem(arma::find_finite(y)).ones();
  }

  arma::mat xmn = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat ymn = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat xsd = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  arma::mat ysd = arma::mat(x.n_rows,x.n_cols).fill(datum::nan);
  // uvec locs = find_finite(x.slice(0));
  // int nlocs = locs.n_elem;
  int nlocs = x.slice(0).n_elem;
#pragma omp parallel for
  for(int loc = 0; loc < nlocs; ++loc){
    // decalre inside as private variables for parallel
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
  else Rcpp::stop("stat should be 'xmn' or 'ymn', 'xsd', 'ysd'!");
}

