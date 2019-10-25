#' get and set attribute of stars object
#'
#' \code{get_attr} retrieves an attribute from stars, and rearranges to wide-form.
#' \code{set_attr}, on the contrary, puts a wide-form matrix to stars as an attribute.
#' @name get_set_attr
#' @param x object of class \code{stars}
#' @param attr charter indicating attribute
#' @importFrom stars st_dimensions
#' @importFrom magrittr %>%
#' @return \code{get_attr} returns a wide-from matrix
#' @export
#' @examples
#' library(stars)
#' x <- matrix(1:18, nrow=6)
#' st <- x %>% st_as_stars() %>% setNames('v')
#' st %>% image(text_values = TRUE)
#' st %>% get_attr('v')
get_attr <- function(x, attr) {
  ar_long <- x[[ attr ]]
  attr(dim(ar_long), 'names') <- NULL
  ndims <- length(dim(ar_long))
  y_is_neg <- st_dimensions(x)[[2]]$delta < 0

  if (ndims == 2) {
    if (y_is_neg) ar_long %>% aperm()
    else ar_long %>% aperm() %>% apply(2, rev)
  }
  else if (ndims == 3) {
    if (y_is_neg) ar_long %>% aperm(c(2,1,3))
    else ar_long %>% aperm(c(2,1,3)) %>% apply(2:3, rev)
  }
  else stop('Only 2-d or 3-d array are supported')
}

#' @name get_set_attr
#' @param ar_wide a wide-from matrix
#' @return \code{set_attr} returns a stars object
#' @export
set_attr <- function(x, attr, ar_wide) {
  attr(dim(ar_wide), 'names') <- NULL
  ndims <- length(dim(ar_wide))
  y_is_neg <- st_dimensions(x)[[2]]$delta < 0

  if (ndims == 2){
    if (y_is_neg) ar_long <- ar_wide %>% aperm()
    else ar_long <- ar_wide %>% apply(2, rev) %>% aperm()
  }
  else if (ndims == 3){
    if (y_is_neg) ar_long <- ar_wide %>% aperm(c(2,1,3))
    else ar_long <- ar_wide %>% apply(2:3, rev) %>% aperm(c(2,1,3))
  }
  else stop('Only 2-d or 3-d array are supported')
  dim(ar_long) <- dim(x)
  x[[ attr ]] <- ar_long
  x
}
