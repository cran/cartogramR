#' Coerce a cartogramR to a sf object
#'
#' Coerce a cartogramR to a sf object returning the sf object used to
#'   construct the cartogram with the cartogram as geometry and some
#'   more attributes
#'
#' @param x a cartogramR object
#' @param   \\dots arguments passed to or from other methods.
#'
#' @return a sf object including all the data (attributes) contained
#'   in the original sf object used to construct the cartogram and
#'   - original areas of region (`orig_area`)
#'   - final/deformed areas of region (`final_area`)
#'   - target areas of region (`target_area`)
#'   - original centers (`x_orig_centers` and `y_orig_centers`)
#'   - final centers (`x_final_centers` and `y_final_centers`)
#'
#' @rdname as.sf.cartogramR
#' @export
#'
#' @md
as.sf.cartogramR <- function(x, ...) {
    if  (!inherits(x, "cartogramR")) stop(paste(deparse(substitute(x)), "must be a cartogramR object"))
    return(x$cartogram)
}
