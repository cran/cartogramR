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
    namevarall <- names(x$initial_data)
    namegeom <- attr(x$initial_data, "sf_column")
    namesvar <- namevarall[namevarall!=namegeom]
    data <- data.frame(x$initial_data)[,namesvar]
    y <- cbind.data.frame(data, x$orig_area, x$final_area, x$target_area,
                          x$orig_centers, x$final_centers)
    p <- ncol(y)
    names(y)[(p-6):p] <- c("orig_area", "final_area", "target_area", "x_orig_centers",
                    "y_orig_centers", "x_final_centers", "y_final_centers")
    y <- sf::st_as_sf(cbind.data.frame(y,  geometry=x$cartogram))
    return(y)
}
