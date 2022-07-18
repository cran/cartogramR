#' Errors of a cartogram object
#'
#' @param object a cartogramR object
#' @param \\dots arguments passed to or from other methods. The following
#'   arguments are available:
#'     - type; a character string giving the type of residuals (see details)
#'       - "relative error"
#'       - "absolute error"
#'       - "symmetric difference"
#'     - center ; a character string giving the type of center:
#'       - "point_on_surface" ([st_point_on_surface] applied on original
#'          and on deformed/cartogram region).
#'       - "deformed_center" (the center function, see
#'          [cartogramR_options], is applied on region and
#'          this center follows the deformation  giving
#'          the center on the deformed/cartogram region)
#'       - "centroid" (centroid of original and
#'          deformed/cartogram region).
#'       Can be abbreviated.
#'
#' @return A numeric vector which contains for each region observed
#'     area minus theorical area
#' @details The error vector contains the values of the differences
#'     between actual area of regions in the cartogram and theorical area
#'     (obtained with conservation of total area and constant density over
#'     region in the final cartogram)
#'
#'     Relative error are the error vector divided by the theorical area
#'
#'     Symmetric difference are the symmetric difference between
#'     actual area of regions in the cartogram and the original
#'     area. Each region is scaled to have an area equal to 1 and centered
#'     around the chosen center.
#'
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   residuals(carto)
#' }
#'
#' @md
residuals.cartogramR <- function(object, ...) {
  if  (!inherits(object, "cartogramR")) stop(paste(deparse(substitute(object)), "must be a cartogramR object"))
  default_options <- list(type="relative error", center = "point_on_surface")
  argsup <- eval(substitute(alist(...)))
  req_options <- names(argsup)
  default_options[req_options] <- argsup
  type <- match.arg(default_options$type, choices=c("relative error", "absolute error", "symmetric difference"))
  if (type=="absolute error") {
    r <- object$final_area - object$target_area
    return(r)
  }
  else if (type=="relative error") {
    rr <- (object$final_area/object$target_area -1)
    return(rr)
  }
  else if (type=="symmetric difference") {
    center <- match.arg(default_options$center, choices=c("point_on_surface", "centroid", "deformed_center"))
    if ((center=="deformed_center") &
        (object$details["method"]=="dcn" )) stop ('center type "deformed_center" not available for DCN algorithm')
    n_reg <- length(object$orig_area)
    delta <- rep(0,n_reg)
    y_geom <- sf::st_geometry(object$initial_data)
    standardize <- function(mpoly,center) {
      center <- sf::st_point(center)
      area <- sf::st_area(mpoly)
      (mpoly - center)/sqrt(area)
    }
    if (sf::sf_extSoftVersion()["GEOS"] < "3.8.0") {
      if (!requireNamespace("lwgeom", quietly = TRUE))
        stop("update GEOS (version >= 3.8.0 needed) OR install package lwgeom")
      make_valid_sfg <- function(x) {
        valpol <- lwgeom::lwgeom_make_valid(sf::st_geometry(x))[[1]]
      }
    } else {
      make_valid_sfg <- function(x) {
        valpol <- sf::st_make_valid(x)
      }
    }
    for (i in 1:n_reg) {
      center_orig <-
        switch(center,
               point_on_surface=sf::st_coordinates(
                                      sf::st_point_on_surface(
                                            make_valid_sfg(sf::st_geometry(object$initial_data)[[i]])))[1,],
               deformed_center=object$orig_centers[i,],
               centroid=sf::st_coordinates(
                              sf::st_centroid(
                                    make_valid_sfg(sf::st_geometry(object$initial_data)[[i]])))[1,])

               center_final <-
        switch(center,
               point_on_surface=sf::st_coordinates(
                                      sf::st_point_on_surface(
                                            make_valid_sfg( object$cartogram[[i]])))[1,],
               deformed_center=object$final_centers[i,],
               centroid=sf::st_coordinates(
                              sf::st_centroid(
                                    make_valid_sfg(object$cartogram[[i]])))[1,])
      original_mp <- standardize(y_geom[[i]],center_orig)
      final_mp <- standardize(object$cartogram[[i]],center_final)
      valpol <- make_valid_sfg(final_mp)
      delta[i] <- 1 - sf::st_area(sf::st_intersection(original_mp, valpol))
    }
    return(delta)
  }
}
