#' Apply the deformation used to build a cartogram to a set of
#' simple geometry  coordinates
#'
#' Please use [cartogramR::warp_features()] instead.
#'
#' @param sfgeom a sf or a sfc object which contains simple feature geometry of
#'     types in the following `POINT`, `MULTIPOINT`, `LINESTRING`,
#'     `MULTILINESTRING`, `POLYGON`, `MULTIPOLYGON
#' @param carto a cartogramR object
#' @param verbose a boolean object to set on verbose mode (default to `FALSE`)
#'
#' @return a sf or a sfc object which contains simple feature geometry transformed
#'
#'
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   LA <- sf::st_sfc(sf::st_point(c(-118.243685, 34.052234)))
#'   sf::st_crs(LA) <- 4326
#'   moregeom <- geom_cartogramR(LA, carto)
#'   plot(carto)
#'   plot(moregeom, add=TRUE, col=2, pch=15)
#' }
#'
#' @rdname cartogramR-deprecated
#'
#' @export
geom_cartogramR <- function(sfgeom, carto, verbose=FALSE) {
  .Deprecated("warp_features")
  warp_features(sfgeom, carto, verbose)
}
