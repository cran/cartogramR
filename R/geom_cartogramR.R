#' Apply the deformation used to build a cartogram to a set of
#' simple geometry  coordinates
#'
#' Apply the deformation used to build a cartogram to a set of
#' simple geometry coordinates. The resulting  simple geometry object can
#' be used to add geometry features on the cartogram.
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
#' @export
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
#' @md
geom_cartogramR <- function(sfgeom, carto, verbose=FALSE) {
  if (!inherits(carto, "cartogramR"))
    stop("carto not inherits from cartogramR class")
  if (carto$options$options["gridexport"]==0)
    stop("carto does not include grid\n Rerun cartogramR with options=list(grid=TRUE)")
  if (inherits(sfgeom, "sf")) {
    sfobject <- TRUE
    sfgeomb <- sf::st_geometry(sfgeom)
  } else {
    if (!inherits(sfgeom, "sfc")) stop("sfgeom not inherits from sfc class")
    sfobject <- FALSE
    sfgeomb <- sfgeom
  }
  reqcrs <- sf::st_crs(carto$cartogram)
  geomcrs <-  sf::st_crs(sfgeomb)
  changeCRS <- FALSE
  if (reqcrs!=geomcrs) {
   if (is.na(geomcrs)) {
    warning(paste("No coordinate reference system !\n Function proceeds assuming CRS of cartogram object:",reqcrs$input,"\n Setting coordinate reference system could be a good idea, see sf::st_crs"))
    sf::st_crs(sfgeomb) <- reqcrs
   } else {
    sfgeomb <- sf::st_transform(sfgeomb, reqcrs)
    changeCRS <- TRUE
    if (verbose) warning(paste("CRS set to CRS of cartogram object:",reqcrs$input))
   }
  }
  LL <-  carto$options$paramsint[1]
  padding <- carto$options$paramsdouble[3]
  n <- length(sfgeomb)
  bbox <- sf::st_bbox(carto$initial_data)
  sfgeombbox <- sf::st_bbox(sfgeomb)
  if (!(((bbox["xmin"] <= sfgeombbox["xmin"]) & (bbox["ymin"] <= sfgeombbox["ymin"])) &
        ((sfgeombbox["xmax"]  <= bbox["xmax"]) & (sfgeombbox["ymax"] <= bbox["ymax"]))))
   stop("Some part of the sfgeom is outside the bounding box of the cartogram.")
  nombb <- names(bbox)
  multipolygons <- rep(NA, n)
  sp <- unlist(lapply(sfgeomb, function(x) inherits(x, "POINT")))
  multipolygons[sp] <- 0L
  mp <- unlist(lapply(sfgeomb, function(x) inherits(x, "MULTIPOINT")))
  multipolygons[mp] <- 1L
  sp <- unlist(lapply(sfgeomb, function(x) inherits(x, "LINESTRING")))
  multipolygons[sp] <- 2L
  mp <- unlist(lapply(sfgeomb, function(x) inherits(x, "MULTILINESTRING")))
  multipolygons[mp] <- 3L
  sp <- unlist(lapply(sfgeomb, function(x) inherits(x, "POLYGON")))
  multipolygons[sp] <- 4L
  mp <- unlist(lapply(sfgeomb, function(x) inherits(x, "MULTIPOLYGON")))
  multipolygons[mp] <- 5L
  if (any(is.na(multipolygons))) stop("one or more components are not in
       POINT, MULTIPOINT, LINESTRING, MULTILINESTRING, POLYGON, MULTIPOLYGON")
  attributs <- attributes(sfgeomb)
  results <- .Call(carto_geomcarto,
                   sfgeomb, multipolygons, carto$gridx, carto$gridy,
                   padding, LL, bbox[c(1,3,2,4)], as.integer(verbose))
  bbox <- attr(results,"bbox")
  names(bbox) <- nombb
  attributes(results) <- attributs
  attr(results,"bbox") <- bbox
  if (sfobject) {
    res <- sf::st_set_geometry(sfgeom, results)
  } else {
    return(results)
  }
}
