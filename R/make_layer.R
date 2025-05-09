#' Make a layer
#'
#' Create a sfc object containing final centers, original centers,
#'   centers displacement, original graticule or final graticule.
#'
#' @param x a cartogramR object
#' @param type  a character string giving the type of layer:
#'     - "final_centers": if `method` is `dcn`, [sf::st_centroid]
#'        is applied on deformed/cartogram region ; if `method` is
#'        `gsm` or `gn` (ie flow based), initial "centers" are calculated
#'        and the cartogram deformation is applied on
#'        these "centers" giving the final_centers.
#'     - "original_centers" if `method` is `dcn`, [sf::st_centroid]
#'        is applied on original regions); if `method` is
#'        `gsm` or `gn` (ie flow based), initial "centers" are calculated
#'        using cartogramR `center` option see [cartogramR_options].
#'     - "centers_translation" linestring giving the movement of
#'        centers due to the deformation used to have the cartogram
#'     - "final_graticule" (method `gsm` or `gn`) graticule
#'        obtained by the cartogram algorithm
#'     - "original_graticule" (method `gsm` or `gn`) graticule
#'       used by the cartogram algorithm
#'
#' @return a sfc object
#'
#' @export
#'
#' @md
make_layer <- function(x, type=c("final_centers", "original_centers", "centers_translation", "final_graticule", "original_graticule")) {
  if  (!inherits(x, "cartogramR")) stop(paste(deparse(substitute(x)), "must be a cartogramR object"))
  type <- match.arg(type)
  if (type=="final_centers") {
    return(x$final_centers)
  }
  if (type=="original_centers") {
    return(x$orig_centers)
  }
  if (type=="original_graticule") {
    if (!(attr(x, "method") %in% c("gsm", "gn"))) stop("cartogram method should be either 'gsm' or 'gn'")
    bbox <- attr(x, "initial_bbox")
    LL <-  attr(x, "options")$paramsint[1]
    pf  <- attr(x, "options")$paramsdouble[3]
    graticule <- .Call(carto_makeoriggraticule, pf, as.integer(LL), bbox)
    sf::st_crs(graticule) <- sf::st_crs(x$cartogram)
    return(graticule)
  }
  if (type=="final_graticule") {
    if (!(attr(x, "method") %in% c("gsm", "gn"))) stop("cartogram method should be either 'gsm' or 'gn'")
    bbox <- attr(x, "initial_bbox")
    if (attr(x, "options")$options["gridexport"]==0)
      stop("cartogram does not include grid\n Rerun cartogramR with options=list(grid=TRUE)")
    LL <-  attr(x, "options")$paramsint[1]
    pf  <- attr(x, "options")$paramsdouble[3]
    graticule <- .Call(carto_makefinalgraticule, pf, as.integer(LL), bbox, x$gridx, x$gridy)
    sf::st_crs(graticule) <- sf::st_crs(x$cartogram)
    return(graticule)
  }
  if (type=="centers_translation") {
    matorig <- sf::st_coordinates(x$orig_centers)[,1:2]
    matfinal <- sf::st_coordinates(x$final_centers)[,1:2]
     coordLine <-lapply(1:nrow(matorig), function(n) { sf::st_linestring(rbind(matorig[n,],matfinal[n,]))})
    movement <- sf::st_sfc(coordLine)
    sf::st_crs(movement) <- sf::st_crs(x$cartogram)
    return(movement)
  }
}
