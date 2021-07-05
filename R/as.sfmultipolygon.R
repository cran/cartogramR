#' Transform a sf object with several rows (polygons) by region to an
#' sf object with one row by region and thus one multipolygon by region
#'
#' @param data a sf object
#' @param idregion a character string which indicates the name of the
#'     column (in `data` object) which contains the region
#'     identifier.
#' @param closepolygon a boolean (default to `FALSE`) if `TRUE` it
#'     controls if polygons are closed and if not add the first
#'     vertice at the end.
#'
#' @return a sf object with one row by region and one multipolygon by
#'     region.
#'
#' @export
#'
#' @md
as.sfmultipolygon <- function(data, idregion, closepolygon=FALSE) {
  if (! inherits(data,"sf")) stop(paste(deparse(substitute((data))), "is not an object of class sf"))
  dataf <- sf::st_drop_geometry(data)
  if (! idregion%in%names(dataf)) stop(paste("no variable R",deparse(substitute((idrgion))),"in sf",deparse(substitute((data)))))
  region <- dataf[,idregion]
  ureg <- unique(region)
  final <- as.list(ureg)
  ## loop on region (unique)
  for (i in 1:length(ureg)) {
    ## extraction of data and coordinates
    don <- data[region==ureg[i],]
    coord <- sf::st_coordinates(don)
    ## separation
    if (any(names(coord)=="L3")) {
      ff <- paste(coord[,"L1"], coord[,"L2"], coord[,"L3"], sep="_")
    } else {
      ff <- paste(coord[,"L1"], coord[,"L2"], coord[,"L3"], sep="_")
    }
    liste <- split(coord[,1:2],list(factor(ff)))
    ## transform to  multipolygon
    if (closepolygon) {
      final[[i]] <- sf::st_multipolygon(lapply(liste,FUN=function(x) {
        m <- matrix(x,ncol=2)
        if (any(m[1, ] != m[nrow(m), ]))
          m <- rbind(m, m[1,,drop=FALSE])
        list(m) }))
    }  else {
      final[[i]] <- sf::st_multipolygon(lapply(liste,FUN=function(x) list(matrix(x,ncol=2))))
    }
  }
  ## to sfc
  datagm <- sf::st_as_sfc(final)
  ## remove duplicated in data
  newdataf <- dataf[!duplicated(region),]
  ## to sf
  datam <- sf::st_sf(data.frame(newdataf,geometry=datagm))
  ## CRS
  sf::st_crs(datam) <- sf::st_crs(data)
  return(datam)
}
