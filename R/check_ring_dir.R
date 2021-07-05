#' Polygon rings directions are checked and corrected if asked. 
#' 
#' Polygon ring are seen from above: exterior ring counter clockwise,
#' holes clockwise
#'
#' @param polygons a sfc object which contains simple feature geometry of
#'     types `POLYGON` or `MULTIPOLYGON`
#'     
#' @param check.only a boolean which indicates if the function only
#'     checks the ring direction (`check.only=TRUE`) or checks and
#'     corrects the polygon direction (`check.only=FALSE`)
#' @return Either a logical vector which indicates if line i of polygons
#'     is in the right direction (`TRUE`) or not or the corrected sfc object
#'  
#' @export
#' @examples
#'   data(usa)
#'   all(check_ring_dir(sf::st_geometry(usa), check.only=TRUE))
#'
#'
#' @md
check_ring_dir <- function(polygons, check.only=TRUE) {
    if (!inherits(polygons, "sfc")) stop("polygons not inherits from sfc class")
    n <- length(polygons)
    multipolygons <- rep(NA, n)
    sp <- unlist(lapply(polygons, function(x) inherits(x, "POLYGON")))
    multipolygons[sp] <- 0L
    mp <- unlist(lapply(polygons, function(x) inherits(x, "MULTIPOLYGON")))
    multipolygons[mp] <- 1L
    if (any(is.na(multipolygons))) stop("one or more components are neither
       POLYGON nor MULTIPOLYGON")
    results <- .Call(carto_checkring, polygons, as.integer(multipolygons), as.integer(!check.only))
    if (check.only) return(as.logical(results)) else return(results)  
}
