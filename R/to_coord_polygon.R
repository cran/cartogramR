#' Transform from coordinates system used in flow based cartogram to coordinates
#' system used in the polygons
#' 
#' Apply the mapping from the coordinates system used in flow based cartogram
#'   to the coordinates system used in the polygons (caracterised by the CRS)
#'
#' @param coord a vector of length 2 or a two columns matrix containing
#'    xy  coordinates to transform
#' @param carto a cartogramR object
#'
#' @return a vector of length 2 or a two columns matrix containing xy
#'   coordinates in the coordinate systems of polygons used to build
#'   the cartogram
#'
#' @export
#' @examples
#'  \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   to_coord_polygon(c(256,256), carto)
#'  }
#'
#'
#' @md
to_coord_polygon <- function(coord, carto) {
    if (!inherits(carto, "cartogramR"))
        stop("carto not inherits from cartogramR class")
    if (!is.numeric(coord))
        stop("coordinates must be numeric")
    if (is.vector(coord)) {
        if (length(coord)!=2) 
            stop("coordinates must be a vector of length 2 or a two columns matrix containing xy coordinates")
        coord <- matrix(coord, ncol=2, nrow=1)
    }
    if (is.matrix(coord)) {
        if (ncol(coord)!=2) 
            stop("coordinates must be a vector of length 2 or a two columns matrix containing xy coordinates")
    }
    LL <-  carto$options$paramsint[1]
    padding <- carto$options$paramsdouble[3]
    bbox <- sf::st_bbox(carto$initial_data)
    Delta <- c(diff(bbox[c(1,3)]), diff(bbox[c(2,4)]))
    gg <- c(sum(bbox[c(1,3)]), sum(bbox[c(2,4)]))/2
    mm <- gg - Delta/2 * padding
    MM <- gg + Delta/2 * padding
    biggest <- which.max(MM-mm)
    smallest <- which.min(MM-mm)
    scale <- (MM[biggest]-mm[biggest])/LL
    newmmB <- gg[biggest] - 0.5*LL*scale
    lxy <- 2^(ceiling(log2((MM[smallest]-mm[smallest])/scale)))
    newmmS <- gg[smallest] - 0.5*lxy*scale
    coord <- sweep(coord,2,rep(scale,2),FUN="*")
     if (biggest==1) {
         coord <- sweep(coord,2,c(newmmB,newmmS),FUN="+")
     } else {
         coord <- sweep(coord,2,c(newmmS,newmmB),FUN="+")
     }
    return(coord)
}
