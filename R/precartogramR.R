#' Make a pre cartogram analysis
#'
#' @param data a sf object which contains at least two columns:
#'     obviously a geometry column (giving the map) and a column which
#'     contains a count by region (leading to a density by region,
#'     density to be equalized by deformation). Each row of `data` is a region and
#'     contains the simple feature geometry of type `POLYGON` or
#'     `MULTIPOLYGON`. Polygon ring directions are not checked but
#'     exterior ring must counter clockwise and holes clockwise (use
#'     option `check_ring_dir` of [sf::st_read] to achieve the right
#'     orientation of ring direction on import or use [check_ring_dir]
#'     function)
#'
#' @param method the method to be used, can be one of the following:
#'        `gsm` or `GastnerSeguyMore` (default), `gn` or
#'        `GastnerNewman`, `dcn` or `DougenikChrismanNiemeyer`.
#' @param gridpower2 a vector of exponent (to be raised at the power of 2) that
#'    gives the log2(size) of the grid (default to `8:11`);
#'  meaningful for method `gsm` or `GastnerSeguyMore` (default), `gn` or
#'  `GastnerNewman`
#' @param pf Determines space between map and boundary (default to 1.5);
#'  meaningful for method `gsm` or `GastnerSeguyMore` (default), `gn` or
#'  `GastnerNewman`
#' @param verbose a boolean object to set on verbose mode (default to `FALSE`);
#'  meaningful for method `gsm` or `GastnerSeguyMore` (default), `gn` or
#'  `GastnerNewman`
#' @return either a `dbv.cartogramR` object (if method
#'     is `dcn` or `DougenikChrismanNiemeyer`) see [dist_between_vertices]
#'     for details or
#'     a `gridanalysis.cartogramR` (if method
#'     is `gsm` or `GastnerSeguyMore` (default), `gn` or
#'        `GastnerNewman`) see [grid_analysis] for details
#'
#' @export
#' @examples
#' \donttest{
#'   data(usa)
#'   precarto <- precartogramR(usa)
#'   plot(precarto)
#'   summary(precarto)
#' }
#'
#' @references
#'
#'  * Dougenik, J., Chrisman, R. &  Niemeyer, D. (1985).
#'    An algorithm to construct continuous area cartograms.
#'    Professional Geographer **37**: 75-81.
#'  * Gastner, M. & Newman, M.E.J. (2004). Diffusion-based
#'    method for producing density equalizing
#'    maps. *Proc. Natl. Acad. Sci. USA*, **101**:7499-7504
#'  * Gastner, M., Seguy, V. & More, P. (2018). Fast flow-based
#'    algorithm for creating density-equalizing map
#'    projections. *Proceedings of the National Academy of Sciences
#'    USA*, **115**:E2156-E2164
#'
#' @md
precartogramR <- function(data, method=c("gsm", "gn", "dcn",
                                        "GastnerSeguyMore", "GastnerNewman",
                                        "DougenikChrismanNiemeyer"),
                         gridpower2=8:11, pf=1.5, verbose=FALSE) {
    method <- match.arg(method)
    if (method=="DougenikChrismanNiemeyer") method <- "dcn"
    if (method=="GastnerNewman") method <- "gn"
    if (method=="GastnerSeguyMore") method <- "gsm"
    if (!inherits(data, "sf")) stop(paste(deparse(substitute((data))),"not inherits from sf class"))
    if (is.na(sf::st_crs(data))) {
        warning("No coordinate reference system !\n Function proceeds assuming planar coordinates\n Remarks:\n 1. This function use planar coordinates to calculate areas and thus does not give correct answers for longitude/latitude data\n 2. Setting coordinate reference system could be a good idea, see sf::st_crs")
    } else {
        if (sf::st_is_longlat(data)) stop(paste(deparse(substitute((data))),"is an unprojected map.\n This function use planar coordinates to calculate areas and thus does not give correct answers for longitude/latitude data:\n => Use \"st_transform()\" to transform longitude/latitude coordinates to another coordinate system."))
    }
    if (method=="dcn") {
        return(dist_between_vertices(data))
    } else {
        return(grid_analysis(data, gridpower2, pf, verbose))
    }
}
