#' Make a continuous cartogram (density equalizing maps)
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
#' @param count a character string which indicates the name of the
#'     column (in `data` object) which contains the count by
#'     region.
#' @param method the method to be used, can be one of the following:
#'        `gsm` or `GastnerSeguyMore` (default), `gn` or
#'        `GastnerNewman`, `dcn` or `DougenikChrismanNiemeyer`.
#' @param options a named list given to [cartogramR_options] function
#'     which process options see [cartogramR_options] for
#'     details. Default to `NULL`.
#' @return A list with the following components:
#' - cartogram: a sf object (in the same order of `data` or sorted by `idregion`
#'    see reordered argument) which contains the cartogram
#'    (ie the initial polygons after deformation)
#' - orig_area: original areas of regions
#' - final_area: final areas of regions in the cartogram
#' - orig_centers: the initial centers calculated with [st_point_on_surface]
#' - final_centers: the centers after deformation
#' - gridx: (for flow-based method) final grid (x-axis) if requested
#'         (see [cartogramR_options] for details).
#' - gridy: (for flow-based method) final grid (y-axis) if requested
#'         (see [cartogramR_options] for details).
#' - count: the count by region
#' - target_area: target areas of regions
#' - initial_data: the initial sf object
#' - details: names of original data, idcount variable, algorithm
#' - options: values of options
#'
#' @export
#' @import data.table
#' @useDynLib cartogramR, .registration = TRUE, .fixes = "carto_"
#' @examples
#' \donttest{
#'   data(usa)
#'   carto <- cartogramR(usa, "electors64")
#'   plot(carto)
#'   summary(carto)
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
#'    USA*, **115**:E2156-E2164, website:
#'   [go-cart](https://go-cart.io/)
#'
#' @md
cartogramR <- function(data, count,
                       method=c("gsm", "gn", "dcn",
                                "GastnerSeguyMore", "GastnerNewman", "DougenikChrismanNiemeyer"), options=NULL) {
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
n_reg <- nrow(data)
y_geom <- sf::st_geometry(data)
if (length(count)>1) stop("Only one variable for count")
if (is.character(count)) {
  if(is.na(match(count,names(data)))) stop("Names of count variable is wrong")
} else {
  if (!(count%in%(1:ncol(data)))) stop("Column number for count variable is wrong")
}
countVar <- as.double(data[[count]])
if (!is.numeric(countVar) || any(countVar<0) || any(is.na(countVar)))
  stop("Count must be numeric, non negative, NA are not allowed")
    multipolygons <- rep(NA, n_reg)
    sp <- unlist(lapply(y_geom, function(x) inherits(x, "POLYGON")))
    multipolygons[sp] <- 0L
    mp <- unlist(lapply(y_geom, function(x) inherits(x, "MULTIPOLYGON")))
    multipolygons[mp] <- 1L
    if (any(is.na(multipolygons))) stop("one or more geometries are neither
       POLYGON nor MULTIPOLYGON")
    if ((method=="dcn") & (any(multipolygons==0)))
      stop("Dougenik, Chrisman & Niemeyer algorithm requires multipolygons")
    currentoptions <- cartogramR_options(options, method)
    if (currentoptions$check.ring.dir) {
        if (currentoptions$check.only) {
            if (any(!check_ring_dir(y_geom, currentoptions$check.only)))
              stop("at least one polygon is not oriented correctly, please use check.only=FALSE to correct")
        } else {
            y_geom <- check_ring_dir(y_geom, currentoptions$check.only)
        }
    }
    centers <- unlist(lapply(y_geom, currentoptions$center))
    centersx <- centers[seq(1, n_reg*2 - 1, by=2)]
    centersy <- centers[seq(2, n_reg*2, by=2)]
    minim <- min(countVar[countVar > 0])* currentoptions$paramsdouble["mp"]
    trimCountVar <- countVar
    trimCountVar[trimCountVar==0] <- minim
if (method=="dcn") {
  L1 <- L2 <- L3 <- X <- Y <- index <- rdup <- rlast <- NULL # due to NSE notes in R CMD check
  coord <- unique(data.table::data.table(data.frame(sf::st_coordinates(y_geom))))
  coord[,":="("index"= .I - 1)]
  fperm <- function(x) if (length(x)>1) return(c(x[-1],x[1])) else -1
  flast <- function(x) {
    nx <- length(x)
    return(rep(x[nx],nx)) }
  coord[,":="("rdup"=fperm(index), "rlast"=flast(index)),by=list(X,Y)]
  paramsint <- currentoptions$paramsint
  paramsdouble <- currentoptions$paramsdouble
  options <- currentoptions$paramsint
  results <- .Call(carto_dcn, y_geom, coord[,X], coord[,Y],
                   trimCountVar, as.integer(coord[,L1]),
                   as.integer(coord[,L2]), as.integer(coord[,L3]),
                   as.integer(coord[,rdup]),
                   as.integer(coord[,rlast]), as.integer(paramsint), paramsdouble, as.integer(currentoptions$option) )
  names(results) <- c("cartogram","orig_area","final_area","orig_centers","final_centers")
} else {
    niveaux <- 1:n_reg
    ff <- factor(niveaux, levels=niveaux)
    max_id <- n_reg
    min_id <- 1
    nombb <- names(attr(y_geom,"bbox"))
    target_area <- tapply(trimCountVar,ff,unique)
    if (!is.numeric(target_area)) stop("check the variable count by region (not unique by region)")
    if (length(target_area)!=n_reg) stop("check the variable count by region (not unique by region)")
    n_polycorn <- rapply(y_geom,function(x) if (is.matrix(x)) nrow(x))
   nbep <- rep(NA, n_reg)
   ## par polygone "sf"/ligne combien de polygones deplies
   nbep[sp] <- unlist(lapply(y_geom[sp],length))
   ## combien de polygones "deplies" par multipolygone (par ligne)
   nbep[mp] <-  unlist(lapply(y_geom[mp],function(z) sum(unlist(lapply(z,length)))))
    facteur <- rep(ff, nbep)
    nb_polyinreg <- table(facteur)
    attributs <- attributes(y_geom)
    bbox <- sf::st_bbox(data)[c(1,3,2,4)]
    results <- .Call(carto_cartogramR, centersx, centersy, y_geom, as.double(target_area),
              as.integer(nb_polyinreg),  as.integer(n_polycorn),
              as.integer(c(length(facteur), n_reg, min_id, max_id)),
              bbox, currentoptions$paramsdouble, as.integer(currentoptions$paramsint),
              as.integer(currentoptions$options), as.integer(multipolygons), PACKAGE="cartogramR")
    if (currentoptions$options["gridexport"]) {
      names(results) <- c("cartogram","orig_area","final_area","orig_centers","final_centers","gridx","gridy")
      } else {
        results <- results[-(6:7)]
        names(results) <- c("cartogram","orig_area","final_area","orig_centers","final_centers")
        }
    }
    results$final_centers <- matrix(c(results$orig_centers, results$final_centers), ncol=2)
    results$orig_centers <- matrix(centers, ncol=2, byrow=TRUE)
    results$count <- trimCountVar
    results$target_area <- trimCountVar/sum(trimCountVar)*sum(results$orig_area)
    results$initial_data <- data
    algo <- switch(method,"gsm"= "Gastner, Seguy & More (2018) fast flow-based algorithm", "gn"="Gastner & Newman (2004) diffusion algorithm", "dcn" = "Dougenik, Chrisman &  Niemeyer (1985) rubber band algorithm")
    results$details <- c(initial_data_name=deparse(substitute(data)),
                         initial_count_name=count,
                         method=method,algorithm=algo)
    results$options <- currentoptions
    class(results) <- c("cartogramR", "list")
    return(results)
 }
