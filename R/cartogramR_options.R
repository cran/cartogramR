#' Set the options of [cartogramR] in the correct format
#'
#' @param options a named list with some (or all) the following components:
#' - absrel:  (all method)  boolean, if `TRUE` relative convergence
#'      if `FALSE` absolute convergence (default to `TRUE`)
#' - abserror: (all method) Areas on cartogram differ at most by an
#'     (absolute value of) error of abserror. That is,
#'   \eqn{max_{polygons} |area_on_cartogram - target_area| <= abserror}
#'   (default to 10000)
#' - abstol:  (`"dcn"`) the absolute convergence error tolerance:
#'    \eqn{max_\{polygons\} |area(i) - area(i-1)|}
#'    default to 1000
#' - relerror: (all method) Areas on cartogram differ at most by an
#'     (absolute value of) relative error of relerror. That is,
#'   \eqn{max_\{polygons\} |area_on_cartogram / target_area - 1| <= relerror}
#'   (default to 0.01)
#' - reltol:  (`"dcn"`) the absolute convergence tolerance:
#'    \eqn{max_\{polygons\} abs(area(i) - area(i-1))/area(i-1)}
#'    default to 1e-3
#' - L: (`"gsm" or "gn"`) integer, gives the value of `L` (default
#'     is 512), must be a power of two (for fftw)
#' - maxit:  (all method) the maximum number of iterations,
#'       default to 50 for `"gsm" or "gn"` and 5000 for `"dcn"`.
#' - maxit_internal:  (`"gsm" or "gn"`) the maximum number of internal
#'       iterations, default to 10000.
#' - min_deltat:  (`"gsm" or "gn"`) the time step minimum, default to 1e-14.
#' - mp: (all method) if a region contains exactly zero population, it will be
#'   replaced by mp times the smallest (strictly) positive population in any
#'   region (default to 0.2)
#' - pf: (`"gsm" or "gn"`) Determines space between map and boundary (default to 1.5)
#' - sigma: (`"gsm" or "gn"`) Width of Gaussian blur to smoothen the density (default to 5)
#' - maxinc: (`"gsm" or "gn"`) number of times algorithm is allowed to
#'      increase (default to 3)
#' - center: (`"gsm" or "gn"`) either a character string
#'    (only possible choices are `"centroid"`
#'     or `"point_on_surface"`) or a function. If the
#'    object is a function, it  will be used to
#'    calculate the "center" of polygons; `"point_on_surface"`
#'      will use the function [sf::st_point_on_surface]
#'      while `"centroid"` (the default) will use [sf::st_centroid].
#' - grid: (`"gsm" or "gn"`) boolean, if `TRUE` export the final
#'      grid from flow algorithm (default to `TRUE`). Setting to `FALSE`
#`      reduce memory allocation but it will be impossible to add further
#`      geometry to the cartogram (with [carto_geom])
#' - check.ring.dir: (all method) boolean, if `TRUE` controls polygons orientation
#'   (default to `TRUE`)
#' - check.only: (all method) boolean, if `TRUE` control only polygons orientation
#' - check.duplicated: (all method) boolean, if `TRUE` controls
#'      if polygons have duplicated points (default to `TRUE`)
#' - verbose: (all method) integer giving the verbosity level
#'           (default to `0`, not verbose)
#' - num_threads: (all method) int, number of threads (default to 1).
#'                For `"dcn"` used for centroid calculations.
#'                If multipolygons/regions have huge
#'                 number of points it can be useful to set it to a small
#'                 number (2 or 3 ?). For `"gsm"` useful for large grid.
#'                 `num_threads=-1` means that open MP will choose the
#'                 number of threads using maximum.
#' @param method the method to be used, can be one of the following:
#'        `gsm` or `GastnerSeguyMore` (default), `gn` or
#'        `GastnerNewman`, `dcn` or `DougenikChrismanNiemeyer`.
#' @return a list to be processed by [cartogramR]
#' @export
#'
#' @examples
#' \donttest{
#'   data(usa)
#'   carto1 <- cartogramR(usa, "electors64", options=list(verbose=1, L=256))
#'   plot(carto1)
#' }
#'
#' @references
#'
#'  * Dougenik, J., Chrisman, R. &  Niemeyer, D. (1985).
#'    An algorithm to construct continuous area cartograms.
#'    Professional Geographer **37**: 75-81.
#'  * Gastner, M. & Newman, M. E. J. (2004). Diffusion-based
#'    method for producing density equalizing
#'    maps. *Proc. Natl. Acad. Sci. USA*, **101**:7499-7504
#'  * Gastner, M., Seguy, V. & More, P. (2018). Fast flow-based
#'    algorithm for creating density-equalizing map
#'    projections. *Proceedings of the National Academy of Sciences
#'    USA*, **115**:E2156-E2164
#'
#' @md
cartogramR_options <- function(options,
                               method = c("gsm", "gn", "dcn",
                                          "GastnerSeguyMore", "GastnerNewman",
                                          "DougenikChrismanNiemeyer")) {
  method <- match.arg(method)
  if (method=="DougenikChrismanNiemeyer") method <- "dcn"
  if (method=="GastnerNewman") method <- "gn"
  if (method=="GastnerSeguyMore") method <- "gsm"
  resd <- list("relerror"=1e-2, "reltol"=1e-3,
               "abserror"=1e4, "abstol"=1e3, absrel=TRUE, "L"=512L,"maxit"=ifelse(method=="dcn", 5000L, 50L), "maxit_internal"=10000L, "min_deltat"=1e-14,
               "mp" = 0.2, "pf"=1.5 , "sigma"= 5, "maxinc"=3, "grid"=TRUE, check.ring.dir=TRUE, check.only=FALSE, check.duplicated=TRUE,
               center="centroid", "verbose"=0,
               num_threads=1)
  select <- names(options) %in% names(resd)
  if (!all(select)) stop(paste("the following name(s)",names(options)[!select],"does/do not match the options name"))
  resd[names(options)] <- options
  if ((!is.character(resd$center))&&(!is.function(resd$center))) stop("the center component must be either a character or a function that aims to calculate a center for each multipolygon")
  if (is.character(resd$center)) {
    prov <- match.arg(resd$center,c("centroid","point_on_surface"))
    if (prov=="centroid")
      resd$center <- sf::st_centroid else
                                               resd$center <- sf::st_point_on_surface
  }
  if (all(unlist(lapply(resd[c("maxit","maxit_internal","min_deltat", "relerror","reltol","abserror","abstol","L", "mp", "pf", "sigma", "maxinc", "verbose", "num_threads")],is.numeric)))) {
    if (resd$maxit<1) stop("at least 1 iteration")
    if (resd$maxit_internal<2) stop("at least 2 for max of iterations")
    if (resd$min_deltat<=0) stop("delta_t must be positive (and usually small)")
    if (resd$abserror<=0) stop("absolute error must be strictly positive")
    if (resd$abstol<=0) stop("absolute tolerance must be strictly positive")
    if (resd$relerror<=0) stop("relative error must be strictly positive")
    if (resd$reltol<=0) stop("relative tolerance must be strictly positive")
    if (log2(resd$L)!=floor(log2(resd$L))) {
      stop("L must be a power of two")
    }
    if (resd$mp<=0) stop("mp must be positive")
    if (resd$pf<1) stop("pf must greater or equal to 1")
    if (resd$sigma<=0) stop("sigma must be positive")
    if (resd$maxinc<=0) stop("maxinc must be positive")
    if (resd$verbose<0) stop("verbose must be non negative")
    if ((resd$num_threads==0)|(resd$num_threads< -1)) stop("num_threads must be either greater or equal to 1 or equal to -1 (automatic)")
  } else {
    stop("one or several parameters in the following list:\nmaxit, maxit_internal, relerror, reltol, abserror, abstol,  L,  mp, pf, sigma, verbose, num_threads\nare not numeric")
  }
  if (!is.logical(resd$absrel)) stop("absrel must be boolean")
  if (!is.logical(resd$grid)) stop("grid must be boolean")
  if (!is.logical(resd$check.ring.dir)) stop("check.ring.dir must be boolean")
  if (!is.logical(resd$check.only)) stop("check.only must be boolean")
  if (!is.logical(resd$check.duplicated)) stop("check.duplicated must be boolean")
  if (method=="dcn") {
    resd$center <- sf::st_centroid
    return(list(paramsdouble=c(ifelse(resd$absrel, resd$relerror, resd$abserror),
                               ifelse(resd$absrel, resd$reltol, resd$abstol),
                               mp=resd$mp),
                paramsint=c(maxit=resd$maxit),
                options=c(verbose=resd$verbose,
                          absrel=as.integer(resd$absrel),
                          num_threads=as.integer(resd$num_threads),
                          maxinc=as.integer(resd$maxinc)),
                check.ring.dir=resd$check.ring.dir, check.only=resd$check.only, center=resd$center,
                check.duplicated=resd$check.duplicated))
  } else {
    return(list(paramsdouble=c(ifelse(resd$absrel, resd$relerror,
                                             resd$abserror),
                                      mp=resd$mp,
                                      pf=resd$pf,sigma=resd$sigma),
                       paramsint=c(resd$L,maxit=resd$maxit),
                       options=c(verbose=resd$verbose,
                                 diff=ifelse(method=="gn",1,0),
                                 gridexport=as.numeric(resd$grid),
                                 absrel=as.integer(resd$absrel),
                                 maxitint=as.integer(resd$maxit_internal),
                                 mindeltat=as.integer(-log10(resd$min_deltat)),
                          num_threads=as.integer(resd$num_threads)),
                       check.ring.dir=resd$check.ring.dir, check.only=resd$check.only, center=resd$center, check.duplicated=resd$check.duplicated))
  }
}
