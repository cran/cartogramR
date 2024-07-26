#' Print a cartogram object
#'
#' @param x a cartogramR object
#' @param   \\dots arguments passed to or from other methods.
#' @return No return value, called for side effects
#' @rdname print.cartogramR
#' @export
#' @md
print.cartogramR <- function(x, ...) {
    if  (!inherits(x, "cartogramR")) stop(paste(deparse(substitute(x)), "must be a cartogramR object"))
  cat("object of class cartogramR\n")
  cat(paste("  Distorsion of polygons contained in", attr(x, "initial_data_name"),"\n"))
  cat(paste("  based on variable", attr(x, "initial_count_name"), "\n"))
  cat(paste("  using", attr(x, "algorithm"), "\n"))
    if (attr(x, "options")$options["absrel"]==1) {
      critname <- "Max of (Abs) Relative Errors: "
      crit <- max(abs(residuals.cartogramR(x, type="relative error")))
      } else {
      critname <- "Max of (Abs) Absolute Errors: "
      crit <- max(abs(residuals.cartogramR(x, type="absolute error")))
      }
  cat(paste0("  ", critname, crit, "\n"))
  return(invisible(NULL))
}
