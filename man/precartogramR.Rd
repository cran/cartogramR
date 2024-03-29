% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precartogramR.R
\name{precartogramR}
\alias{precartogramR}
\title{Make a pre cartogram analysis}
\usage{
precartogramR(
  data,
  method = c("gsm", "gn", "dcn", "GastnerSeguyMore", "GastnerNewman",
    "DougenikChrismanNiemeyer"),
  gridpower2 = 8:11,
  pf = 1.5,
  verbose = FALSE
)
}
\arguments{
\item{data}{a sf object which contains at least two columns:
obviously a geometry column (giving the map) and a column which
contains a count by region (leading to a density by region,
density to be equalized by deformation). Each row of \code{data} is a region and
contains the simple feature geometry of type \code{POLYGON} or
\code{MULTIPOLYGON}. Polygon ring directions are not checked but
exterior ring must counter clockwise and holes clockwise (use
option \code{check_ring_dir} of \link[sf:st_read]{sf::st_read} to achieve the right
orientation of ring direction on import or use \link{check_ring_dir}
function)}

\item{method}{the method to be used, can be one of the following:
\code{gsm} or \code{GastnerSeguyMore} (default), \code{gn} or
\code{GastnerNewman}, \code{dcn} or \code{DougenikChrismanNiemeyer}.}

\item{gridpower2}{a vector of exponent (to be raised at the power of 2) that
gives the log2(size) of the grid (default to \code{8:11});
meaningful for method \code{gsm} or \code{GastnerSeguyMore} (default), \code{gn} or
\code{GastnerNewman}}

\item{pf}{Determines space between map and boundary (default to 1.5);
meaningful for method \code{gsm} or \code{GastnerSeguyMore} (default), \code{gn} or
\code{GastnerNewman}}

\item{verbose}{a boolean object to set on verbose mode (default to \code{FALSE});
meaningful for method \code{gsm} or \code{GastnerSeguyMore} (default), \code{gn} or
\code{GastnerNewman}}
}
\value{
either a \code{dbv.cartogramR} object (if method
is \code{dcn} or \code{DougenikChrismanNiemeyer}) see \link{dist_between_vertices}
for details or
a \code{gridanalysis.cartogramR} (if method
is \code{gsm} or \code{GastnerSeguyMore} (default), \code{gn} or
\code{GastnerNewman}) see \link{grid_analysis} for details
}
\description{
Make a pre cartogram analysis
}
\examples{
\donttest{
  data(usa)
  precarto <- precartogramR(usa)
  plot(precarto)
  summary(precarto)
}

}
\references{
\itemize{
\item Dougenik, J., Chrisman, R. &  Niemeyer, D. (1985).
An algorithm to construct continuous area cartograms.
Professional Geographer \strong{37}: 75-81.
\item Gastner, M. & Newman, M.E.J. (2004). Diffusion-based
method for producing density equalizing
maps. \emph{Proc. Natl. Acad. Sci. USA}, \strong{101}:7499-7504
\item Gastner, M., Seguy, V. & More, P. (2018). Fast flow-based
algorithm for creating density-equalizing map
projections. \emph{Proceedings of the National Academy of Sciences
USA}, \strong{115}:E2156-E2164
}
}
