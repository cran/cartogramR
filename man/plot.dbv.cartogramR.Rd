% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.dbv.cartogramR.R
\name{plot.dbv.cartogramR}
\alias{plot.dbv.cartogramR}
\title{Plot a dbv.cartogram object}
\usage{
\method{plot}{dbv.cartogramR}(x, which = 1:2, ask = TRUE, key = TRUE, last = 10, probminx = 0.9, ...)
}
\arguments{
\item{x}{a dbv.cartogram object}

\item{which}{if a subset of the plots is required, specify a subset of the
numbers \code{1:2}}

\item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see
\link{par}\code{(ask=.)}}

\item{key}{logical; if \code{TRUE}, a legend is drawn}

\item{last}{draw the density of distance between vertices for the \code{last}
coordinates}

\item{probminx}{the sample quantiles (of distance between vertices)
corresponding to the probability is used as a minimum
of x-axis for the density plot (used only if \code{last} is \code{NULL})}

\item{\\dots}{arguments passed to or from other methods.}
}
\value{
No return value, called for side effects
}
\description{
Plot a dbv.cartogram object
}
\details{
The first plot is the density of distance between
consecutive vertice by region. Only the upper quantiles
are shown. The second plot is a barplot by region of the
number of vertice divided by the perimeter of the region
}
\examples{
  data(usa)
  precarto <- precartogramR(usa, method="dcn")
  plot(precarto)

}
