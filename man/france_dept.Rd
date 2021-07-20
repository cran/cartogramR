\name{france_dept}
\docType{data}
\alias{france_dept}
\title{Map of the population of mainland France (year 2018)}
\description{
  This data set is a basemap of mainland France with the population in 2018 (\code{pop2018}), the number of physicians in 2018 (\code{n_physicians}), the number of general practitioner  in 2018 (\code{n_gp})  and the number of general practitioner for 100000 inhabitants  in 2018 (\code{n_gp_per100000}) in each department (\code{dept_name} or \code{id}).
}
\usage{data(france_dept)}
\source{\url{https://www.data.gouv.fr/fr/datasets/admin-express/#_}}
        \references{
           \itemize{
\item \url{https://www.insee.fr/fr/statistiques/2012713#tableau-TCRD_004_tab1_departements}
\item \url{https://www.insee.fr/fr/statistiques/2012677#tableau-TCRD_068_tab1_departements}}}
\format{A sf object containing 8 columns of data and the geometry
        Projected CRS: RGF93 / Lambert-93 (EPSG : 2154)}
\keyword{datasets}
