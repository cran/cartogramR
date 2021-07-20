## Installation

Install cartogramR from CRAN with:

``` r
install.packages("cartogramR")
```
For windows users, the CRAN distributes a compiled package and there is no need to install it from source. If you install package from sources, you will need [FFTW](http://www.fftw.org/). This library is fairly common thus a package for it is usually available (linux: see fftw-devel (rpm), fftw-dev (deb); Mac OS X: fftw (brew)).
## Usage
 
For example cartogram for the number of electors in the USA in 1964:
1. Load package, data and verify that at least one grid point fall in
   each region (the smallest region is Washington DC)

   ``` r
   library(cartogramR)
   data(usa)
   plot(precartogramR(usa))
   ```
2. The default grid size $L=512$ is OK, thus we can run the cartogram for 
   variable `electors64` and plot the cartogram:
   ``` r
   carto <- cartogramR(usa, "electors64")
   plot(carto)
   ```


## Acknowledgements
  - Many many thanks to [Michael T. Gastner](https://www.yale-nus.edu.sg/about/faculty/michael-t-gastner/) for his
   [flow based cartogram](https://github.com/Flow-Based-Cartograms/go_cart) programs on which cartogramR is (heavily) based
  -  Timothée GIRAUD UMS 2414 RIATE / CNRS Paris, for suggestions and careful reading
