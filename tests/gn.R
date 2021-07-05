## gn

library(cartogramR)
data(usa)
options(digits=5)
carto <-  cartogramR(usa, "electors64", method="gn", options=list(L=128, relerror=2))
carto$final_centers
carto$final_area
carto$cartogram[[15]]
summary(carto)
