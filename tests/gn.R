## gn

library(cartogramR)
data(usa)
oldoptions <- options(digits=5)
usa2 <- usa[usa$region=="4",]
carto <-  cartogramR(usa2, "electors64", method="gn", options=list(L=128, relerror=1.5))
carto$final_centers
carto$final_area
carto$cartogram[[10]]
summary(carto)
options(oldoptions)
