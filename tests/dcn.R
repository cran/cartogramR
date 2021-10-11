## dcn

library(cartogramR)
data(usa)
oldoptions <- options(digits=5)
usa2 <- usa[usa$region=="4",]
precarto <- precartogramR(usa2, method="dcn")
summary(precarto)

carto <-  cartogramR(usa2, "electors64", method="dcn", options=list(relerror=2))
carto$final_centers
carto$final_area
carto$cartogram[[10]]
summary(carto)
options(oldoptions)
