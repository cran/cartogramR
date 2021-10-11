## gsm

library(cartogramR)
data(usa)
usa2 <- usa[usa$region=="4",]
oldoptions <- options(digits=5)
precarto <- precartogramR(usa2, gridpower2=5:8)
summary(precarto)

carto <-  cartogramR(usa2, "electors64", method="gsm", options=list(L=128, relerror=1.5))
carto$final_centers
carto$final_area
carto$cartogram[[10]]
summary(carto)
options(oldoptions)
