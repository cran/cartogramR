## gsm

library(cartogramR)
data(usa)
options(digits=5)
precarto <- precartogramR(usa, gridpower2=5:8)
summary(precarto)

carto <-  cartogramR(usa, "electors64", method="gsm", options=list(L=128, relerror=2))
carto$final_centers
carto$final_area
carto$cartogram[[15]]
summary(carto)
