## dcn

library(cartogramR)
data(usa)
options(digits=5)
precarto <- precartogramR(usa, method="dcn")
summary(precarto)

carto <-  cartogramR(usa, "electors64", method="dcn", options=list(relerror=4))
carto$final_centers
carto$final_area
carto$cartogram[[15]]
summary(carto)

library(cartogramR)
data(france_dept)
carto2 <- cartogramR(france_dept, "n_physicians", method="dcn", options=list(relerror=1e-5))
fr <- sf::st_geometry(france_dept)
plot(fr[c(92,75,78)+1])
plot(carto2$cartogram[c(92,75,78)+1])
