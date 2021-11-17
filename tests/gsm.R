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

poiY <- sf::st_sfc(sf::st_point(c(-1946151, -538753)))
sf::st_crs(poiY) <- sf::st_crs(usa2)
poiYb <- geom_cartogramR(poiY, carto)
print(sf::st_coordinates(poiY))

Y <- matrix(c(-1901626,-1935314,-1929683,-1898426,-1901626,-444869,-454442,-484786,-460618,-444869),5,2)
polY <- sf::st_sfc(sf::st_polygon(list(Y)), check_ring_dir = TRUE)
sf::st_crs(polY) <- sf::st_crs(usa2)
polYb <- geom_cartogramR(polY, carto)
print(sf::st_coordinates(polYb))

Y <- matrix(c(-1944542,-1927687,-1890629,-1821631,-1799218,-1783629,-1766393,-1745082,-539743,-517879,-485451,-458782,-443586,-424518,-417270,-414516),8,2)
linY <- sf::st_sfc(sf::st_linestring(Y))
sf::st_crs(linY) <- sf::st_crs(usa2)
linYb <- geom_cartogramR(linY, carto)
print(sf::st_coordinates(linYb))

delta <- c(43000, 50000)
p1 <- matrix(c(-1944542,-1901542,-1815542,-1858542,-1901542,-1944542,-539743,-539743,-439743,-339743,-339743,-539743), ncol=2)
p2 <- matrix(c(-1901542,-1901542,-1858542,-1901542,-489743,-439743,-439743,-489743), ncol=2)
p3 <- matrix(c(-1815542,-1772542,-1772542,-1815542,-1815542,-539743,-539743,-489743,-489743,-539743), ncol=2)
p4 <- matrix(c(-1802642,-1802642,-1781142,-1781142,-1802642,-524743,-499743,-499743,-524743,-524743), ncol=2)
p5 <- matrix(c(-1815542,-1772542,-1772542,-1815542,-389743,-439743,-389743,-389743), ncol=2)
mpolY <- as.list(1:2)
mpolY[[1]] <- sf::st_multipolygon(list(list(p1,p2), list(p3,p4), list(p5)))
mpolY[[2]] <- sf::st_multipolygon(list(list(p1+matrix(4*delta,nrow(p1),ncol(p1), byrow=TRUE),
                                        p2+matrix(4*delta,nrow(p2),ncol(p2), byrow=TRUE)),
                                   list(p3+matrix(4*delta,nrow(p3),ncol(p3), byrow=TRUE),
                                        p4+matrix(4*delta,nrow(p4),ncol(p4), byrow=TRUE))))
mpolY <- sf::st_sfc(mpolY,check_ring_dir = TRUE)
sf::st_crs(mpolY) <- sf::st_crs(usa2)
mpolYb <- geom_cartogramR(mpolY, carto)
print(sf::st_coordinates(mpolYb))

options(oldoptions)
