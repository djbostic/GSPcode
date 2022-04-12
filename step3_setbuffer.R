# load packages
library(tidyverse) # general purpose data science toolkit
library(sp)        # spatial objects
library(sf)
library(raster)    # for raster objects
library(readr)
library(rgeos)

# take output from previous script into this one - we want domestic wells intersected with gsps
domesticwells_in_gsps

miles_to_meters <- function(x){x * 1609.34}
# township level accuracy (36 square miles)
bt <- miles_to_meters(3.5)
# calculate buffers: write a fucntion and apply it
buffer_intersect_sf <- function(x) {
  # calculate the buffer
  x = st_buffer(x, dist = bt)
  # trim the buffer to the cv
  x = st_intersection(x, gsps)
  return(x) # return result
}

blist <- st_buffer(domesticwells_in_gsps, dist=bt)
blist1 <- st_union(blist)

st_write(blist1, "Boundaries/domesticwellbuffer.shp")

ggplot() +
  geom_sf(data=blist1, col="grey10")+
  geom_sf(data=domesticwells_in_gsps, col="cornflowerblue", cex=.1) +
  theme_void()

# subset points to gsp area
dcv <- ds[gsps, ] 

plot(gsp)
points(dcv, cex=.2)

xy_dcv <- coordinates(dcv)
xy_dcv <- unique(xy_dcv)
dim(xy_dcv)

area_cv <- area(gsp)
density_dcv <- nrow(xy_dcv)/area_cv
cat("Point density =", density_dcv)

dist_dcv <- dist(xy_dcv)
distmatrix_dcv <- as.matrix(dist_dcv)
diag(distmatrix_dcv) <- NA

dmin <- apply(distmatrix_dcv, 1, min, na.rm=TRUE)
head(dmin)
wdmin <- apply(distmatrix_dcv, 1, which.min)

ord <- rev(order(dmin))

far100 <- ord[1:150]

ord_new <- ord[-1:-100]
xy_new <- xy_dcv[ord_new, ]
xydf <- data.frame(xy_new)
spxydf <- SpatialPointsDataFrame(coords = xydf, data = xydf, proj4string = merc)
neighbors <- wdmin[far100]

plot(gsp) 
points(dcv, cex=0.1)
points(xy_dcv[far100, ], col='blue', pch=20)
points(xy_dcv[neighbors, ], col='red')

# Add 6 mile buffers around the points. Townships are 36 square mile grids, with a maximum uncertainty in position of $\sqrt{18}$ miles.  
# convert miles to meters
miles_to_meters <- function(x){x * 1609.34}
# township level accuracy (36 square miles)
bt <- miles_to_meters(3.5)
# calculate buffers: write a fucntion and apply it
buffer_intersect <- function(x) {
  # calculate the buffer
  x = gBuffer(x, width = bt, capStyle = "ROUND", joinStyle = "BEVEL")
  # trim the buffer to the cv
  x = raster::intersect(x, gsp)
  return(x) # return result
}

blist <- buffer_intersect(dcv)

ca <- shapefile("Boundaries/ca/CA_State_TIGER2016.shp")
ca <- spTransform(ca, merc) 

#plot(ca)
plot(gsp, main = "Coverage of Minimum Threshold Network", col= "#be9897", border = "black")
plot(blist, add = T, col = "grey90", lwd = 1)

uall <- raster::union(blist)

all_d <- dcv 
# remove duplicate locations
ad <- remove.duplicates(all_d)
ad_distinct_links <- ad@data %>% distinct(Well_ID) %>% pull(Well_ID) # distinct links
ad <- ad[match(ad_distinct_links, ad@data$Well_ID), ] # distinct wells in terms of link

bl <- gBuffer(gsp, width = bt, capStyle = "ROUND", joinStyle = "BEVEL") # township level buffer list
blc <- raster::intersect(bl, blist) # buffer list cropped to cv

# calulate buffers
buff_ts <- gBuffer(dcv, width = bt, capStyle = "ROUND", joinStyle = "BEVEL")# township level buffer
#buff_ts <- raster::intersect(buff_ts, gsp)     # trim excess buffer to cv
pc <- gArea(buff_ts) / gArea(gsp)      # percent coverage

buff_ts$Area_sqmi <- gArea(buff_ts, byid=TRUE)
gsp$Area_sqmi <- gArea(gsp, byid=TRUE)

newdf <- buff_ts$Area_sqmi / gsp$Area_sqmi


buff_ts <- spTransform(buff_ts, merc)

#write_rds(buff_ts, "Output/interpolation_buffer_ncvc.rds")

# visualize and export
plot(ca, col="white")
#plot(gsp, add = T, col="#FCFCE4", cex.main = 1.5, cex.sub = 1.5,
   #  main = "Areal Coverage", 
    # sub = paste0("Coverage = ", round(pc*100, 2), "%", "\n", 
     #             "(nwell = ", formatC(nrow(ad), big.mark = ","),")"))
plot(gsp, col="grey90")
plot(buff_ts, col = "#BDD9FA", add = T) # the buffer
#plot(uni, col = "grey90", add = T) # the buffer
#plot(ad, add=T, pch = 19, cex = .1, col = "cornflowerblue")
plot(dcv, add=T, pch = 19, cex = .1, col = "#073B78") # points that made the buffer
#legend(-13400000,  4820000, legend = c("Well", "Covered", "Uncovered"), pch = 19, col = c("yellow", "grey90", "orange"), border = FALSE, cex = 1.5)

bt <- as(buff_ts, "sf")
a <- as(ad, "sf")

ggplot()+
  geom_sf(data=gsp) +
  #scale_fill_manual("Region", values = c("San Joaquin Valley" = "#70AB53", "Coastal" = "#535EAB", "Central" = "#AB6353")) +
  geom_sf(data=bt, col="grey90") +
  geom_sf(data=a, col="cornflowerblue", cex=.1) +
  theme_void(base_size = 16)
