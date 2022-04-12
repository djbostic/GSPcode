# load packages
library(tidyverse) # general purpose data science toolkit
library(sp)        # spatial objects
library(raster)    # for raster objects
library(here)
library(sf)

here()
# set working directory
setwd(here())

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
gsps <- st_read(here("Boundaries/GSP_submitted/GSA_Master.shp"))
gsps <- st_transform(gsps, merc)

# load MTs
mn <- read.csv("2020_2022_MTs.csv") %>% filter(is.na(Long)==FALSE & is.na(Lat)==FALSE)
mts <- st_as_sf(mn, coords = c("Long","Lat"), crs=4326) %>% st_transform(., crs=merc)

# plot to make sure everything looks okay
ggplot()+
  geom_sf(data=gsps, col="blue")+
  geom_sf(data=mts, col="red")

mt <- as_Spatial(mts)
