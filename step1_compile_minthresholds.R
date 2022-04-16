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

# all gsp boundaries
gsps <- st_read("Boundaries/GSP_submitted/GSA_Master.shp") %>% st_transform(., merc) %>% st_make_valid(.)

# add dictionary of coastal, CV, and southern basins
gsps$region <- ifelse(grepl("5-",gsps$BASIN), "Central Valley",
                      ifelse(grepl("4-", gsps$BASIN), "Ventura",
                             ifelse(grepl("6-", gsps$BASIN), "Indian Wells",
                                    ifelse(grepl("3-", gsps$BASIN), "Cuyama",
                                           ifelse(grepl("2-", gsps$BASIN), "Napa", "Other")))))

ca <- st_read("Boundaries/ca/CA_State_TIGER2016.shp")
ca <- st_transform(ca, merc)
ggplot()+
  geom_sf(data=ca, fill = "transparent") +
  geom_sf(data=gsps, aes(fill = gsps$region)) +
  scale_fill_manual("Region", values = c("Central Valley" = "#70AB53", "Ventura" = "#535EAB", "Cuyama" = "#AB6353", "Indian Wells" = "red", "Napa"="blue", "Other"="gray90")) +
  theme_void(base_size = 16)


sjvgsp <- filter(gsps, region == "Central Valley")
coastalgsp <- filter(gsps, region== "Ventura")
indianwellsgsp <- filter(gsps, region == "Indian Wells")
cuyamagsp <- filter(gsps, region=="Cuyama")

gsp <- as_Spatial(gsps)
