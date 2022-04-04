# load packages
library(tidyverse) # general purpose data science toolkit
library(sp)        # spatial objects
library(raster)    # for raster objects

# set working directory
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/Volumes/GoogleDrive/My Drive/Graduate School/GSP_Analy/Organized_CnD/")

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
gsps <- st_read("Boundaries/all_cobs/All_CO_GSPs.shp")
gsps <- st_transform(gsps, merc)

# csv names
monitoring_network <- list.files("1Compile_MTs/Data/MTs/", pattern = "\\.shp$")

# remove non-CV GSPs - Borrego, Cuyama, Santa Cruz, PasoRobles, Salinas, Oxnard, Pleasant Valley, Indian Wells
#monitoring_network <- monitoring_network[-c(3, 8, 14:16, 21, 36, 39, 40)]

# remove only Borrego
monitoring_network <- monitoring_network[-3]

# data list
mn <- list() 

for(i in 1:length(monitoring_network)){
  mn[[i]] = st_read(file.path("1Compile_MTs/Data/MTs/", monitoring_network[i]), stringsAsFactors = FALSE)
}

for(i in 1:length(mn)){
  mn[[i]]$MT_dtw <- as.numeric(as.character(mn[[i]]$MT_dtw))
}

for(i in 1:length(mn)){
  mn[[i]] <- subset(mn[[i]], select = c("Well_ID", "MT_dtw"))
}

# remove null geometry from Kern River
#mn[[17]] <- mn[[17]][-36,]

# add column with GSP name
for(i in 1:length(mn)){
  mn[[i]]$GSP_Name <- monitoring_network[i]
  mn[[i]]$GSP_Name <- gsub("_.*","", mn[[i]]$GSP_Name)
}

mt_merc <- list()
MTs <- mn
for(t in 1:length(MTs)){
  mt_merc[[t]] <- st_transform(MTs[[t]], merc)
}
# idk why mt_merc[[4]]'s coordinates are different, but here's a fix
#mt_merc[[4]] <- st_zm(mt_merc[[4]], drop=TRUE)

# plot to make sure everything looks okay
plot(gsps$geometry)
for(t in 1:length(mt_merc)){
  plot(mt_merc[[t]]$geometry, add=T, cex=.2, col="cornflowerblue") 
}

# combine all MTs into one shapefile
all_MTs <- do.call(rbind, mt_merc)
all_MTs <- st_zm(all_MTs, drop=TRUE)
q <- all_MTs[!st_is_empty(all_MTs),drop=FALSE]
mt <- as_Spatial(q)

# transform to mercator
ds <- spTransform(mt, merc)

# output MTs
#write_rds(ds, "Output/all_mts.rds")

dd <- as(ds, "sf")
st_write(hey, "1Compile_MTs/Output/MinimumThresholdWells.shp")

hey <- filter(hey, region=="Coastal")

ox <- filter(hey, GSP_Name.1 == "Santa Cruz Mid-County")
ggplot() +
  geom_sf(data=hey, aes(col=MT_dtw))
