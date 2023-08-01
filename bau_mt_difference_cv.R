# load packages needed
library(tidyverse)
library(sf)
library(dplyr)
library(mapview)
library(ggplot2)
library(automap)
library(raster) 
library(rgeos)

# set working directory
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("/Volumes/GoogleDrive/My Drive/Graduate School/GSP_Analy/Organized_CnD/")

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
# GSP boundaries
gsp <- gsps
gsp <- st_transform(gsp, merc)

mt <- read_rds("InterpolationGWLevels/mt_raster.rds")
mt <- mt$Prediction # minimum threshold raster

cv.outline <- readr::read_rds("bau_rp_erl/b_trim_fill.rds") # Outline of study area
gsp.outlines <- gsps %>% st_zm(.,drop = TRUE) %>% as_Spatial(.)

cgwl <- read_rds("InterpolationGWLevels/bau15.rds") # 2019 current groundwater level raster
cgwl <- mask(cgwl, cv.outline) # crop raster to study area

bau10 <- raster("bau_rp_erl/r2040_10.tif")
bau15 <- raster("bau_rp_erl/r2040_15.tif")
bau20 <- raster("bau_rp_erl/r2040_20.tif")

ss10 <- raster("bau_rp_erl/r2020_10.tif")
ss15 <- raster("bau_rp_erl/r2020_15.tif")
ss20 <- raster("bau_rp_erl/r2020_20.tif")

# reset all cells less than 0 to NA
bau10[bau10 < 0] <- NA
bau15[bau15 < 0] <- NA
bau20[bau20 < 0] <- NA

ss10[ss10 < 0] <- NA
ss15[ss15 < 0] <- NA
ss20[ss20 < 0] <- NA

origin(cgwl) <- origin(bau10) # reset the origin of the gwl so they match
origin(mt) <- origin(bau10)
#domwells <- read_rds("Output/domesticwells_allcobs_90.rds")

#dw_gsp <- st_intersection(dw_gsp, gsp.outlines)
#drywells1 <- as(pldry, "sf") 
#drywells <- sf::as_Spatial(drywells1)

#drywells1 <- drywells1[!is.na(drywells1$pump_loc), ]
#dww <- subset(drywells1, year > 1990)
#### Look at properties of rasters ####

cellStats(bau10, 'max') # max = 978 ft --> should decrease be a larger number than increase?
cellStats(bau10, 'min') # min = -263 ft --> this leads me to think that BAU is showing change in gwl with decline positive and increase as negative

cellStats(ss10, 'max') # max = 602
cellStats(ss10, 'min') # min = -42

cellStats(cgwl, 'max') # max = 619 ft
cellStats(cgwl, 'min') # min = 4 ft

#### Look at differences between BAU and Strict Sustainability (SS) ####
checkBAUandSS <- function(rast1,rast2,rast3, rast4, rast5, rast6)
{ set1_rasters <- list(rast1,rast2,rast3) # set 1 is 3 BAU scenarios
set2_rasters <- list(rast4, rast5, rast6) # set 2 is 3 SS scenarios

type <- c("10", "15", "20") # names of the 3 scenarios 
merc_rast1 <- list()
merc_rast2 <- list()
diff_rast <- list()
df_rast <- list()

for (i in c(1:3)){
  merc_rast1[[i]] <- projectRaster(set1_rasters[[i]], 
                                   crs=merc) # reproject both rasters
  merc_rast2[[i]] <- projectRaster(set2_rasters[[i]], 
                                   crs=merc)
  
  diff_rast[[i]] <- merc_rast1[[i]] - merc_rast2[[i]] # difference the rasters
  names(diff_rast[[i]]) <- "z" 
  
  df_rast[[i]] <- as.data.frame(diff_rast[[i]], xy = TRUE) %>% 
    dplyr::mutate(type=type[i]) %>% 
    dplyr::filter(!is.na(z))
}
rast_df <- dplyr::bind_rows(df_rast[[1]], df_rast[[2]], df_rast[[3]]) # combine 3 scenarios into one data frame
return(rast_df)
}

l1 <- checkBAUandSS(bau10, bau15, bau20, ss10, ss15, ss20)

# this graph is showing the areas where BAU is LESS than SS (where changes in the groundwater in BAU are smaller than changes in SS).. not sure how informative this is
ggplot(mapping = aes(x = x, y = y, fill = z)) + 
  geom_raster(data = l1, na.rm = TRUE) +
  geom_polygon(data=cv.outline, aes(x=long, y=lat, group=group),
               col='black', fill="NA") + 
  scale_fill_distiller("BAU - SS (ft)", palette = "Spectral", limits=c(-300,0))+ 
  facet_wrap(~type) +
  theme_void() +
  coord_fixed(1.3)

#### Compare MT to BAU/SS Scenarios #### 
gw_levels <- function(rast1,rast2,rast3, rast4, rast5)
{ set1_rasters <- list(rast1, rast2, rast3) # set 1 is 3 scenarios
merc_rast2 <- projectRaster(rast4, 
                            crs=merc) # change crs of cgwl raster
merc_rast3 <- projectRaster(rast5, 
                            crs=merc) # change crs of MT raster

type <- c("10", "15", "20") # names of the 3 scenarios 
merc_rast1 <- list()
scenariosurface_rast <- list()
diff_rast <- list()
df_rast <- list()

# create theoretical groundwater surface by adding the scenario to the 2017 groundwater level
for (i in c(1:3)){
  merc_rast1[[i]] <- projectRaster(set1_rasters[[i]], 
                                   crs=merc) # reproject scenario rasters
  
  scenariosurface_rast[[i]] <- merc_rast1[[i]] # difference the rasters merc_rast2 +
  names(scenariosurface_rast[[i]]) <- "z" 
  
  scenariosurface_rast[[i]][scenariosurface_rast[[i]] < 0] <- 0 # reassign negative values to 0 (prevent 'ponding')
}

# Difference the hypothetical scenario surface and the MT surface to see how close they are to each other
for (i in c(1:3)){
  diff_rast[[i]] <- scenariosurface_rast[[i]] - merc_rast3 # difference the scenario and MT scenariosurface_rast[[i]] - 
  names(diff_rast[[i]]) <- "z" 
  
  df_rast[[i]] <- as.data.frame(diff_rast[[i]], xy = TRUE) %>% 
    dplyr::mutate(type=type[i]) %>% 
    dplyr::filter(!is.na(z))
}

rast_df <- dplyr::bind_rows(df_rast[[1]], df_rast[[2]], df_rast[[3]]) # combine 3 scenarios into one data frame
return(rast_df)
}

bau_level <- gw_levels(bau10, bau15, bau20, cgwl, mt)

# reassign values > 200 or < 200 ft to see differences more clearly
bau_level$w <- ifelse(bau_level$z > 200, 201, ifelse(bau_level$z < -200, -201, bau_level$z))

# this is showing the difference between the BAU scenario and the MT surface
# POSITIVE = BAU > MT, NEGATIVE = MT > BAU
ggplot(mapping = aes(x = x, y = y, fill = z)) + 
  geom_raster(data = bau_level, na.rm = TRUE) +
  geom_polygon(data=gsp.outlines, aes(x=long, y=lat, group=group),
               col='black', fill="NA") + 
  scale_fill_distiller("BAU Surface - MT Surface (ft) \nPositive = BAU > MT \nNegative = MT > BAU", palette = "RdBu")+ 
  facet_wrap(~type) +
  labs(title="  Difference between BAU and MT Scenarios", subtitle = "  > is deeper") + 
  theme_void() +
  coord_fixed(1.3)


ss_level <- gw_levels(ss10, ss15, ss20, cgwl, mt)
# reassign values > 200 or < 200 ft to see differences more clearly
ss_level$w <- ifelse(ss_level$z > 200, 201, ifelse(ss_level$z < -200, -201, ss_level$z))

ggplot(mapping = aes(x = x, y = y, fill = w)) + 
  geom_raster(data = ss_level, na.rm = TRUE) +
  geom_polygon(data=gsp.outlines, aes(x=long, y=lat, group=group),
               col='black', fill="NA") + 
  scale_fill_distiller("SS Surface - MT Surface (ft) \nPositive = SS > MT \nNegative = MT > SS", palette = "RdBu")+ 
  facet_wrap(~type) +
  labs(title="  Difference between SS and MT Scenarios", subtitle = "  > is deeper") + 
  theme_void() +
  coord_fixed(1.3)



#### Dry well analysis #### 
scenariosurface <- function(rast1,rast2,rast3, rast4)
{ set1_rasters <- list(rast1, rast2, rast3) # set 1 is 3 scenarios
merc_rast2 <- projectRaster(rast4, 
                            crs=merc) # change crs of cgwl raster

type <- c("10", "15", "20") # names of the 3 scenarios 
merc_rast1 <- list()
scenariosurface_rast <- list()
diff_rast <- list()
df_rast <- list()

# create theoretical groundwater surface by adding the scenario to the 2017 groundwater level
for (i in c(1:3)){
  merc_rast1[[i]] <- projectRaster(set1_rasters[[i]], 
                                   crs=merc) # reproject scenario rasters
  
  scenariosurface_rast[[i]] <- merc_rast1[[i]] # difference the rasters
  names(scenariosurface_rast[[i]]) <- "z" 
  
  scenariosurface_rast[[i]][scenariosurface_rast[[i]] < 0] <- 0 # reassign negative values to 0 (prevent 'ponding')
}
for (i in c(1:3)){
  names(scenariosurface_rast[[i]]) <- "z" 
}

return(scenariosurface_rast)
}

bau <- scenariosurface(bau10, bau15, bau20, cgwl)

drywell <- function(scenariogwl, wells){
  int <- raster::extract(scenariogwl, wells)
  i2 <- cbind(wells, int)
  #names(i2@data)[72] <- "int"
  i <- i2[!is.na(i2$int), ]
  
  
  #i$upper_dry <- ifelse(i$mean_ci_upper <= i$int, "Failing", "Active")
  #i$lower_dry <- ifelse(i$mean_ci_lower <= i$int, "Failing", "Active")
  i$pumloc_dry <- ifelse(i$pump_loc <= i$int, "Failing", "Active")
  i$tcd_dry <- ifelse(i$TotalCompletedDepth <= i$int, "Failing", "Active")
  
  
  #upperdry <- i[i$upper_dry == "Failing", ]
  #lowerdry <- i[i$lower_dry == "Failing", ]
  pldry <- i[i$pumloc_dry == "Failing", ]
  tcddry <- i[i$tcd_dry == "Failing", ]
  
  #print(length(unique(upperdry$WCRNumber)) / length(unique(i$WCRNumber)))
  #print(length(unique(lowerdry$WCRNumber)) / length(unique(i$WCRNumber)))
  dw <- c(length(unique(pldry$WCRNumber)), length(unique(pldry$WCRNumber)) / length(unique(i$WCRNumber)), length(unique(tcddry$WCRNumber)), length(unique(tcddry$WCRNumber)) / length(unique(i$WCRNumber)))
  ret <- list(dw, pldry)
  return(ret)
}



b10 <- drywell(bau[[1]], drywells1)[[1]]
b15 <- drywell(bau[[2]], drywells1)[[1]]
b20 <- drywell(bau[[3]], drywells1)[[1]]
mtd <- drywell(mt$Prediction, drywells1)[[1]]
a <- as.data.frame(rbind(b10, b15,b20,mtd))
colnames(a) <- c("Count_PL","Perc_PL", "Count_TCD", "Perc_TCD")
rownames(a) <- c("bau10", "bau15", "bau20", "mt")

plot(gsp.outlines)
plot(dww, add=T, cex=0.1)
plot(b, add=T, col="red", cex=.1)
ggplot() + 
  geom_polygon(data=gsp.outlines, aes(x=long, y=lat, group=group),
               col='black', fill="NA") +
  geom_polygon(data=dww, aes(x=long, y=lat, col='grey20')) + 
  geom_polygon(data=b, aes(x=long, y=lat,col='red')) + 
  labs(title="Dry wells under BAU scenario") + 
  theme_void() +
  coord_fixed(1.3)

#####################
mt <- readr::read_rds("Output/minthreshinterpolation_CV.rds")
mt <- mt$Prediction # minimum threshold raster
mt <- projectRaster(mt, 
              crs=merc)

cv.outline <- readr::read_rds("Boundaries/cv/cvoutline.rds") # Outline of study area
gsp.outlines <- gsp %>% st_zm(.,drop = TRUE) %>% as_Spatial(.)
bau15 <- raster("bau_rp_erl/r2040_15.tif")
bau15 <- projectRaster(bau15, 
                    crs=merc)

plot(bau15)
plot(gsp.outlines, add=T)
plot(mt)

bau15[bau15 < 0] <- 0
bau15 <- mask(x=bau15, mask=gsp.outlines)
bau15 <- intersect(bau15, gsp.outlines)
names(bau15) <- "z"
names(mt) <- "z"

bau <- as.data.frame(bau15, xy = TRUE) %>% 
  dplyr::mutate(type="Business As Usual") %>% 
  dplyr::filter(!is.na(z))

mt <- as.data.frame(mt, xy = TRUE) %>% 
  dplyr::mutate(type="Minimum Threshold") %>% 
  dplyr::filter(!is.na(z))

rast_df <- dplyr::bind_rows(bau, mt)  
  dplyr::filter(!is.na(Prediction))# combine 3 scenarios into one data frame

rast_df$z <- ifelse(rast_df$z > 800, 801, ifelse(rast_df$z < 50, 51, rast_df$z))

baumt <- ggplot(mapping = aes(x = x, y = y, fill = z)) + 
  geom_raster(data = rast_df, na.rm = TRUE) +
  geom_polygon(data=gsp.outlines, aes(x=long, y=lat, group=group),
               col='black', fill="NA") + 
  scale_fill_distiller("Depth to Water (ft)", palette = "RdBu")+ 
  facet_wrap(~type) +
  labs(title="", subtitle = "") + 
  theme_void(base_size = 14) +
  theme(legend.position = c(0.98, 0.70))+
  coord_fixed(1.3)
#ggsave(plot=baumt, filename="Results/Images/BAU_MT_Comp.png", bg = "transparent", type = "cairo", device = "png")

###### BAU MT Fig 1 #####
#### FIGURE 1 - GW LEVELS ####
bau_at_mts <- raster::extract(bau15, mts) # intersect to get value of current water level at well points
mtz <- cbind(mts, bau_at_mts) 

mt_at_mts <- raster::extract(mt_raster, mtz) # compare interpolated MT water level with MTs set at same location
mtzz <- cbind(mtz, mt_at_mts) 

mtzz$mt_cgwl <- as.numeric(mtzz$MT_dtw) - mtzz$bau_at_mts
mtzz$decline <- mtzz$bau_at_mts-as.numeric(mtzz$MT_dtw)
mtzz$BASIN <- substring(mtzz$BASIN, 9)

# a map of the Central Valley aquifer system with each MT location as a point colored by change in groundwater level as a binned value on a diverging color scale (>50 ft, 0 to 50 ft, 0 to -50 ft, -50 to -100 ft, -100 to -150, -150 to -200, < -200).

mtzz$bin <- ifelse(mtzz$decline >= 50, "> 50 ft",
                   ifelse(mtzz$decline < 50 & mtzz$decline >= 0, "0 to 50 ft",
                          ifelse(mtzz$decline < 0 & mtzz$decline >= -50, "-50 to 0 ft",
                                 ifelse(mtzz$decline < -50 & mtzz$decline >= -100, "-100 to -50 ft",
                                        ifelse(mtzz$decline < -100 & mtzz$decline >= -150, "-150 to -100 ft",
                                               ifelse(mtzz$decline < -150 & mtzz$decline >= -200, "-200 to -150 ft",
                                                      ifelse(mtzz$decline < -200, "< -200 ft",
                                                             "error")))))))
mtzz$bins <- factor(mtzz$bin, levels=c("> 50 ft", "0 to 50 ft", "-50 to 0 ft", "-100 to -50 ft", "-150 to -100 ft","-200 to -150 ft","< -200 ft"))
mtzz <- filter(mtzz, is.na(mtzz$bins)==FALSE)

my_colors <- rev(brewer.pal(9, "RdBu")[c(1:4,7:9)])

gwl_map <- ggplot()+
  geom_sf(data=basins, fill="gray95", col=NA)+
  geom_point(data=mtzz, alpha=.85,
             aes(color = bins, geometry = geometry),
             stat = "sf_coordinates", show.legend = TRUE)+
  scale_color_manual("Change in Groundwater \nLevel (CGWL-MT)",
                     values=my_colors
  )+ 
  geom_sf(data=basins, fill=NA, col="gray40", lwd=.3)+
  theme_void()+
  theme(legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.position = c(0.98, 0.70)) +
  guides(col = guide_legend(override.aes = list(size=3.5)))
gwl_map

inset <- readPNG("Results/Images/InsetMap2.png", native=TRUE)

gwl_map2 <- gwl_map + inset_element(p = inset,
                                    left = 0,
                                    bottom = 0,
                                    right = .5,
                                    top = .3)

# Estimated groundwater elevation changes (x-axis) at subbasin (y-axis) under proposed MTs
# add median of groundwater levels
mts_basin <- mtzz %>% group_by(BASIN) %>% summarise(MedMT = median(decline)) 
mts_basin$bin <- ifelse(mts_basin$MedMT >= 50, "> 50 ft",
                        ifelse(mts_basin$MedMT < 50 & mts_basin$MedMT >= 0, "0 to 50 ft",
                               ifelse(mts_basin$MedMT < 0 & mts_basin$MedMT >= -50, "-50 to 0 ft",
                                      ifelse(mts_basin$MedMT < -50 & mts_basin$MedMT >= -100, "-100 to -50 ft",
                                             ifelse(mts_basin$MedMT < -100 & mts_basin$MedMT >= -150, "-150 to -100 ft",
                                                    ifelse(mts_basin$MedMT < -150 & mts_basin$MedMT >= -200, "-200 to -150 ft",
                                                           ifelse(mts_basin$MedMT < -200, "< -200 ft",
                                                                  "error")))))))
mts_basin$bins <- factor(mts_basin$bin, levels=c("> 50 ft", "0 to 50 ft", "-50 to 0 ft", "-100 to -50 ft", "-150 to -100 ft","-200 to -150 ft","< -200 ft"))


elchmt <- mtzz %>% group_by(BASIN)  %>% st_join(., mts_basin)

mycolors2 <- c(rev(brewer.pal(9, "RdBu")[c(1:3,7:8)]), "black", "black")
elevchange <- 
  ggplot()+
  geom_jitter(data=elchmt, aes(x=reorder(BASIN.x,-MedMT, na.rm=TRUE), y=decline, col=bins.x), fill = 4, alpha = 0.5, size = 1.5, show.legend = FALSE)+
  geom_point(data=mts_basin, aes(x=reorder(BASIN, -MedMT, na.rm=TRUE), y=MedMT, fill=bins), pch=21, alpha=1, size=3, col="black", show.legend = FALSE)+
  scale_color_manual("Change in Groundwater \nLevel (BAU-MT)",
                     values=my_colors)+ 
  scale_fill_manual("Change in Groundwater Level \n(BAU-MT)",
                    values=my_colors
  )+ 
  coord_flip(ylim=c(100,-350))+
  theme_light()+
  xlab("Groundwater Subbasin") +
  ylab("Change in Groundwater Level (BAU - MT) (ft)")+
  theme(axis.text.x =element_text(size=12), 
        axis.text.y =element_text(size=10),
        axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))+ 
  guides(col = guide_legend(override.aes = list(size=2, fill=NA)), fill = "none")
elevchange

t1 <- plot_grid(baumt, elevchange, align="hv", labels = c("a.", "b."), label_size=20)
ggsave(plot=t1, filename="Results/Images/baumt2.png", bg = "white", type = "cairo", device = "png", height = 7, width = 14.5)
