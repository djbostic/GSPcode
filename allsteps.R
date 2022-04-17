# load packages
library(tidyverse) # general purpose data science toolkit
library(sp)        # spatial objects
library(raster)    # for raster objects
library(here)
library(sf)
library(automap)
library(lubridate)
library(gstat)
library(RColorBrewer)


# load data
ca <- st_read("Boundaries/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>% filter(STUSPS == "CA")
mts <- st_read("MTs/sjvMTs.shp")
mt_sp <- as_Spatial(mts)

mtjoindata <- mts %>% st_drop_geometry(.) %>% dplyr::select(GSP_ID, GSP_Nam) %>% unique(.)
gsps <- st_read("Boundaries/CentralValleyGSPs.shp")
gsps <- left_join(gsps, mtjoindata, by=c("GSP_ID"))
gsps_sp <- as_Spatial(gsps)
domesticwells <- st_read("DomesticWells/dw_gsp.shp") 
domesticwells <- domesticwells %>% filter(year > 1975 & year < 2017)
dw <- st_intersection(domesticwells, gsps)

#### Create Buffer ####
# Add 6 mile buffers around the points. Townships are 36 square mile grids, with a maximum uncertainty in position of $\sqrt{18}$ miles.  
#Therefore we want 3.5*1609 = 5632.69m for width
bufferdistance = 5632.69
buffer <- st_buffer(mts, dist=bufferdistance, joinStyle = "BEVEL", endCapStyle = "ROUND")
buffer <- st_union(buffer)

# calculate percent coverage
st_area(buffer) / st_area(gsps)

# plot
ggplot()+
  geom_sf(data=ca, fill="grey90")+
  geom_sf(data=gsps, fill="#468189", lwd=0)+
  geom_sf(data=buffer, fill="#0A2463", lwd=0)+
  theme_void()+
  labs(title="Coverage of the Minimum Threshold Representative \nMonitoring Network")+
  theme(plot.title = element_text(size=16, hjust=0.5))


#### Interpolate #### 
buffer_sp <- as_Spatial(buffer)
r <- raster(buffer_sp)           # create a template raster to interpolate over
res(r) <- 5000            # > township resolution: 6 miles = 9656.06 meters
g <- as(r, "SpatialGrid") # convert raster to spatial grid object

# subset pts to the central valley polygon
mt_sp$MT_dtw <- as.numeric(mt_sp$MT_dtw) 
mt_sp@data <- filter(mt_sp@data, MT_dtw > 0)

subset_gsp_outline <- function(x){x[buffer_sp, ]}
mt_gsp_outline <- subset_gsp_outline(mt_sp) 

# get sets of overlapping points
get_set <- function(x, y){zerodist(x)[, y]}
s1_mt <- get_set(mt_gsp_outline, 1)      # index of set 1: wells wtih an overlapping observation
s2_mt <- get_set(mt_gsp_outline, 2)      # index of set 2: wells wtih an overlapping observation

# get parallel minima of overlapping points
min_list_mt = pmin(mt_gsp_outline[s1_mt,]$MT_dtw, mt_gsp_outline[s2_mt,]$MT_dtw)

# replace DGBS of set 2 wells wtih average of set 1 and 2
mt_gsp_outline[s2_mt, "MT_dtw"] <- min_list_mt

# remove set 1 wells
mt_gsp_outline <- mt_gsp_outline[-s1_mt, ]

# fix incorrect values: remove NAs
no_na <- function(x){x[!is.na(mt_gsp_outline$MT_dtw),]}
mt_gsp_outline <- no_na(mt_gsp_outline)

# log transform Depth Below Ground Surface 
mt_gsp_outline@data$MT_dtw <- log(mt_gsp_outline@data$MT_dtw)

# plot to ensure all is working
title <- "Township Coverage \nMinimum Threshold Wells"
st <- formatC(nrow(mt_gsp_outline), big.mark = ",")
plot(buffer_sp, main = title, sub = paste0("Monitoring Wells Used: ", st))
plot(mt_gsp_outline, add = T, pch = 16, cex = .2, col = "blue")


gs_mt <- gstat(formula = MT_dtw ~ 1, # spatial data, so fitting xy as idp vars
               locations = mt_gsp_outline)        # groundwater monitoring well points 

v_mt <- variogram(gs_mt,              # gstat object
                  width = 1000)    # lag distance

plot(v_mt)
plot(autofitVariogram(MT_dtw~1, mt_gsp_outline, "Exp"))

fve_mt <- fit.variogram(v_mt,         # takes `gstatVariogram` object
                        vgm(.84,   # partial sill: semivariance at the range
                            "Exp",     # linear model type
                            38992,    # range: distance where model first flattens out
                            0.34))      # nugget

#fve_mt <- autofitVariogram(MT_dtw~1, mt_gsp_outline, "Exp")
# plot variogram and fit
plot(v_mt, fve_mt, main="Minimum Threshold Variogram")

# ordinary kriging 
kp_mt <- krige(MT_dtw ~ 1, mt_gsp_outline, g, model = fve_mt)

# backtransformed
bt_mt <- exp( kp_mt@data$var1.pred + (kp_mt@data$var1.var / 2) )

# means of backtransformed values and the sampled values
mu_bt_mt <- mean(bt_mt)
mu_original_mt <- mean(mean(exp(mt_gsp_outline$MT_dtw)))
(mu_original_mt/mu_bt_mt)

# these means differ by > 5%, thus we make another correction
#btt_mt <- bt_mt * (mu_original_mt/mu_bt_mt)
kp_mt@data$var1.pred <- bt_mt                    # overwrite w/ correct vals 
kp_mt@data$var1.var  <- exp(kp_mt@data$var1.var)  # exponentiate the variance

# covert to raster brick and crop to buff_ts
ok_mt <- brick(kp_mt)                          # spatialgrid df -> raster brick obj.
ok_mt <- mask(ok_mt, buffer_sp)                       # mask to gsp_outline extent
names(ok_mt) <- c('Prediction', 'Variance') # name the raster layers in brick

ok_mt$ci_upper <- ok_mt$Prediction + (1.96 * sqrt(ok_mt$Variance))
ok_mt$ci_lower <- ok_mt$Prediction - (1.96 * sqrt(ok_mt$Variance))

plot(ok_mt$Prediction)

#### Calculate Coverage of DW ####
# current gwl raster
cgwl_raster <- read_rds("InterpolationGWLevels/cgwl_raster.rds")

# minimum threshold raster
mt_raster <- ok_mt$Prediction

# compare how many wells are within GSP area and interpolation boundary
intial_data <- dw %>% group_by(GSP_Nam, region) %>% summarise(Int=length(unique(WCRNUMB)), TCDa = mean(TOTALCO, na.rm=TRUE))
inti_data <- st_set_geometry(intial_data, NULL)

buff <- st_intersection(dw, buffer)
buff_data <- buff %>% group_by(GSP_Nam, region) %>% summarise(Buffer=length(unique(WCRNUMB)), TCDb = mean(TOTALCO, na.rm=TRUE))
buff_data <- st_set_geometry(buff_data, NULL)

notinboth <- dw[!(dw$WCRNUMB %in% buff$WCRNUMB),]

bi <- left_join(inti_data, buff_data, by="GSP_Nam")
bi$percentcoverage <- bi$Buffer / bi$Int
bif <- filter(bi, percentcoverage < .60)
bi$pc <- ifelse(bi$percentcoverage <=0.6, "< 0.6", 
                ifelse(bi$percentcoverage > 0.6 & bi$percentcoverage <= 0.8, "0.6 - 0.8", 
                       ifelse(bi$percentcoverage > 0.8 & bi$percentcoverage <= 0.9, "0.8 - 0.9", 
                              ifelse(bi$percentcoverage >= 0.9, "0.9 - 1", "NA"))))
percentmap <- left_join(gsps, bi, by="GSP_Nam")
bi <- na.omit(bi)

ggplot() + 
  geom_sf(data=gsps, col='transparent') +
  #geom_sf(data=dw, cex=.1, col="cornflowerblue") +
  geom_sf(data=percentmap, aes(fill = pc)) +
  theme_void(base_size = 16) +  
  scale_fill_manual(values = c("#F1E5E8", "#D4AFB9","#B87A8B", "#9F5668"), name= "Percent Covered") + 
  ggtitle("Percent of Domestic Wells Covered by GSP MTs")

#### Remove wells above current groundwater level ####
cgwl_at_dw <- raster::extract(cgwl_raster, buff) # intersect to get value of current water level at well points
cad <- cbind(dw, cgwl_at_dw) 
cad <- cad[!is.na(cad$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
cad$tcd_dry <- ifelse(cad$TOTALCO <= cad$cgwl_at_dw, "Failing", "Active")
activedw <- cad[cad$tcd_dry == "Active", ]
activedw <- activedw[!is.na(activedw$WCRNUMB), ]

length(unique(activedw$WCRNUMB)) / length(unique(cad$WCRNUMB)) # percent of wells whose TCD is below current groundwater levels (useable wells)


#### assuming we dont have cgwl for coastal and non-cv central basins ####
activedw <- buff
#### dry well analysis ####
mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
mad <- cbind(activedw, mt_at_dw) 
mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
mad$tcd_dry <- ifelse(mad$TOPOFPE <= mad$mt_at_dw, "Failing", "Active")
tcddry <- mad[mad$tcd_dry == "Failing", ]
tcddry <- tcddry[!is.na(tcddry$WCRNUMB),]

length(unique(tcddry$WCRNUMB)) / length(unique(mad$WCRNUMB)) # percent of wells that go dry in whole data set based on tcd
