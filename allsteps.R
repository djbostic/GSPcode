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

#### COLORS ####
logoblue <- "#007dc5"
darkblue <- "#062e4e"
green <- "#67862f"
rust <- "#de761c"
earth <- "#9b4b08"
tan <- "#ebc668"

# load data
ca <- st_read("Boundaries/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>% filter(STUSPS == "CA")
mts <- st_read("MTs/sjvMTs.shp")
mt_sp <- as_Spatial(mts)

mtjoindata <- mts %>% st_drop_geometry(.) %>% dplyr::select(GSP_ID, GSP_Nam) %>% unique(.)
gsps <- st_read("Boundaries/CentralValleyGSPs.shp")
gsps <- left_join(gsps, mtjoindata, by=c("GSP_ID"))
gsps_sp <- as_Spatial(gsps)
domesticwells <- st_read("DomesticWells/dw_gsp.shp") 
dw <- st_intersection(domesticwells, gsps)
buff <- st_intersection(dw, buffer)

dw_screens <- buff %>% filter(year > 1990 & year < 2017 &
                       BOTTOMO < TOTALCO & 
                       TOPOFPE < TOTALCO &
                       TOPOFPE < BOTTOMO) %>%
  mutate(BottomMinusTwenty = BOTTOMO-20)

dw_tcdonly <- buff %>% filter(year > 1990 & year < 2017 &
                              is.na(TOTALCO)==FALSE)


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

buff_data <- buff %>% group_by(GSP_Nam, region) %>% summarise(Buffer=length(unique(WCRNUMB)), TCDb = mean(TOTALCO, na.rm=TRUE))
buff_data <- st_set_geometry(buff_data, NULL)

notinboth <- dw[!(dw$WCRNUMB %in% buff$WCRNUMB),]

bi <- left_join(inti_data, buff_data, by="GSP_Nam")
bi$percentcoverage <- bi$Buffer / bi$Int
bif <- filter(bi, percentcoverage < .60)
bi$pc <- ifelse(bi$percentcoverage <=0.6, "< 60", 
                ifelse(bi$percentcoverage > 0.6 & bi$percentcoverage <= 0.8, "60 - 80", 
                       ifelse(bi$percentcoverage > 0.8 & bi$percentcoverage <= 0.9, "80 - 90", 
                              ifelse(bi$percentcoverage >= 0.9, "90 - 100", "NA"))))
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
cgwl_at_dw <- raster::extract(cgwl_raster, dw_tcdonly) # intersect to get value of current water level at well points
cad <- cbind(dw_tcdonly, cgwl_at_dw) 
cad <- cad[!is.na(cad$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
cad$tcd_dry <- ifelse(cad$TOTALCO <= cad$cgwl_at_dw, "Failing", "Active")
activedw <- cad[cad$tcd_dry == "Active", ]
activedw <- activedw[!is.na(activedw$WCRNUMB), ]

length(unique(activedw$WCRNUMB)) / length(unique(cad$WCRNUMB)) # percent of wells whose TCD is below current groundwater levels (useable wells)

#### dry well analysis ####
mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
mad <- cbind(activedw, mt_at_dw) 
mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)

mad$dry <- ifelse(mad$TOTALCO <= mad$mt_at_dw, "Failing", "Active")
mad$bos_dry <- ifelse(mad$TOPOFPE <= mad$mt_at_dw & mad$BottomMinusTwenty > mad$mt_at_dw, "Partially Dewatered",
                      ifelse(mad$BottomMinusTwenty <= mad$mt_at_dw, "Fully Dewatered", "Active"))


bosdry <- mad[mad$dry == "Failing", ]
tosdry <- mad[mad$bos_dry == "Partially Dewatered", ]
active <- mad[mad$dry == "Active",]

dry <- mad[mad$bos_dry == "Fully Dewatered" | mad$bos_dry == "Partially Dewatered",]

length(unique(bosdry$WCRNUMB)) / length(unique(mad$WCRNUMB)) # percent of wells that go dry in whole data set based on TCD
length(unique(tosdry$WCRNUMB)) / length(unique(mad$WCRNUMB)) # percent of wells that go dry in whole data set based on screened interval
length(unique(dry$WCRNUMB)) / length(unique(mad$WCRNUMB))

mad <- mad %>% dplyr::select(WCRNUMB, TOTALCO, TOPOFPE, BOTTOMO, year, BASIN, GSP_Nam, cgwl_at_dw, mt_at_dw, dry, tcd_dry)
export <- st_drop_geometry(mad) %>% write.csv(., "drywells_tcdonly_1990.csv")

ggplot()+
  geom_sf(data=ca)+
  geom_sf(data=gsps)+
  geom_sf(data=mad, aes(col=dry))+
  scale_color_manual(values=c(earth, tan))





ggplot(data=mad, aes(TOTALCO, fill=factor(dry)))+
  geom_histogram(lwd=1, binwidth = 50)+  
  scale_fill_manual(name= "Well Status", values=c("cornflowerblue", logoblue, darkblue))+
  xlim(0,1000)+
  xlab("Total Completed Depth")+
  ylab("Number of Wells")+
  theme_bw()

ggplot()+
  stat_bin(data=active, aes(TOTALCO, color="1. Active"),geom="step", lwd=1, binwidth = 50)+
  stat_bin(data=bosdry, aes(TOTALCO, color="2. Dewatered"),geom="step", lwd=1, binwidth = 50)+
  #stat_bin(data=tosdry, aes(TOPOFPE, color="2. Partially Dewatered"),geom="step", lwd=1, binwidth = 50)+
  #stat_bin(data=bosdry, aes(TOPOFPE, color="3. Fully Dewatered"),geom="step", lwd=1, binwidth = 50)+
  scale_color_manual(name= "Well Status", values=c(green, earth))+
  xlab("\nTotal Completed Depth (ft bgs)")+
  ylab("Number of Wells\n")+
  theme_minimal()  + 
  coord_cartesian(xlim = c(0, 1000)) +
  theme(legend.position = c(0.75, 0.5))+
  theme(axis.title.x = element_text(color="grey20", size=14),
        axis.title.y = element_text(color="grey20", size=14),
        axis.text = element_text(color="grey20", size=12),
        legend.text = element_text(color="grey20", size=12),
        legend.title =  element_text(color="grey20", size=14))


gddw <- mad %>% group_by(GSP_Nam) %>% summarise(AllWells = length(unique(WCRNUMB)),
                                                Active = length(unique(WCRNUMB[dry=="Active"])),
                                                Dry = length(unique(WCRNUMB[dry=="Failing"])))
gddw$percfail <- 100*(gddw$Dry/gddw$AllWells)
gddw$percactive <- 100*(gddw$Active / gddw$AllWells)

bplot <- as.data.frame(st_drop_geometry(gddw)) %>% 
  filter(percfail >0) %>%
  dplyr::select(., GSP_Nam, Active=percactive, Dewatered=percfail) %>%
  gather(WellStatus, value, -c(GSP_Nam))

bplot$WellStatus_f <- factor(bplot$WellStatus , levels=c("Dewatered","Active"))

gddw <- filter(gddw, percfail>0)
ggplot(gddw, aes( y=percfail, x=reorder(GSP_Nam,percfail))) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  theme_minimal()+
  labs(x ="Groundwater Sustainability Plan", y = "Percent of Wells")+
  theme(axis.title.x = element_text(color="grey20", size=14),
        axis.title.y = element_text(color="grey20", size=14),
        axis.text = element_text(color="grey20", size=12),
        legend.text = element_text(color="grey20", size=12),
        legend.title = element_text(color="grey20", size=12))+
  theme(legend.position = c(0.75, 0.5))

#### raster map ####
# rasterize dry wells
r <- raster(buffer_sp)
res(r) <- 1610*6 # township-level grids (1,609 meters in a mile)
r <- rasterize(buffer_sp, r)
plot(r, col="grey90", lwd=10)
quads <- as(r, "SpatialPolygons")

test_spdf <- as(r, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


tcddry <- bosdry[!is.na(bosdry$WCRNUMB),]
pld <- as_Spatial(tcddry)
pld <- rasterize(coordinates(pld), r, fun="count")
pld <- as(pld, "SpatialPolygonsDataFrame")
pld <- st_as_sf(pld)
pld$dw_pl <- ifelse(pld$layer <=10, "1 - 10", 
                    ifelse(pld$layer > 10 & pld$layer <= 100, "10 - 100", 
                           ifelse(pld$layer > 100 & pld$layer <= 500, "100 - 500", 
                                  ifelse(pld$layer > 500, "500 - 1500", "NA"))))


ints <- st_intersection(pld, st_make_valid(gsps))
ints$dw_pl <- factor(ints$dw_pl, levels=c("1 - 10", "10 - 100", "100 - 500", "500 - 1500"))
ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y), fill="grey80", col="grey60")+
  geom_sf(data=ints, aes(fill=dw_pl), col=NA)  +
  geom_sf(data=gsps, fill=NA, col="grey50", lwd=.1) +
  #geom_sf(data=active, cex=.1, col="grey80") +
  scale_fill_brewer(palette = "OrRd") + 
  theme_void(base_size = 16)+ 
  labs(fill = "Number of dry wells\nper 36 square miles\n(township grids)")
