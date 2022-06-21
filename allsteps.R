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
red3 <- "#9b2226"
red2 <- "#ae2012"
red1 <- "#bb3e03"
orange3 <- "#ca6702"
orange2 <- "#ee9b00"
orange1 <- "#e9d8a6"
blue1 <- "#94d2bd"
blue2 <- "#0a9396"
blue3 <- "#005f73"
black <- "#001219"


#a <- st_read("MTs/MinimumThresholdWells/MinimumThresholdWells.shp")

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
ca <- st_read("Boundaries/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>% filter(STUSPS == "CA") %>% st_transform(., crs=merc)
allgwbasins <- st_read("Boundaries/i08_B118_CA_GroundwaterBasins/i08_B118_CA_GroundwaterBasins.shp") %>% st_transform(., crs=merc)
sgmabasins <- st_read("Boundaries/GSP_submitted/GSA_Master.shp") %>% st_transform(., crs=merc)

gsps <- filter(sgmabasins, grepl("5-", BASIN))
#gsps$GSP.NAME <- ifelse(is.na(gsps$Local.ID)==TRUE, gsps$BASIN, gsps$Local.ID)
gsps$GSP.ID <- as.numeric(gsps$GSP.ID)
gspcoda <- read_csv("Boundaries/gsp_coda.csv")
gsps <- left_join(gsps, gspcoda, by="GSP.ID")
gsps <- filter(gsps, is.na(GSP.NAME)==FALSE)

dacs <- st_read("Boundaries/census_data_disadvantaged_communities_2018/DAC_Pl18.shp") %>% filter(DAC18 == "Y") %>% st_transform(., merc)
dcs <- st_intersection(dacs, gsps)

ggplot()+
  geom_sf(data=ca, col=black, fill=NA)+
  geom_sf(data=allgwbasins, col=NA, fill=blue1)+
  geom_sf(data=sgmabasins, col=orange2, fill=orange1)+
  geom_sf(data=gsps, col=red1, fill=orange3)+
  theme_void()


# create shapefile from csv of MTs
# load MTs
mn <- read.csv("MTs/CentralValleyMTs.csv") %>% filter(is.na(Long)==FALSE & is.na(Lat)==FALSE & is.na(MT_dtw)==FALSE)
mts <- st_as_sf(mn, coords = c("Long","Lat"), crs=4326) %>% st_transform(., crs=merc)

# plot to make sure everything looks okay
ggplot()+
  geom_sf(data=gsps, col="red")+
  geom_sf(data=dcs, col=blue1)
  geom_sf(data=mts, col="blue")

mt_sp <- as_Spatial(mts)

mts <- st_intersection(mts, gsps)
mtjoindata <- mts %>% st_drop_geometry(.) %>% dplyr::select(GSP_ID, GSP_Nam) %>% unique(.)


gsps <- left_join(gsps, mtjoindata, by=c("GSP_ID"))
gsps_sp <- as_Spatial(gsps)

domesticwells <- st_read("DomesticWells/dw_gsp.shp")  %>% st_transform(., st_crs(gsps))
dw1 <- st_intersection(domesticwells, gsps)
dw <- dw1 %>% filter(year >= 1990 & year <= 2017 &
                      is.na(TOTALCO)==FALSE)
hi <- dw %>% group_by(GSP_Nam) %>% summarise(`Number of Domestic Wells` = length(unique(WCRNUMB)),`Average TCD` = mean(TOTALCO, na.rm=TRUE), length(unique(WCRNUMB)))


richsdws <- read_rds("DomesticWells/domcv6_mean_gw_with_beta_GF_CI.rds") %>% st_as_sf(.) %>% st_transform(., st_crs(gsps)) %>% filter(year >=1990)
dw <- st_intersection(richsdws, gsps)
hi <- dw %>% group_by(GSP.NAME) %>% summarise(`Number of Domestic Wells` = length(unique(WCRNumber)),`Average TCD` = mean(TotalCompletedDepth, na.rm=TRUE), `Average Pump Depth` = mean(pump_loc), `Fraction of All DWs`=100*length(unique(WCRNumber))/nrow(dw)) %>% st_drop_geometry(.)

write_csv(hi, "DomesticWells/dwbygsp.csv")

#### Create Buffer ####
# Add 6 mile buffers around the points. Townships are 36 square mile grids, with a maximum uncertainty in position of $\sqrt{18}$ miles.  
#Therefore we want 3.5*1609 = 5632.69m for width
bufferdistance = 5632.69
buffer <- st_buffer(mts, dist=bufferdistance, joinStyle = "BEVEL", endCapStyle = "ROUND")
buffer <- st_union(buffer)

# calculate percent coverage
sum(st_area(buffer)) / sum(st_area(gsps))

# plot
ggplot()+
  #geom_sf(data=ca, fill="grey90")+
  geom_sf(data=gsps, fill=blue3, col=blue2, lwd=.1)+
  geom_sf(data=buffer, fill=red1, col=orange1, lwd=.5)+
  theme_void()+
  #labs(title="Coverage of the Minimum Threshold Representative \nMonitoring Network")+
  theme(plot.title = element_text(size=16, hjust=0.5))

buff <- st_intersection(dw, buffer)
dac_winbuff <- st_intersection(dcs, buffer)

notcovered <- dcs[!(dcs$PLACEFP %in% dac_winbuff$PLACEFP),]

ggplot()+
  geom_sf(data=gsps, fill=blue2, col=blue3, lwd=.1)+
  geom_sf(data=buffer, fill=orange1, col=orange3, lwd=.5)+
  geom_sf(data=dcs, fill=red3, col=red2)+
  geom_sf(data=dac_winbuff, fill=black, col=black, lwd=.5)+
  theme_void()

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
plot(autofitVariogram(MT_dtw~1, mt_gsp_outline, "Cir"))

fve_mt <- fit.variogram(v_mt,         # takes `gstatVariogram` object
                        vgm(.5,   # partial sill: semivariance at the range
                            "Exp",     # linear model type
                            450000,    # range: distance where model first flattens out
                            0.2))      # nugget

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

# these means differ by < 5%, thus we make another correction
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

write_rds(ok_mt$Prediction, "InterpolationGWLevels/mt_raster.rds")

#### Calculate Coverage of DW ####
# current gwl raster
cgwl_raster <- read_rds("InterpolationGWLevels/cgwl_raster.rds")

# minimum threshold raster
mt_raster <- ok_mt$Prediction

# compare how many wells are within GSP area and interpolation boundary
intial_data <- dw %>% group_by(GSP.NAME) %>% summarise(Int=length(unique(WCRNumber)), TCDa = mean(TotalCompletedDepth, na.rm=TRUE), PL = mean(pump_loc, na.rm=TRUE))
inti_data <- st_set_geometry(intial_data, NULL)

buff_data <- buff %>% group_by(GSP.NAME) %>% summarise(Buffer=length(unique(WCRNumber)), TCDb = mean(TotalCompletedDepth, na.rm=TRUE), PLb = mean(pump_loc, na.rm=TRUE))
buff_data <- st_set_geometry(buff_data, NULL)

notinboth <- dw[!(dw$WCRNumber %in% buff$WCRNumber),]

bi <- left_join(inti_data, buff_data, by="GSP.NAME")
bi$percentcoverage <- bi$Buffer / bi$Int
bif <- filter(bi, percentcoverage < .60 & is.na(percentcoverage)==FALSE)
bi$pc <- ifelse(bi$percentcoverage <=0.6, "< 60", 
                ifelse(bi$percentcoverage > 0.6 & bi$percentcoverage <= 0.8, "60 - 80", 
                       ifelse(bi$percentcoverage > 0.8 & bi$percentcoverage <= 0.9, "80 - 90", 
                              ifelse(bi$percentcoverage >= 0.9, "90 - 100", "NA"))))
percentmap <- left_join(gsps, bi, by="GSP.NAME")
bi <- na.omit(bi)
percentmap <- filter(percentmap, is.na(pc)==FALSE)

ggplot() + 
  #geom_sf(data=gsps, col='transparent') +
  #geom_sf(data=dw, cex=.1, col="cornflowerblue") +
  geom_sf(data=percentmap, aes(fill = pc), col="gray20") +
  theme_void(base_size = 20) +
  #scale_fill_brewer(palette = "BlGr", name="Percent Coverage") +   
  scale_fill_manual(values = c(red3, orange3, blue2, blue3), name= "Percent Covered")

diffraster <- ok_mt$Prediction - cgwl_raster$layer
extr <- raster::extract(diffraster, mts)
abc <- cbind(mts, extr)
abc <- abc %>% group_by(GSP_Nam) %>% mutate(NumWells = length(unique(Well_ID)))

abc %>% filter(NumWells > 35) %>% 
  ggplot(., aes(x=extr, y=reorder(GSP_Nam, extr, na.rm=TRUE))) +
  geom_boxplot() + 
  theme_minimal(base_size=20) +
  labs(y="Groundwater Susatinability Plan", x="Change in Groundwater Level (ft)")

#### Remove wells above current groundwater level ####
cgwl_at_dw <- raster::extract(cgwl_raster, buff) # intersect to get value of current water level at well points
cad <- cbind(buff, cgwl_at_dw) 
cad <- cad[!is.na(cad$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
cad$tcd_dry <- ifelse(cad$TotalCompletedDepth <= cad$cgwl_at_dw, "Failing", "Active")
activedw <- cad[cad$tcd_dry == "Active", ]
activedw <- activedw[!is.na(activedw$WCRNumber), ]

length(unique(activedw$WCRNumber)) / length(unique(cad$WCRNumber)) # percent of wells whose TCD is below current groundwater levels (useable wells)

#### dry well analysis ####
mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
mad <- cbind(activedw, mt_at_dw) 
mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)

mad$dry <- ifelse(mad$TotalCompletedDepth <= mad$mt_at_dw, "Failing", "Active")
mad$bos_dry <- ifelse(mad$pump_loc <= mad$mt_at_dw & mad$TotalCompletedDepth > mad$mt_at_dw, "Partially Dewatered",
                      ifelse(mad$TotalCompletedDepth <= mad$mt_at_dw & mad$pump_loc <= mad$mt_at_dw, "Fully Dewatered", "Active"))


bosdry <- mad[mad$bos_dry == "Fully Dewatered", ]
tosdry <- mad[mad$bos_dry == "Partially Dewatered", ]
active <- mad[mad$bos_dry == "Active",]

dry <- mad[mad$bos_dry == "Fully Dewatered" | mad$bos_dry == "Partially Dewatered",]

length(unique(bosdry$WCRNumber)) / length(unique(mad$WCRNumber)) # percent of wells that go dry in whole data set based on TCD
length(unique(tosdry$WCRNumber)) / length(unique(mad$WCRNumber)) # percent of wells that go dry in whole data set based on screened interval
length(unique(dry$WCRNumber)) / length(unique(mad$WCRNumber))

mad <- mad %>% dplyr::select(WCRNumber, TotalCompletedDepth, pump_loc, top, bot, year, BASIN, GSP.NAME, cgwl_at_dw, mt_at_dw, dry, bos_dry)
export <- st_drop_geometry(mad) %>% write.csv(., "drywells_tcdonly_1990.csv")

abdf <- mad %>% group_by(GSP.NAME, bos_dry) %>% summarize(Count=length(unique(WCRNumber)))
d <- abdf %>% spread(bos_dry, Count) %>% mutate(totals = sum(Active, `Fully Dewatered`, `Partially Dewatered`, na.rm=TRUE))

d$percFully <- d$`Fully Dewatered`/d$totals
d$percPartial <- d$`Partially Dewatered`/d$totals
write.csv(st_drop_geometry(d), "DomesticWells/DRYWELLTOTALS.csv")

ggplot(data=mad, aes(TotalCompletedDepth, fill=factor(dry)))+
  geom_histogram(lwd=1, binwidth = 50)+  
  scale_fill_manual(name= "Well Status", values=c("cornflowerblue", blue2, blue3))+
  xlim(0,1000)+
  xlab("Total Completed Depth")+
  ylab("Number of Wells")+
  theme_bw()

ggplot()+
  stat_bin(data=active, aes(TotalCompletedDepth, color="1. Active"),geom="step", lwd=1, binwidth = 50)+
  stat_bin(data=bosdry, aes(TotalCompletedDepth, color="2. Fully Dewatered"),geom="step", lwd=1, binwidth = 50)+
  stat_bin(data=tosdry, aes(TotalCompletedDepth, color="3. Partially Dewatered"),geom="step", lwd=1, binwidth = 50)+
  #stat_bin(data=bosdry, aes(TOPOFPE, color="3. Fully Dewatered"),geom="step", lwd=1, binwidth = 50)+
  scale_color_manual(name= "Well Status", values=c(blue3, blue2, orange1))+
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


gddw <- mad %>% group_by(GSP.NAME) %>% summarise(AllWells = length(unique(WCRNumber)),
                                                Active = length(unique(WCRNumber[bos_dry=="Active"])),
                                                PartiallyDry = length(unique(WCRNumber[bos_dry=="Partially Dewatered"])),
                                                FullyDry = length(unique(WCRNumber[bos_dry=="Fully Dewatered"])))
gddw$percfail <- 100*((gddw$PartiallyDry+gddw$FullyDry)/gddw$AllWells)
gddw$percactive <- 100*(gddw$Active / gddw$AllWells)

bplot <- as.data.frame(st_drop_geometry(gddw)) %>% 
  filter(percfail >0, AllWells > 50) %>%
  dplyr::select(., GSP.NAME, Active=Active, `Partially Dewatered`=PartiallyDry, `Fully Dewatered`=FullyDry) %>%
  gather(WellStatus, value, -c(GSP.NAME))

bplot$WellStatus_f <- factor(bplot$WellStatus , levels=c("Fully Dewatered", "Partially Dewatered","Active"))

gddw <- filter(gddw, percfail>0)
ggplot(gddw, aes( y=percfail, x=reorder(GSP.NAME,percfail))) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  theme_minimal()+
  labs(x ="Groundwater Sustainability Plan", y = "Percent of Dewatered Domestic Wells")+
  theme(axis.title.x = element_text(color="grey20", size=14),
        axis.title.y = element_text(color="grey20", size=14),
        axis.text = element_text(color="grey20", size=12),
        legend.text = element_text(color="grey20", size=12),
        legend.title = element_text(color="grey20", size=12))+
  theme(legend.position = c(0.75, 0.5))

bbplot <- bplot %>% group_by(GSP.NAME) %>% filter(sum(value) > 500)
ggplot(bbplot, aes(y=value, x=reorder(GSP.NAME,value), fill=WellStatus_f)) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  theme_minimal()+
  labs(x ="Groundwater Sustainability Plan", y = "Number of Domestic Wells at Minimum Threshold Groundwater Level", fill="Well Status")+
  scale_fill_manual(values = c(red2, orange3, blue3))+
  theme(axis.title.x = element_text(color="grey20", size=16),
        axis.title.y = element_text(color="grey20", size=16),
        axis.text = element_text(color="grey20", size=14),
        legend.text = element_text(color="grey20", size=14),
        legend.title = element_text(color="grey20", size=14))+
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


activesquares <- as_Spatial(activedw)
activesquares <- rasterize(coordinates(activesquares), r, fun="count")
activesquares <- as(activesquares, "SpatialPolygonsDataFrame")
activesquares <- st_as_sf(activesquares)

trya <- raster::merge(activesquares, pld)

tcddry <- dry[!is.na(dry$WCRNumber),]
pld <- as_Spatial(tcddry)
pld <- rasterize(coordinates(pld), r, fun="count")
pld <- as(pld, "SpatialPolygonsDataFrame")
pld <- st_as_sf(pld)
pld$dw_pl <- ifelse(pld$layer <=50, "1 - 50", 
                    ifelse(pld$layer > 50 & pld$layer <= 100, "51 - 100", 
                           ifelse(pld$layer > 100 & pld$layer <= 300, "101 - 300", 
                                  ifelse(pld$layer > 300, "301 - 630", "NA"))))


ints <- st_intersection(pld, st_make_valid(gsps))
ints$dw_pl <- factor(ints$dw_pl, levels=c("1 - 50", "51 - 100", "101 - 300", "301 - 630"))
ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y), fill="grey80", col="grey60")+
  geom_sf(data=ints, aes(fill=dw_pl), col=NA)  +
  geom_sf(data=gsps, fill=NA, col="grey50", lwd=.1) +
  #geom_sf(data=active, cex=.1, col="grey80") +
  scale_fill_brewer(palette = "OrRd") + 
  theme_void(base_size = 16)+ 
  labs(fill = "Number of partially and \nfully dewatered wells\nper 36 square miles\n(township grids)")
