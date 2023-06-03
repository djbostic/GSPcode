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

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
# ca
ca <- st_read("Boundaries/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>% filter(STUSPS == "CA") %>% st_transform(., crs=merc)

# gsps
allgwbasins <- st_read("Boundaries/i08_B118_CA_GroundwaterBasins/i08_B118_CA_GroundwaterBasins.shp") %>% st_transform(., crs=merc)
sgmabasins <- st_read("Boundaries/GSP_submitted/GSA_Master.shp") %>% st_transform(., crs=merc)
gspcoda <- read_csv("Boundaries/gsp_coda.csv")
gsps <- left_join(sgmabasins, gspcoda, by="GSP.ID")
gsps <- filter(gsps, is.na(GSP.NAME)==FALSE)
gsps$gsp_area <- st_area(gsps)
# dacs
dacs <- st_read("Boundaries/census_data_disadvantaged_communities_2018/DAC_Pl18.shp") %>% filter(DAC18 == "Y") %>% st_transform(., merc) %>% mutate(area = st_area(.))

dcs <- st_intersection(dacs, gsps)
dcs$dcs_area <- st_area(dcs)
dcs$perc_overlap <- dcs$dcs_area / dcs$area
dcs <- filter(dcs, as.numeric(perc_overlap) > .5)

# MTs
mn <- read.csv("MTs/CentralValleyMTs.csv") %>% filter(is.na(Long)==FALSE & is.na(Lat)==FALSE & is.na(MT_dtw)==FALSE)
mts <- st_as_sf(mn, coords = c("Long","Lat"), crs=4326) %>% st_transform(., crs=merc)
mts <- st_intersection(mts, gsps)
mtjoindata <- mts %>% st_drop_geometry(.) %>% dplyr::select(GSP_ID, GSP_Nam) %>% unique(.)

mt_sp <- as_Spatial(mts)

gsps <- left_join(gsps, mtjoindata, by=c("GSP_ID"))
gsps_sp <- as_Spatial(gsps)

# domestic wells
richsdws <- read_rds("DomesticWells/domcv6_mean_gw_with_beta_GF_CI.rds") %>% st_as_sf(.) %>% st_transform(., st_crs(gsps)) %>% filter(year >=1990)
dw <- st_intersection(richsdws, gsps)
hi <- dw %>% group_by(GSP.NAME) %>% summarise(`Number of Domestic Wells` = length(unique(WCRNumber)),`Average TCD` = mean(TotalCompletedDepth, na.rm=TRUE), `Average Pump Depth` = mean(pump_loc), `Fraction of All DWs`=100*length(unique(WCRNumber))/nrow(dw)) %>% st_drop_geometry(.)

# public supply wells
psws <- read.csv("DomesticWells/gama_location_construction_v2.csv") 
psws <-  st_as_sf(psws, coords = c("GM_LONGITUDE","GM_LATITUDE"), crs=4326) %>% st_transform(., crs=merc)
psw <- psws %>% filter(is.na(GM_WELL_DEPTH_FT)==FALSE & GM_WELL_DEPTH_FT >0 & GM_WELL_DEPTH_FT < 3000) %>% st_intersection(., gsps)

# ALL DRINKING WATER WELLS
dws <- dw %>% mutate(type = "domestic") %>% select(WCR = WCRNumber, type, BASIN, GSP.NAME, TCD = TotalCompletedDepth, pump_loc, top, bot)
psww <- psw %>% mutate(pump_loc = NA, type = "public supply") %>% 
  select(WCR = GM_WELL_ID, type, BASIN, GSP.NAME, TCD = GM_WELL_DEPTH_FT, pump_loc, top = GM_TOP_DEPTH_OF_SCREEN_FT, bot = GM_BOTTOM_DEPTH_OF_SCREEN_FT)

dwws <- rbind(dws, psww) %>% filter(TCD > 0 & TCD < 5000)


# groundwater levels
cgwl_raster <- read_rds("InterpolationGWLevels/cgwl_raster.rds")
mt_raster <- read_rds("InterpolationGWLevels/mt_raster.rds")

#### RESULT ONE - WELLS GO DRY ####

cgwl_at_dw <- raster::extract(cgwl_raster, dwws) # intersect to get value of current water level at well points
cad <- cbind(dwws, cgwl_at_dw) 
cad <- cad[!is.na(cad$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
cad$tcd_dry <- ifelse(cad$TCD <= cad$cgwl_at_dw, "Failing", "Active")
activedw <- cad[cad$tcd_dry == "Active", ]
activedw <- activedw[!is.na(activedw$WCR), ]

length(unique(activedw$WCR)) / length(unique(cad$WCR)) # percent of wells whose TCD is below current groundwater levels (useable wells)

#### dry well analysis ####
# TCD
mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
mad <- cbind(activedw, mt_at_dw) 
mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)

mad$dry <- ifelse(mad$TCD <= mad$mt_at_dw, "Failing", "Active")
table(mad$dry)
table(mad$dry)/length(unique(mad$WCR))

# pump location
pls <- filter(mad, is.na(mad$pump_loc)==FALSE & mad$dry == "Active")
pls$pldry <- ifelse(pls$pump_loc <= pls$mt_at_dw, "Failing", "Active")
table(pls$pldry)
table(pls$pldry)/length(unique(pls$WCR))

# bottom of well screen
bots <- filter(mad, is.na(mad$bot)==FALSE & mad$bot >0 & mad$dry == "Active")
bots$botdry <- ifelse(bots$bot <= bots$mt_at_dw, "Failing", "Active")
table(bots$botdry)
table(bots$botdry)/length(unique(bots$WCR))

# fix wells whose BOT is > TCD or < TOP
mad$bot20 <- mad$bot - 20
mad$TCD20 <- mad$TCD - 20 
mad$TCD20 <- ifelse(mad$TCD20 <= mad$bot20 & is.na(mad$bot20)==FALSE, mad$bot20, mad$TCD20)

# all 
mad$dry <- ifelse(mad$mt_at_dw >= mad$bot20 & is.na(mad$bot20)==FALSE, "BOTdry",
                  ifelse(mad$mt_at_dw >= mad$TCD20 & is.na(mad$TCD20)==FALSE, "TCDdry",
                                ifelse(mad$TCD > mad$mt_at_dw , "Active", 
                                       ifelse(is.na(mad$bot)==TRUE | is.na(mad$TCD) ==TRUE, "notavailable", "uhoh"))))
                  
table(mad$dry)
table(mad$DRY)/length(unique(mad$WCR)) 

mad$DRY <- ifelse((is.na(mad$bot20)==FALSE & mad$mt_at_dw >= mad$bot20) | mad$mt_at_dw >= mad$TCD20, "dry", "active")
table(mad$DRY)


dry <- mad[mad$DRY == "dry", ]

#### FIGURES ####
##### 1. DISTRIBUTION OF DRY WELLS #####
# rasterize dry wells
basins <- gsps %>% group_by(BASIN) %>%  st_buffer(100) %>% summarise(geometry = st_union(geometry))
ggplot(data=basins)+
  geom_sf()

r <- raster(gsps_sp)
res(r) <- 1610*6 # township-level grids (1,609 meters in a mile)
r <- rasterize(gsps_sp, r)
plot(r, col="grey90", lwd=10)
quads <- as(r, "SpatialPolygons")

test_spdf <- as(r, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

tcddry <- dry[!is.na(dry$WCR),]

all <- mad[!is.na(mad$WCR),]
a <- as_Spatial(all)
a <- rasterize(coordinates(a), r, fun="count")
a <- as(a, "SpatialPolygonsDataFrame")
a <- st_as_sf(a)

a$pt_count <- lengths(st_intersects(a, tcddry))
a$perc <- 100*a$pt_count/a$layer

a$dw_pl <- ifelse(a$perc == 0, "0",
                  ifelse(a$perc > 0 & a$perc <=25, "1 - 25", 
                    ifelse(a$perc > 25 & a$perc <= 50, "26 - 50", 
                           ifelse(a$perc > 50 & a$perc <= 75, "51 - 75", 
                                  ifelse(a$perc > 75, "76 - 100", "NA")))))
a <- filter(a, dw_pl != "0")
ints <- st_intersection(a, st_make_valid(basins))
ints$dw_pl <- factor(ints$dw_pl, levels=c("1 - 25", "26 - 50", "51 - 75", "76 - 100"))

ggplot() +
  geom_tile(data=test_df, aes(x=x, y=y), fill="grey80", col="grey60")+
  geom_sf(data=ints, aes(fill=dw_pl), col=NA)  +
  geom_sf(data=basins, fill=NA, col="grey50", lwd=.6) +
  #geom_sf(data=active, cex=.1, col="grey80") +
  scale_fill_brewer(palette = "OrRd") + 
  theme_void(base_size = 16)+ 
  labs(fill = "% of fully dewatered wells\nper 36 square miles\n(township grids)")

##### 2. NUMBER OF DRY WELLS PER BASIN ####
mad$topdry <- ifelse(mad$DRY == "dry", "dry", 
                     ifelse(mad$top == 0 & is.na(mad$top)==FALSE, "not available",
                     ifelse(mad$top <= mad$mt_at_dw & is.na(mad$top)==FALSE, "topdry",
                 ifelse(is.na(mad$top)==TRUE, "not available", "active"))))


gddw <- mad %>% distinct(WCR, .keep_all = TRUE) %>% group_by(BASIN) %>% summarise(AllWells = length(unique(WCR)),
                                                 Active = length(unique(WCR[DRY=="active"])),
                                                 PartiallyDry = length(unique(WCR[topdry=="topdry"])),
                                                 FullyDry = length(unique(WCR[DRY=="dry"])))
gddw$Active <- gddw$Active - gddw$PartiallyDry
gddw$percfail <- 100*((gddw$PartiallyDry+gddw$FullyDry)/gddw$AllWells)
gddw$percactive <- 100*(gddw$Active / gddw$AllWells)

bplot <- as.data.frame(st_drop_geometry(gddw)) %>% 
  filter(percfail >0) %>%
  dplyr::select(., BASIN, Active=Active, `Partially Dewatered`=PartiallyDry, `Fully Dewatered`=FullyDry) %>%
  gather(WellStatus, value, -c(BASIN))

bplot$WellStatus_f <- factor(bplot$WellStatus , levels=c("Fully Dewatered", "Partially Dewatered","Active"))

gddw <- filter(gddw, percfail>0)
ggplot(gddw, aes( y=percfail, x=reorder(BASIN,percfail))) + 
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

bbplot <- bplot %>% group_by(BASIN)
ggplot(bbplot, aes(y=value, x=reorder(BASIN,value), fill=WellStatus_f)) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  theme_minimal()+
  labs(x ="Groundwater Basin", y = "Number of Drinking Water Wells", fill="Well Status")+
  scale_fill_manual(values = c(red2, orange3, blue3))+
  theme(axis.title.x = element_text(color="grey20", size=16),
        axis.title.y = element_text(color="grey20", size=16),
        axis.text = element_text(color="grey20", size=14),
        legend.text = element_text(color="grey20", size=14),
        legend.title = element_text(color="grey20", size=14))+
  theme(legend.position = c(0.75, 0.5))

##### 3. GROUNDWATER LEVEL DECLINE ####
mts <- st_intersection(gsps, mts)

cgwl_at_mts <- raster::extract(cgwl_raster, mts)
mtz <- cbind(mts, cgwl_at_mts) 

mt_at_mts <- raster::extract(mt_raster, mtz) # intersect to get value of current water level at well points
mtzz <- cbind(mtz, mt_at_mts) 

mtzz$decline <- as.numeric(mtzz$MT_dtw)- mtzz$cgwl_at_mts

mtzz %>% group_by(BASIN) %>% filter(length(unique(Well_ID)) > 1) %>%
ggplot(data=., aes(x=reorder(BASIN, decline, na.rm=TRUE), y=decline))+
  geom_jitter(color = orange3,
              fill = 4,
              alpha = 0.25)+
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 2,
    color = red3,
    fill = red3,
    alpha = .8,
  )+
  coord_flip(ylim=c(-100, 350))+
  theme_light()+
  xlab("Groundwater Basin") +
  ylab("Change in Groundwater Level (MT - CGWL)")
 
##### 4. DACS COVERED #####
# percent of dacs with at least 50% of their area covered
bufferdistance = 5632.69
buffer <- st_buffer(mts, dist=bufferdistance, joinStyle = "BEVEL", endCapStyle = "ROUND")
buffer <- st_union(buffer)

dinmt <- st_intersection(dcs, buffer)
dinmt$int_area <- st_area(dinmt)
dinmt$percarea <- dinmt$int_area/dinmt$area

wellsindacs <- st_intersection(mad, dinmt)

wid <- wellsindacs %>% group_by(NAMELSAD) %>% summarise(percCoverage = unique(as.numeric(percarea)),
                                                        AllWells = length(unique(WCR)),
                                                               Active = length(unique(WCR[DRY=="active"])),
                                                               PartiallyDry = length(unique(WCR[topdry=="topdry"])),
                                                               FullyDry = length(unique(WCR[DRY=="dry"])))
wid$Active <- wid$Active - wid$PartiallyDry
wid$percfail <- 100*((wid$PartiallyDry+wid$FullyDry)/wid$AllWells)
wid$percactive <- 100*(wid$Active / wid$AllWells)

wid$percactive_cat <- ifelse()

# percent of wells in DACs that were not covered 


# percent of all wells not covered, percent of DACs not covered

