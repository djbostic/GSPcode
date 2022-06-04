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

# dacs
dacs <- st_read("Boundaries/census_data_disadvantaged_communities_2018/DAC_Pl18.shp") %>% filter(DAC18 == "Y") %>% st_transform(., merc)
dcs <- st_intersection(dacs, gsps)

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
dws <- dw %>% mutate(type = "domestic") %>% select(WCR = WCRNumber, type, GSP.NAME, TCD = TotalCompletedDepth, pump_loc, top, bot)
psww <- psw %>% mutate(pump_loc = NA, type = "public supply") %>% 
  select(WCR = GM_WELL_ID, type, GSP.NAME, TCD = GM_WELL_DEPTH_FT, pump_loc, top = GM_TOP_DEPTH_OF_SCREEN_FT, bot = GM_BOTTOM_DEPTH_OF_SCREEN_FT)

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
