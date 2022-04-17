# load packages
library(tidyverse) # general purpose data science toolkit
library(sf)        # spatial objects
library(raster)    # for raster objects
library(readr)
library(here)
library(lubridate)


# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

#### load data ####


# download wells from OSWCR
wells <- read_csv("DomesticWells/oswcr_bulkdatadownload_04112022/OSWCR.csv") 
wellz <- wells %>% mutate(DATEWORKENDED = as.Date(DATEWORKENDED, "%m/%d/%Y"), year=year(DATEWORKENDED)) %>% filter(year>1960 & is.na(DECIMALLONGITUDE)==FALSE & is.na(DECIMALLATITUDE)==FALSE & !grepl("/", DECIMALLATITUDE) & !grepl("/", DECIMALLONGITUDE))

wellz <- st_as_sf(wellz, coords = c("DECIMALLONGITUDE","DECIMALLATITUDE"), crs=4326) %>% st_transform(., crs=merc)
domesticwells_in_gsps <- filter(wellz, grepl("domestic", PLANNEDUSEFORMERUSE, ignore.case = TRUE))

# how many wells have TCDs, TOP of Screens, and Bottom of Screens?
# TCD
dw_tcd <- filter(domesticwells_in_gsps, is.na(TOTALCOMPLETEDDEPTH)==FALSE & TOTALCOMPLETEDDEPTH > 0)
length(unique(dw_tcd$WCRNUMBER))

dw_tos <- filter(domesticwells_in_gsps, is.na(TOPOFPERFORATEDINTERVAL)==FALSE & TOPOFPERFORATEDINTERVAL > 0)
length(unique(dw_tos$WCRNUMBER))

dw_bos <- filter(domesticwells_in_gsps, is.na(BOTTOMOFPERFORATEDINTERVAL)==FALSE & BOTTOMOFPERFORATEDINTERVAL > 0)
length(unique(dw_bos$WCRNUMBER))

dw_all <- filter(domesticwells_in_gsps, is.na(BOTTOMOFPERFORATEDINTERVAL)==FALSE & BOTTOMOFPERFORATEDINTERVAL > 0 &
                   is.na(TOTALCOMPLETEDDEPTH)==FALSE & TOTALCOMPLETEDDEPTH > 0 &
                   is.na(TOPOFPERFORATEDINTERVAL)==FALSE & TOPOFPERFORATEDINTERVAL > 0)
length(unique(dw_all$WCRNUMBER))

dw_all <- filter(dw_all, year > 1975 & year < 2017)

domwells_gsps <- st_intersection(dw_all, gsps)
st_write(domesticwells_in_gsps, "DomesticWells/dw_gsp.shp", append=FALSE)

# see how many wells there are through the years

domsub <- domesticwells_in_gsps %>% filter(year > 1965 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 1970 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 1975 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 1980 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 1985 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub90 <- domesticwells_in_gsps %>% filter(year > 1990 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 1995 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

domsub <- domesticwells_in_gsps %>% filter(year > 2000 & year < 2017)
length(unique(domsub$WCRNUMBER))
mean(domsub$TOTALCOMPLETEDDEPTH, na.rm=TRUE)

hi <- domesticwells_in_gsps %>% group_by(GSP.ID) %>% summarise(mean(TOTALCOMPLETEDDEPTH, na.rm=TRUE), length(unique(WCRNUMBER)))
