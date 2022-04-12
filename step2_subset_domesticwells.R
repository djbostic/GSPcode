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
# all gsp boundaries
gsps <- st_read("Boundaries/GSP_submitted/GSA_Master.shp") %>% st_transform(., merc)

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

# download wells from OSWCR
wells <- read_csv("DomesticWells/oswcr_bulkdatadownload_04112022/OSWCR.csv") 
wellz <- wells %>% mutate(DATEWORKENDED = as.Date(DATEWORKENDED, "%m/%d/%Y"), year=year(DATEWORKENDED)) %>% filter(year>1960 & is.na(DECIMALLONGITUDE)==FALSE & is.na(DECIMALLATITUDE)==FALSE & !grepl("/", DECIMALLATITUDE) & !grepl("/", DECIMALLONGITUDE))

wellz <- st_as_sf(wellz, coords = c("DECIMALLONGITUDE","DECIMALLATITUDE"), crs=4326) %>% st_transform(., crs=merc)
wells_in_gsps <- st_intersection(wellz, gsps)
domesticwells_in_gsps <- filter(wells_in_gsps, grepl("domestic", PLANNEDUSEFORMERUSE, ignore.case = TRUE))

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
