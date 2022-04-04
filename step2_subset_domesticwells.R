# load packages
library(tidyverse) # general purpose data science toolkit
library(sf)        # spatial objects
library(raster)    # for raster objects
library(readr)

# set working directory
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/Volumes/GoogleDrive/My Drive/Graduate School/GSP_Analy/Organized_CnD/")

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

#### load data ####
# all gsp boundaries
gsps <- st_read("Boundaries/all_cobs/All_CO_GSPs.shp") %>% st_transform(., merc) %>% filter(., GSP_Name != "Las Posas")

# add dictionary of coastal, CV, and southern basins
gsps$region <- ifelse(gsps$GSP_Name == "Alpaugh" | 
                      gsps$GSP_Name == "Central Kings" | 
                      gsps$GSP_Name == "Chowchilla" | 
                      gsps$GSP_Name == "Fresno County - Delta Mendota" | 
                      gsps$GSP_Name == "Delano-Earlimart" | 
                      gsps$GSP_Name == "East Kaweah" | 
                      gsps$GSP_Name == "Eastern San Joaquin" | 
                      gsps$GSP_Name == "Eastern Tule" | 
                      gsps$GSP_Name == "Farmers WD" | 
                      gsps$GSP_Name == "Greater Kaweah" | 
                      gsps$GSP_Name == "Northern Central DM" | 
                      gsps$GSP_Name == "Kern Groundwater Authority" | 
                      gsps$GSP_Name == "Grasslands" |
                      gsps$GSP_Name == "Henry Miller" |
                      gsps$GSP_Name == "James" |  
                      gsps$GSP_Name == "Kern River" |
                      gsps$GSP_Name == "Kings River East" |
                      gsps$GSP_Name == "Lower Tule River ID" |
                      gsps$GSP_Name == "Mid-Kaweah" |
                      gsps$GSP_Name == "Madera" |
                      gsps$GSP_Name == "McMullin Area" |
                      gsps$GSP_Name == "Merced" |
                      gsps$GSP_Name == "North Fork Kings" |
                      gsps$GSP_Name == "New Stone" |
                      gsps$GSP_Name == "North Kings" |
                      gsps$GSP_Name == "Olcese" |
                      gsps$GSP_Name == "Pixley ID" |
                      gsps$GSP_Name == "Root Creek" |
                      gsps$GSP_Name == "SJREC" |  
                      gsps$GSP_Name == "South Kings" |
                      gsps$GSP_Name == "TriCounty WA - Tule" |
                      gsps$GSP_Name == "Tulare" |
                      gsps$GSP_Name == "Westlands" |
                      gsps$GSP_Name == "Buena Vista" |
                      gsps$GSP_Name == "Aliso WD", "San Joaquin Valley", 
                ifelse(gsps$GSP_Name == "Oxnard" | 
                       gsps$GSP_Name == "Pleasant Valley" |
                       gsps$GSP_Name == "180_400ftAq - Salinas" |
                       gsps$GSP_Name == "Santa Cruz Mid-County", "Coastal", 
                 ifelse(gsps$GSP_Name == "Cuyama" | 
                        gsps$GSP_Name == "Paso Robles" |
                        gsps$GSP_Name == "Indian Wells" , "Central", "NA")))

ca <- st_read("Boundaries/ca/CA_State_TIGER2016.shp")
ca <- st_transform(ca, merc)
ggplot()+
  geom_sf(data=ca, fill = "transparent") +
  geom_sf(data=gsps, aes(fill = gsps$region)) +
  scale_fill_manual("Region", values = c("San Joaquin Valley" = "#70AB53", "Coastal" = "#535EAB", "Central" = "#AB6353")) +
  theme_void(base_size = 16)
gsps <- st_zm(gsps, drop=TRUE)
st_write(gsps, "Output/gsps_with_region.shp", delete_dsn =TRUE)

# Central Valley domestic wells with estimated pump depths (from Rich)
domcv <- read_rds("2Subset_DWs/Data/domcv6_mean_gw_with_beta_GF_CI.rds")
domcvda <- as(domcv, "sf")

#domcvda <- filter(domcvda, year > 1960 & year < 2017)

# california cleaned domestic wells (from Rich)
noncvgsp <- filter(gsps, region != "San Joaquin Valley")

dom <- st_read("2Subset_DWs/Data/wells/domestic_wells_ca.shp", stringsAsFactors=FALSE) %>% st_transform(., merc)
dom$year <- as.numeric(dom$year)
domsub <- dom %>% st_intersection(.,noncvgsp) 

domcvda <- filter(domcvda, year > 1960 & year < 2017)
domsub <- domsub %>% filter(year > 1960 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1965 & year < 2017)
domsub <- domsub %>% filter(year > 1965 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1970 & year < 2017)
domsub <- domsub %>% filter(year > 1970 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1975 & year < 2017)
domsub <- domsub %>% filter(year > 1975 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1980 & year < 2017)
domsub <- domsub %>% filter(year > 1980 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1985 & year < 2017)
domsub <- domsub %>% filter(year > 1985 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1990 & year < 2017)
domsub90 <- domsub %>% filter(year > 1990 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 1995 & year < 2017)
domsub <- domsub %>% filter(year > 1995 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

domcvda <- filter(domcvda, year > 2000 & year < 2017)
domsub <- domsub %>% filter(year > 2000 & year < 2017)
length(unique(domcvda$WCRNumber)) + length(unique(domsub$WCRNmbr))
mean(mean(domcvda$TotalCompletedDepth, na.rm=TRUE), mean(domsub$TtlCmpD, na.rm=TRUE))

#write_rds(domsub, "Output/domesticwells_allcobs_75.rds")
#write_rds(domsub90, "Output/domesticwells_allcobs_90.rds")

hi <- domsub %>% group_by(GSP_Name) %>% summarise(mean(TtlCmpD, na.rm=TRUE), length(unique(WCRNmbr)))


# OSWCR csvs Darcy collected(no estimated pump depths)
oswcr_files <- list.files("2Subset_DWs/Data/DWs/", pattern = "\\.csv$")
oswcr <- list()
for(i in 1:length(oswcr_files)){
  oswcr[[i]] = read.csv(file.path("2Subset_DWs/Data/DWs/", oswcr_files[i]), stringsAsFactors = FALSE)
}
oswcr_shp <- lapply(oswcr, function(x) st_as_sf(x, coords = c("lon", "lat"), crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
oswcr_merc <- lapply(oswcr_shp, function(x) st_transform(x, merc))
oswcr_compiled <- lapply(oswcr_merc, function(x) x %>% dplyr::select(WCRNumber, year, type, TotalCompletedDepth))
oswcr_compiled <- do.call(rbind, oswcr_compiled)

# filter oswcr data
oswcr_filtered <- oswcr_compiled %>% filter(year > 1975 & type == "domestic")
plot(gsps)
plot(oswcr_filtered$geometry, add=TRUE, col="red")

# export data to be used elsewhere
write_rds(oswcr_filtered, "Output/oswcr_allcobs.rds")
