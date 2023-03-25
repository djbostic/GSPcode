# load packages
library(sp)        # spatial objects
library(raster)    # for raster objects
library(here)
library(sf)
library(automap)
library(lubridate)
library(gstat)
library(RColorBrewer)
library(Cairo)
library(cowplot)
library(magick)
library(tidyverse) # general purpose data science toolkit
library(patchwork)



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

##### LOAD DATA #####
# ca
ca <- st_read("Boundaries/cb_2018_us_state_500k/cb_2018_us_state_500k.shp") %>% filter(STUSPS == "CA") %>% st_transform(., crs=merc)

# gsps
allgwbasins <- st_read("Boundaries/i08_B118_CA_GroundwaterBasins/i08_B118_CA_GroundwaterBasins.shp") %>% st_transform(., crs=merc)

sgmabasins <- st_read("Boundaries/GSP_submitted/GSA_Master.shp") %>% st_transform(., crs=merc)
gspcoda <- read_csv("Boundaries/gsp_coda.csv") %>% mutate(GSP.ID = as.character(GSP.ID))
gsps <- left_join(sgmabasins, gspcoda, by="GSP.ID")
gsps <- filter(gsps, is.na(GSP.NAME)==FALSE)
gsps$gsp_area <- st_area(gsps)
gsps_sp <- as_Spatial(gsps)

basins <- gsps %>% group_by(BASIN) %>%  st_buffer(100) %>% summarise(geometry = st_union(geometry))

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
mtjoindata <- mts %>% st_drop_geometry(.) %>% dplyr::select(GSP.ID, GSP.NAME) %>% unique(.)

# domestic wells
richsdws <- read_rds("DomesticWells/domcv6_mean_gw_with_beta_GF_CI.rds") %>% st_as_sf(.) %>% st_transform(., st_crs(gsps)) %>% filter(year >=1990)
dw <- st_intersection(richsdws, gsps)
hi <- dw %>% group_by(GSP.NAME) %>% summarise(`Number of Domestic Wells` = length(unique(WCRNumber)),`Average TCD` = mean(TotalCompletedDepth, na.rm=TRUE), `Average Pump Depth` = mean(pump_loc), `Fraction of All DWs`=100*length(unique(WCRNumber))/nrow(dw)) %>% st_drop_geometry(.)

# public supply wells
psws <- read.csv("DomesticWells/gama_location_construction_v2.csv") %>% filter(., GM_WELL_CATEGORY == "MUNICIPAL" | GM_WELL_CATEGORY == "WATER SUPPLY, OTHER")
psws <-  st_as_sf(psws, coords = c("GM_LONGITUDE","GM_LATITUDE"), crs=4326) %>% st_transform(., crs=merc)
psw <- psws %>% filter(is.na(GM_WELL_DEPTH_FT)==FALSE & GM_WELL_DEPTH_FT >0 & GM_WELL_DEPTH_FT < 5000) %>% st_intersection(., gsps)

# ALL DRINKING WATER WELLS
dws <- dw %>% mutate(type = "domestic") %>% dplyr::select(WCR = WCRNumber, type, BASIN, GSP.NAME=GSP.NAME, TCD = TotalCompletedDepth, pump_loc, top, bot)
psww <- psw %>% mutate(pump_loc = NA, type = "public supply") %>% 
  dplyr::select(WCR = GM_WELL_ID, type, BASIN, GSP.NAME=GSP.NAME, TCD = GM_WELL_DEPTH_FT, pump_loc, top = GM_TOP_DEPTH_OF_SCREEN_FT, bot = GM_BOTTOM_DEPTH_OF_SCREEN_FT)

dwws <- rbind(dws, psww) %>% filter(TCD > 0 & TCD < 5000)

# groundwater levels
cgwl_raster <- read_rds("InterpolationGWLevels/cgwl_raster.rds")
mt_raster <- read_rds("InterpolationGWLevels/mt_raster.rds")

#### FIGURE 1 - GW LEVELS ####
cgwl_at_mts <- raster::extract(cgwl_raster, mts) # intersect to get value of current water level at well points
mtz <- cbind(mts, cgwl_at_mts) 

mt_at_mts <- raster::extract(mt_raster, mtz) # compare interpolated MT water level with MTs set at same location
mtzz <- cbind(mtz, mt_at_mts) 

mtzz$mt_cgwl <- as.numeric(mtzz$MT_dtw) - mtzz$cgwl_at_mts
mtzz$decline <- mtzz$cgwl_at_mts-as.numeric(mtzz$MT_dtw)
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
#mtzz <- filter(mtzz, is.na(mtzz$bins)==FALSE)

gwl_map <- ggplot()+
  geom_sf(data=basins, fill="gray95", col=NA)+
  geom_point(data=mtzz, alpha=.85,
             aes(color = bins, geometry = geometry),
             stat = "sf_coordinates")+
  scale_color_brewer("Change in Groundwater Level \n(CGWL-MT)",
    palette = "BrBG",
    direction = -1
  )+ 
  geom_sf(data=basins, fill=NA, col="gray40", lwd=.3)+
  theme_void()+
theme(legend.title=element_text(size=14), 
      legend.text=element_text(size=14)) +
  guides(col = guide_legend(override.aes = list(size=3.5)))
gwl_map

library(png)
inset <- readPNG("Results/Images/InsetMap2.png", native=TRUE)

library(patchwork)
gwl_map2 <- gwl_map + inset_element(p = inset,
                        left = 0,
                        bottom = 0,
                        right = .5,
                        top = .3)

#ggsave(plot=gwl_map2, filename="Results/Images/FIGURE1A.png", bg = "transparent", width=6.70, height=5.58, type = "cairo", device = "png")

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

elevchange <- 
  ggplot()+
  geom_jitter(data=elchmt, aes(x=reorder(BASIN.x,-MedMT, na.rm=TRUE), y=decline, col=bins.x), fill = 4, alpha = 0.5, size = 1.5)+
  geom_point(data=mts_basin, aes(x=reorder(BASIN, -MedMT, na.rm=TRUE), y=MedMT, fill=bins), pch=21, alpha=1, size=3, col="gray90", show.legend = FALSE)+
  scale_color_brewer("Change in Groundwater Level \n(CGWL-MT)",
                     palette = "BrBG",
                     direction = -1)+ 
  scale_fill_brewer("Change in Groundwater Level \n(CGWL-MT)",
                     palette = "BrBG",
                     direction = -1
  )+ 
  coord_flip(ylim=c(100,-350))+
  theme_light()+
  xlab("Groundwater Subbasin") +
  ylab("Change in Groundwater Level (CGWL - MT)")+
  theme(axis.text.x =element_text(size=12), 
        axis.text.y =element_text(size=10),
        axis.title = element_text(size=12))+ 
  guides(col = guide_legend(override.aes = list(size=2)))
elevchange

#ggsave(plot=elevchange, filename="Results/Images/FIGURE1B.png", bg = "transparent", width=9, height=5.58, type = "cairo", device = "png")

# plot both next to each other
#p <- plot_grid(gwl_map, elevchange, labels = "AUTO")
#save_plot("Results/Images/Figure1a_1b_ElevChange.png", p, ncol = 2, type="cairo", bg="transparent", base_width=12, base_height=10)

# How many GSPs set groundwater levels over 200 feet below CGWL?
ftbelow200 <- mtzz %>% filter(decline <= -200)
length(unique(ftbelow200$Well_ID))
length(unique(ftbelow200$GSP.NAME))
length(unique(ftbelow200$BASIN))

#### FIGURE 2 - DRY WELLS ####
##### all well analysis ######
wellanalysis <- function(x=dwws){
  cgwl_at_dw <- raster::extract(cgwl_raster, x) # intersect to get value of current water level at well points
  cad_dw <- cbind(x, cgwl_at_dw) 
  cad_dw <- cad_dw[!is.na(cad_dw$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
  cad_dw$tcd_dry <- ifelse(cad_dw$TCD <= cad_dw$cgwl_at_dw, "Failing", "Active")
  activedw <- cad_dw[cad_dw$tcd_dry == "Active", ]
  activedw <- activedw[!is.na(activedw$WCR), ]
  print(length(unique(activedw$WCR)) / length(unique(cad_dw$WCR))) # percent of wells whose TCD is below current groundwater levels (useable wells)
  
  # TCD
  mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
  mad <- cbind(activedw, mt_at_dw) 
  mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
  mad$tcddry <- ifelse(mad$TCD <= mad$mt_at_dw, "Failing", "Active")
  print(table(mad$tcddry))
  print(table(mad$tcddry)/length(unique(mad$WCR)))
  
  # bottom of well screen
  bots <- filter(mad, is.na(mad$bot)==FALSE & mad$bot > 0 & mad$tcddry == "Active")
  bots$botdry <- ifelse(bots$bot <= bots$mt_at_dw, "Failing", "Active")
  print(table(bots$botdry))
  print(table(bots$botdry)/length(unique(bots$WCR)))
  
  # all 
  mad$top <- ifelse(is.na(mad$top)==TRUE, 0, mad$top)
  mad$dry <- ifelse((mad$mt_at_dw >= mad$TCD & mad$TCD > 0), "TCDdry",
                       ifelse((mad$top > 0 & mad$mt_at_dw >= mad$top & mad$mt_at_dw < mad$TCD), "topdry",
                              ifelse(mad$pump_loc > 0 & mad$mt_at_dw >= mad$pump_loc & mad$mt_at_dw < mad$TCD & mad$mt_at_dw < mad$top, "pumpdry", "active")))
  print(table(mad$dry)/length(unique(mad$WCR)))
  return(mad)
}

##### domestic well analysis ######
mad_dw <- wellanalysis(dws)

# all for plot
# rasterize dry wells
r <- raster(gsps_sp)
res(r) <- 1610*6 # township-level grids (1,609 meters in a mile)
r <- rasterize(gsps_sp, r)
plot(r, col="grey90", lwd=10)
quads <- as(r, "SpatialPolygons")
quadss <- st_as_sf(quads)

# DRY DOMESTIC WELLS BAR PLOT
reorgdry <- function(x=mad_dw, y="DW or PSW"){
gddw <- x %>% distinct(WCR, .keep_all = TRUE) %>% group_by(type, BASIN) %>% summarise(AllWells = length(unique(WCR)),
                                                                                                     Active = length(unique(WCR[dry=="active"])),
                                                                                                     PartiallyDry = length(unique(WCR[dry=="topdry"])),
                                                                                                     FullyDry = length(unique(WCR[dry=="TCDdry"])),
                                                                                           PumpDry=length(unique(WCR[dry=="pumpdry"])))

gddw$percfail <- 100*((gddw$PartiallyDry+gddw$FullyDry+gddw$PumpDry)/gddw$AllWells)
gddw$check <- 100*((gddw$Active+gddw$PartiallyDry+gddw$FullyDry+gddw$PumpDry)/gddw$AllWells)

bplot <- as.data.frame(st_drop_geometry(gddw)) %>% 
  filter(percfail >0) %>%
  dplyr::select(., BASIN, Active=Active, `Partially Dewatered`=PartiallyDry, `Fully Dewatered`=FullyDry, `Pump Dewatered`=PumpDry) %>%
  gather(WellStatus, value, -c(BASIN))

bplot$WellStatus_f <- factor(bplot$WellStatus , levels=c("Pump Dewatered", "Fully Dewatered", "Partially Dewatered","Active"))

dryplot <- ggplot(bplot, aes(y=value, x=reorder(BASIN,value, sum), fill=WellStatus_f)) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  theme_minimal()+
  labs(x ="Groundwater Basin", y = y, fill="Well Status")+
  scale_fill_manual(values = c(orange1, red2, orange3, blue3))+
  theme(axis.title.x = element_text(color="grey20", size=16),
        axis.title.y = element_text(color="grey20", size=16),
        axis.text = element_text(color="grey20", size=14),
        legend.text = element_text(color="grey20", size=14),
        legend.title = element_text(color="grey20", size=14))+
  theme(legend.position = c(0.75, 0.5))
return(dryplot)
}


# DOMESTIC WELLS ACTIVE, PD, FD, PD
domwells<- reorgdry(x=mad_dw, y="Domestic Wells")

# filter for type here
blues <- c("#F2EDE4", "lightblue", "cornflowerblue", "blue3")
oranges <- c("#F2EDE4", "#F9AB74", "#EA7A04","#993F00")
reds <- c("#F2EDE4", "#FA6D6D", "#C70039", "#900C3F")
tans <- c("#F2EDE4",  "#D9CAAD", "#BFA575", "#735D36")


wellplot <- function(x=mad_dw, y="active", z=blues, w="Percent of active \ndomestic wells\nper 36 square miles\n(township grids)", u="top") {
activo <- x[x$dry == y, ]
all <- x[!is.na(x$WCR),]
a <- as_Spatial(all)
a <- rasterize(coordinates(a), r, fun="count")
a <- as(a, "SpatialPolygonsDataFrame")
a <- st_as_sf(a)


a$pt_count_act <- lengths(st_intersects(a, activo))
a$perc_act <- 100*a$pt_count_act/a$layer

a$dw_pl_act <- ifelse(a$perc_act == 0, "0",
                  ifelse(a$perc_act > 0 & a$perc_act <=25, "1 - 25%", 
                         ifelse(a$perc_act > 25 & a$perc_act <= 50, "26 - 50%", 
                                ifelse(a$perc_act > 50 & a$perc_act <= 75, "51 - 75%", 
                                       ifelse(a$perc_act > 75, "76 - 100%", "NA")))))
v <- filter(a, dw_pl_act != "0")
ints <- st_intersection(v, st_make_valid(basins))
ints$dw_pl_act <- factor(ints$dw_pl_act, levels=c("1 - 25%", "26 - 50%", "51 - 75%", "76 - 100%"))

dw_map <- ggplot() +
  geom_sf(data=a, fill="grey80", col="grey90")+
  geom_sf(data=ints, aes(fill=dw_pl_act), col=NA)  +
  geom_sf(data=basins, fill=NA, col="grey40", lwd=.4) +
  scale_fill_manual(values = z)+
  ggtitle(w)+
  theme_void(base_size=8)+
  theme(plot.title = element_text(size=12, hjust=.5),
        legend.position = u,
        legend.key.width = unit(.4,"cm"),
        legend.key.height = unit(.4,"cm"),
        legend.spacing = unit(0.25,"cm"),
        legend.title = element_blank(),
        legend.justification = "center") +
  guides(fill=guide_legend(title.position="top", label.position="bottom")) +
  coord_sf()
}

active <- wellplot(z=blues, w="% Active")
topd <- wellplot(x=mad_dw, y="topdry", z=oranges, w="% Partially Dewatered")
tcdd <- wellplot(y="TCDdry", z=reds, w="% Fully Dewatered")
pumpd <- wellplot(y="pumpdry", z=tans, w="% Pump Dewatered")

alldw <- active + topd + tcdd + pumpd +  plot_layout(ncol = 4) + plot_annotation(title="Domestic Wells", theme=theme(plot.title = element_text(size=14, hjust=.5)), tag_levels = 'A')
ggsave(plot=alldw, filename="Results/Images/FIGURE3_DW.png", bg = "transparent", type = "cairo", device = "png")

##### public supply well analysis ######
mad_psw <- wellanalysis(x=psww)
# DRY PUBLIC SUPPLY WELLS BAR PLOT
psws <- reorgdry(x=mad_psw, y="Public Supply Wells")+
  scale_fill_manual(values=c(red2, orange3, blue3),
    limits = c("Fully Dewatered", "Partially Dewatered", "Active"))

# PSW ACTIVE, PD, FD
activeP <- wellplot(x=mad_psw, z=blues, w="% Active")
topdP <- wellplot(x=mad_psw, y="topdry", z=oranges, w="% Partially Dewatered")
tcddP <- wellplot(x=mad_psw, y="TCDdry", z=reds, w="% Fully Dewatered")

allpsw <- activeP + topdP + tcddP + plot_annotation(title="Public Supply Wells", theme=theme(plot.title = element_text(size=14, hjust=.5)), tag_levels = 'A')
ggsave(plot=allpsw, filename="Results/Images/FIGURE3_PSW.png", bg = "transparent", type = "cairo", device = "png")

