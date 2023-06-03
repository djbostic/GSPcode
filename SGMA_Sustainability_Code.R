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
library(png)


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

#dcs <- st_intersection(dacs, gsps)
dcs <- st_intersection(dacs, basins)
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
mtzz <- filter(mtzz, is.na(mtzz$bins)==FALSE)

my_colors <- rev(brewer.pal(9, "BrBG")[c(1:4,7:9)])

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

#ggsave(plot=gwl_map2, filename="Results/Images/FIGURE1A.png", bg = "transparent", width=6.70, height=5.58, type = "cairo", device = "png")

# Estimated groundwater elevation changes (x-axis) at subbasin (y-axis) under proposed MTs
# add median of groundwater levels
mts_basin <- mtzz %>% group_by(GSP.NAME) %>% summarise(MedMT = median(decline)) 
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

mycolors2 <- c(rev(brewer.pal(9, "BrBG")[c(1:3,7:8)]), "black", "black")
elevchange <- 
  ggplot()+
  geom_jitter(data=elchmt, aes(x=reorder(BASIN.x,-MedMT, na.rm=TRUE), y=decline, col=bins.x), fill = 4, alpha = 0.5, size = 1.5, show.legend = FALSE)+
  geom_point(data=mts_basin, aes(x=reorder(BASIN, -MedMT, na.rm=TRUE), y=MedMT, fill=bins), pch=21, alpha=1, size=3, col="black", show.legend = FALSE)+
  scale_color_manual("Change in Groundwater \nLevel (CGWL-MT)",
                     values=my_colors)+ 
  scale_fill_manual("Change in Groundwater Level \n(CGWL-MT)",
                    values=mycolors2
  )+ 
  coord_flip(ylim=c(100,-350))+
  theme_light()+
  xlab("Groundwater Subbasin") +
  ylab("Change in Groundwater Level (CGWL - MT) (ft)")+
  theme(axis.text.x =element_text(size=12), 
        axis.text.y =element_text(size=10),
        axis.title = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))+ 
  guides(col = guide_legend(override.aes = list(size=2, fill=NA)), fill = "none")
elevchange

t1 <- plot_grid(gwl_map2, elevchange, align="hv", labels = c("a.", "b."), label_size=20)
ggsave(plot=t1, filename="Results/Images/FIGURE1B.png", bg = "white", type = "cairo", device = "png", height = 7, width = 14.5)

# How many GSPs set groundwater levels over 200 feet below CGWL?
ftbelow200 <- mtzz %>% filter(decline <= -100)
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
  mad$bot <- ifelse(is.na(mad$bot)==TRUE, 0, mad$bot)
  mad$pump_loc <- ifelse(is.na(mad$pump_loc)==TRUE, 0, mad$pump_loc)
  mad$dry <- ifelse((mad$mt_at_dw >= mad$TCD & mad$TCD > 0), "TCDdry",
                       ifelse(mad$top > 0 & mad$mt_at_dw >= mad$top & mad$mt_at_dw < mad$TCD, "topdry",
                              ifelse(mad$pump_loc > 0 & mad$mt_at_dw >= mad$pump_loc & mad$mt_at_dw < mad$TCD & mad$mt_at_dw < mad$top, "pumpdry", "active")))
  print(table(mad$dry)/length(unique(mad$WCR)))
  return(mad)
}
reorgdry <- function(x=mad_dw, y="DW or PSW", n="Title"){
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
  
  bplot$WellStatus_f <- factor(bplot$WellStatus , levels=c("Pump Dewatered", "Fully Dewatered", "Partially Dewatered", "Active"))
  bplot$BASIN <- substring(bplot$BASIN, 9)
  
  dryplot <- ggplot(bplot, aes(y=value, x=reorder(BASIN,value, sum), fill=WellStatus_f)) + 
    geom_bar(position="stack", stat="identity")+
    coord_flip()+
    theme_minimal(base_size=18)+
    labs(x ="Groundwater Basin", y = y, fill="Well Status")+
    ggtitle(n)+
    scale_fill_manual(values = c(orange1, red2, orange3, blue3), guide=guide_legend(reverse=T))+
    theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
          axis.title.x = element_text(color="grey20", size=18),
          axis.title.y = element_text(color="grey20", size=18),
          axis.text = element_text(color="grey20", size=18),
          legend.text = element_text(color="grey20", size=18),
          legend.title = element_text(color="grey20", size=20))+
    theme(legend.position = c(0.75, 0.5))
  return(dryplot)
}
wellplot <- function(x=mad_dw, y="active", z=blues, w="Percent of active \ndomestic wells\nper 36 square miles\n(township grids)") {
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
    geom_sf(data=ints, aes(fill=dw_pl_act), col="white", lwd=.3)  +
    geom_sf(data=basins, fill=NA, col="grey35", lwd=.4) +
    scale_fill_manual(values = z)+
    ggtitle(w)+
    theme_void()+
    theme(plot.title = element_text(size=18, hjust=.5, face="italic"),
          legend.position = "top",
          legend.key.width = unit(.4,"cm"),
          legend.key.height = unit(.4,"cm"),
          legend.spacing = unit(0.25,"cm"),
          legend.title = element_blank(),
          legend.text = element_text(color="grey20", size=12),
          legend.justification = "center") +
    guides(fill=guide_legend(title.position="top", label.position="bottom")) +
    coord_sf()
}

# all for plot
# rasterize dry wells
r <- raster(gsps_sp)
res(r) <- 1610*6 # township-level grids (1,609 meters in a mile)
r <- rasterize(gsps_sp, r)
plot(r, col="grey90", lwd=10)
quads <- as(r, "SpatialPolygons")
quadss <- st_as_sf(quads) %>% st_crop(., gsps)

##### domestic well analysis ######
mad_dw <- wellanalysis(dws)

# DRY DOMESTIC WELLS BAR PLOT
domwells<- reorgdry(x=mad_dw, y="Number of Domestic Wells", n="a. Domestic Wells")
domwells
# DOMESTIC WELLS ACTIVE, PD, FD, PD
blues <- c("#F2EDE4", "lightblue", "#597980", "#005f73")
oranges <- c("#F2EDE4", "#F9AB74", "#EA7A04","#993F00")
reds <- c("#F2EDE4", "#FA6D6D", "#C70039", "#900C3F")
tans <- c("#F2EDE4",  "#D9CAAD", "#BFA575", "#735D36")

active <- wellplot(z=blues, w="% Active")
topd <- wellplot(x=mad_dw, y="topdry", z=oranges, w="% Partially Dewatered")
tcdd <- wellplot(y="TCDdry", z=reds, w="% Fully Dewatered")
pumpd <- wellplot(y="pumpdry", z=tans, w="% Pump Dewatered")

b1 <- active + topd + tcdd + pumpd +  plot_layout(ncol = 4)
b2 <- plot_grid(domwells, b1, labels = c("1.", "2."), nrow=2)

ggsave(plot=b2, filename="Results/Images/FIGURE3_DWa.png", bg = "transparent", type = "cairo", device = "png", width=16, height=18)

##### public supply well analysis ######
mad_psw <- wellanalysis(x=psww)
# DRY PUBLIC SUPPLY WELLS BAR PLOT
psws <- reorgdry(x=mad_psw, y="Number of Public Supply Wells", n="b. Public Supply Wells")+
  scale_fill_manual(values=rev(c(red2, orange3, blue3)),
    limits = rev(c("Fully Dewatered", "Partially Dewatered", "Active")))

# PSW ACTIVE, PD, FD
activeP <- wellplot(x=mad_psw, z=blues, w="% Active")
topdP <- wellplot(x=mad_psw, y="topdry", z=oranges, w="% Partially Dewatered")
tcddP <- wellplot(x=mad_psw, y="TCDdry", z=reds, w="% Fully Dewatered")

a2 <- (activeP + topdP + tcddP)
a3 <- plot_grid(psws, a2, labels = c("1.", "2."), nrow=2)

ggsave(plot=a3, filename="Results/Images/FIGURE3_PSWb.png", bg = "transparent", type = "cairo", device = "png", width=16, height=18)

#### FIGURE 3 - COVERAGE OF WELLS ####
# MT Buffer Network
bufferdistance = 5632.69 #km
buffer <- st_buffer(mts, dist=bufferdistance, joinStyle = "BEVEL", endCapStyle = "ROUND")
buffer <- st_union(buffer) %>% st_sf(.) %>% st_crop(., basins)

#all psw and dws
wellsinMN <- st_intersection(dwws, buffer)

#number of wells within buffer
100*length(unique(wellsinMN$WCR))/length(unique(dwws$WCR))

dwws$winout <- ifelse(dwws$WCR %in% wellsinMN$WCR, "Within", "Outside")

wid <- dwws %>% group_by(BASIN) %>% summarise(AllWells=length(unique(WCR)),
                                              Outside=length(unique(WCR[winout=="Outside"])),
                                                                    Within=length(unique(WCR[winout=="Within"])))

wid$percout <- wid$Outside/wid$AllWells
wid$percin <- wid$Within/wid$AllWells

wid$`% of Drinking Water Wells Covered` <- ifelse(wid$percin <= .25, "0 - 25 %",
                                      ifelse(wid$percin >.25 & wid$percin <= .5, "26 - 50 %",
                                             ifelse(wid$percin >.5 & wid$percin <= .75, "51 - 75 %",
                                                    ifelse(wid$percin >.75 , "76 - 100 %", "error"))))

wdc <- wid %>% 
  tidyr::pivot_longer(., cols=c('% of DAC Area Covered', '% of Wells Partially or Fully Dry'), names_to='variable', values_to="value") %>% 
  select(NAMELSAD, BASIN, variable, value)

wdcc <- wdc %>% group_by(variable, value) %>% summarise(num_dac = length(unique(NAMELSAD)),
                                                        num_basin = length(unique(BASIN)))

wellcoverage <- ggplot()+
  geom_sf(data=dwws, aes(col=winout), size=.05)+
  #geom_sf(data=buffer, fill=NA, col="#032e6b", lwd=.5)+
  geom_sf(data=basins, fill=NA, col="gray20", lwd=.5)+
  scale_color_manual(values=c(orange2, blue1), name="Domestic and Public Supply \nWells Within or Outside \nof Monitoring Network")+
  theme_void(base_size = 16)+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.position = c(.80, 0.70))
wellcoverage

#ggsave(plot=wellcoverage, filename="Results/Images/FIGURE4_WellsMN.png", bg = "transparent", type = "cairo", device = "png")



#### FIGURE 4 - DAC - COVERAGE AND LOCATION OF IMPACTED WELLS ####
# plot DACs, 
dinmt <- st_intersection(dcs, buffer)
dinmt$int_area <- st_area(dinmt)
dinmt$percarea <- as.numeric(dinmt$int_area/dinmt$area)
dinmt$greater50 <- ifelse(dinmt$percarea >= .5, "Within", "Outside")
#dd <- dinmt %>% group_by(BASIN) %>% summarise(DACarea = sum(area), INTarea = sum(int_area), Perc=INTarea/DACarea)

dacscover <- dinmt %>% st_drop_geometry(.) %>% dplyr::select(GEOID, NAMELSAD, area, int_area, percarea, greater50)  %>% left_join(dcs,.)

dacscover$winout <- ifelse(is.na(dacscover$greater50)==TRUE, "Outside", dacscover$greater50)

# number of dacs and GSPs that fall outside of monitoring networks
dacscover %>% filter(winout=="DACs Outside Monitoring Network") %>% group_by(BASIN) %>% summarise(length(unique(NAMELSAD)))
length(unique(dacscover$GEOID[dacscover$winout=="DACs Outside Monitoring Network"]))


daccoverage <- ggplot()+
  geom_sf(data=basins, fill="#e3e6e8", col=NA)+
  geom_sf(data=buffer, fill="white", col="#032e6b", lwd=.3)+
  geom_sf(data=dacscover, aes(fill=winout, col=winout),lwd=.4)+
  geom_sf(data=basins, fill=NA, col="gray20", lwd=.3)+
  scale_fill_manual(values=c("#d1adab",blue2), name="DACs\nWithin or Outside of Monitoring Network", labels = function(x) str_wrap(x, width = 25))+
  scale_color_manual(values=c(red2,blue2))+
  theme_void(base_size = 14)+
  theme(legend.position = c(1, 0.70))+
  guides(col = "none")

daccoverage

#ggsave(plot=daccoverage, filename="Results/Images/FIGURE4_DACsMap.png", bg = "transparent", type = "cairo", device = "png", height=10,width=16)


# bar chart of % of impacted wells within dacs
mad_all <- wellanalysis(dwws) %>% st_drop_geometry(.) %>% dplyr::select(WCR,mt_at_dw,dry)
allwells <- left_join(dwws, mad_all)
wellsindacs <- st_intersection(allwells, dacscover) %>% select(WCR, type, BASIN, TCD, pump_loc, top, bot, mt_at_dw, dry, NAMELSAD, dcs_area, int_area, percarea)
wellsindacs$total <- ifelse(is.na(wellsindacs$dry)==FALSE, wellsindacs$dry,
                            ifelse(wellsindacs$percarea>=0.9999 & is.na(wellsindacs$dry)==TRUE, "active", "Out of Network"))

wellsindacs$total <- ifelse(is.na(wellsindacs$total)==TRUE, "Out of Network", wellsindacs$total)

###### disproportionate impact eval #####
length(unique(wellsindacs$WCR[wellsindacs$total=="Out of Network"]))

# filter wells not in dacs
k <- wellsindacs %>% filter(., !WCR %in% wellsindacs$WCR) %>% group_by(total) %>% summarise(n=length(unique(WCR)))
# (4027+5405+4189)/31488

# number and percent of wells within DACs that are not covered
wellsindacs %>% group_by(winout) %>% summarise(NumWells=length(unique(WCR)), NumDACs=length(unique(NAMELSAD)), NumGSP=length(unique(GSP.NAME)))

wellsindacs2 <- wellsindacs

wellsindacs <- st_intersection(wellsindacs, buffer)

View(wellsindacs %>% group_by(winout, winout.1, GSP.NAME) %>% summarise(NumWells=length(unique(WCR)), NumDACs=length(unique(NAMELSAD))))

# of the 3,484 wells within the 193 DACs in SJV, 508, or nearly 15% of domestic and public supply wells in 31 DACs fall outside of monitoring networks

da <- dacscover %>% group_by(BASIN) %>% summarise(ak1=length(unique(NAMELSAD))) %>% dplyr::select(BASIN, ak1) %>% st_drop_geometry(.)
DACS2 <- dacscover %>% group_by(BASIN) %>% summarise(ak2=length(unique(NAMELSAD)), PERC=as.numeric(sum(int_area, na.rm=TRUE)/sum(area, na.rm=TRUE))) %>% left_join(.,da)
#DACS2$percent <- 100*DACS2$ak2/DACS2$ak1
DACS2$BASIN <- substring(DACS2$BASIN, 9)

DACS2 %>% group_by(BASIN) %>%
  #mutate(within1 = sum(percent[winout=="Within"], na.rm=T)/sum(percent, na.rm=T))%>%
  ggplot(aes(x = reorder(BASIN,PERC, decreasing=T),
             y = PERC)) +
  geom_bar(stat = "identity") +
  #scale_y_continuous(labels = scales::percent)+
  coord_flip()


we <- wellsindacs %>% group_by(BASIN) %>% summarise(we1=length(unique(WCR))) %>% dplyr::select(BASIN, we1) %>% st_drop_geometry(.)
wellsindacs$status2 <- ifelse(wellsindacs$total == "Out of Network", "Wells Outside Monitoring Network",
                              ifelse(wellsindacs$total=="active", "Active Wells in Monitoring Network","Impacted Wells in Monitoring Network"))

kk <- wellsindacs %>% group_by(BASIN,status2) %>% summarise(we2=length(unique(WCR))) %>% left_join(.,we) %>% st_drop_geometry(.)
kk$percent_wells <- 100*kk$we2/kk$we1
kk$BASIN <- substring(kk$BASIN, 9)

###### DAC and WELL IMPACT FIG ######
DACS2$type <- "Coverage of DACs"
kk$type <- "Drinking Water Well Status Within DACs"
DACS2$Status <- "DAC Area Covered"
DACS2$percent <- 100*DACS2$PERC

IMPACTINDACS <- left_join(DACS2, kk, by="BASIN") %>% dplyr::select(BASIN,type.x, percent, type.y, status2, percent_wells)

d2 <- DACS2 %>% dplyr::select(BASIN, type, Status, percent) %>% st_drop_geometry(.)
w2 <- kk %>% dplyr::select(BASIN, type, Status=status2, percent=percent_wells)
o2 <- rbind(d2, w2)

o2$Levelz <- factor(o2$Status, levels=c("DAC Area Covered", "Active Wells in Monitoring Network", "Impacted Wells in Monitoring Network","Wells Outside Monitoring Network" ))

o2$basinlevles <- factor(o2$BASIN, levels=c(" SOLANO",
                                            " ENTERPRISE",
                                            " EASTERN SAN JOAQUIN",
                                            " RED BLUFF", 
                                            " VINA",
                                            " COLUSA",
                                            " WYANDOTTE CREEK",
                                            " KERN COUNTY",
                                            " MERCED",
                                            " DELTA-MENDOTA",
                                            " MADERA",
                                            " CHOWCHILLA",
                                            " NORTH AMERICAN",
                                            " YOLO",
                                            " TULE",
                                            " KINGS",
                                            " WESTSIDE",
                                            " MODESTO",
                                            " NORTH YUBA",
                                            " SOUTH YUBA",
                                            " BOWMAN",
                                            " SOUTH AMERICAN",
                                            " SUTTER",
                                            " CORNING",
                                            " KAWEAH",
                                            " BUTTE",
                                            " TULARE LAKE",
                                            " ANDERSON",
                                            " TURLOCK",
                                            " LOS MOLINOS",
                                            " EAST CONTRA COSTA",
                                            " COSUMNES"))

ggplot(data=o2, aes(x = basinlevles, y = percent, fill = Status), group = BASIN) +
  geom_bar(stat="identity", width=.5, position = "fill")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()

bardac <- ggplot() +
  geom_bar(data=o2, aes(x = basinlevles,
               y = percent,
               fill = Levelz),
           stat = "identity") +
  facet_grid(.~type)+
  #scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values=c("DAC Area Covered"= blue2, "Wells Outside Monitoring Network"="#b4bfb9", "Impacted Wells in Monitoring Network"=orange3, "Active Wells in Monitoring Network"=blue3), labels = function(x) str_wrap(x, width = 5))+
  labs(x ="Groundwater Basin", y = "", fill="Status")+
  coord_flip()+
  theme_bw(base_size=12)+
  theme(legend.position = 'bottom', 
        #legend.spacing.x = unit(1.0, 'cm'),
        #legend.text = element_text(margin = margin(b=10))
        ) +
  guides(fill = guide_legend(title = "Percent of:",
                             label.position = "bottom",
                             title.position = "top", title.hjust = .5)) 
bardac

l <- plot_grid(daccoverage, bardac, align="hv", labels = c("a.", "b."), label_size=20)
ggsave(plot=l, filename="Results/Images/FIGURE_DACS.png", bg = "white", type = "cairo", device = "png", height = 7, width = 14.5)

#ggsave(plot=bardac, filename="Results/Images/FIGURE5_DACsBar.png", bg = "transparent", type = "cairo", device = "png", height=10,width=16)


