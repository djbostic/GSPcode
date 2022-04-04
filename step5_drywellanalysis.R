# load packages
library(tidyverse) # general purpose data science toolkit
library(sp)        # spatial objects
library(raster)    # for raster objects
library(sf)
library(readr)

# set working directory
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("/Volumes/GoogleDrive/My Drive/Graduate School/GSP_Analy/Organized_CnD/")

# set coordinate reference system
merc <- crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
            +k=1.0 +units=m +nadgrids=@null +no_defs")

# load data
# GSP boundaries

gsp <- st_read("Boundaries/all_cobs/fixedgeometry_gspswithregion.shp") %>% filter(., region == "Coastal" & GSP_Name != "Las Posas")
gsps <- st_read("Boundaries/all_cobs/fixedgeometry_gspswithregion.shp")
gsp <- st_transform(gsp, merc) %>% filter(., region == "Coastal" & GSP_Name != "Las Posas")


# domestic wells
domwells_oswcr <- read_rds("Output/domesticwells_allcobs_90.rds")
dw_gsp <- st_transform(domwells_oswcr, merc)
dw_gsp <- st_intersection(dw_gsp, gsp)
#names(dw_gsp)[names(dw_gsp) == "WCRNumber"] <- "WCRNmbr"

# buffer outline
buffer_outline <- st_read("Output/buff_ts_coastal.shp")
buffer_outline <- as(buffer_outline, "sf")

# current gwl raster
cgwl_raster <- raster("Output/1819avginterp_allcobs.tif")
cgwl_raster <- projectRaster(cgwl_raster, crs=merc)

# minimum threshold raster
mt_raster <- read_rds("Output/minthreshinterpolation_allcobs.rds")
mt_raster <- mt_raster$Prediction
mt_raster <- projectRaster(mt_raster, crs=merc)

# compare how many wells are within GSP area and interpolation boundary
intial_data <- dw_gsp %>% group_by(GSP_Name, region) %>% summarise(Int=length(unique(WCRNmbr)), TCDa = mean(TtlCmpD, na.rm=TRUE))
inti_data <- st_set_geometry(intial_data, NULL)

buff <- st_intersection(dw_gsp, buffer_outline)
buff_data <- buff %>% group_by(GSP_Name, region) %>% summarise(Buffer=length(unique(WCRNmbr)), TCDb = mean(TtlCmpD, na.rm=TRUE))
buff_data <- st_set_geometry(buff_data, NULL)

notinboth <- dw_gsp[!(dw_gsp$WCRNmbr %in% buff$WCRNmbr),]

bi <- left_join(inti_data, buff_data, by="GSP_Name")
bi$percentcoverage <- bi$Buffer / bi$Int
bif <- filter(bi, percentcoverage < .60)

ggplot(bi, aes(x = reorder(GSP_Name, percentcoverage), y = percentcoverage, fill=factor(region.x))) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  labs(title = "Percent of domestic wells covered by GSP - Coastal Basins") + 
  xlab("GSP") + 
  ylab("Percent Covered") + 
  scale_fill_manual("Region", values = c("San Joaquin Valley" = "#11854D", "Coastal" = "#111F85", "Central" = "#853911")) + 
  theme_minimal()



bi$pc <- ifelse(bi$percentcoverage <=0.6, "< 0.6", 
                ifelse(bi$percentcoverage > 0.6 & bi$percentcoverage <= 0.8, "0.6 - 0.8", 
                       ifelse(bi$percentcoverage > 0.8 & bi$percentcoverage <= 0.9, "0.8 - 0.9", 
                              ifelse(bi$percentcoverage >= 0.9, "0.9 - 1", "NA"))))
percentmap <- left_join(gsp, bi, by="GSP_Name")
bi <- na.omit(bi)

library(RColorBrewer)
ggplot() + 
  geom_sf(data=gsps, col='transparent') +
  geom_sf(data=gsp, col="transparent") + 
  geom_sf(data=dw_gsp, cex=.1, col="cornflowerblue") +
  #geom_sf(data=tcddry, cex=.1, col="red")
  geom_sf(data=percentmap, aes(fill = pc)) +
  #geom_sf(data=notinboth, cex=.1, col = "red") + 
  theme_void(base_size = 16) +  
  scale_fill_manual(values = c("#bdc9e1", "#02818a","#f6eff7"), name= "Percent Covered") + 
  ggtitle("Percent of Domestic Wells Covered by GSP MTs")

,"#67a9cf",
# now that we know coverage, we move on to calculating what percent of domestic wells go dry 

#### Remove wells above current groundwater level ####
cgwl_at_dw <- raster::extract(cgwl_raster, dw_gsp) # intersect to get value of current water level at well points
cad <- cbind(dw_gsp, cgwl_at_dw) 
cad <- cad[!is.na(cad$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
cad$tcd_dry <- ifelse(cad$TtlCmpD <= cad$cgwl_at_dw, "Failing", "Active")
activedw <- cad[cad$tcd_dry == "Active", ]
activedw <- activedw[!is.na(activedw$WCRNmbr), ]

length(unique(activedw$WCRNmbr)) / length(unique(cad$WCRNmbr)) # percent of wells whose TCD is below current groundwater levels (useable wells)


#### assuming we dont have cgwl for coastal and non-cv central basins ####
activedw <- buff
#### dry well analysis ####
mt_at_dw <- raster::extract(mt_raster, activedw) # intersect to get value of current water level at well points
mad <- cbind(activedw, mt_at_dw) 
mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
mad$tcd_dry <- ifelse(mad$TtlCmpD <= mad$mt_at_dw, "Failing", "Active")
tcddry <- mad[mad$tcd_dry == "Failing", ]
tcddry <- tcddry[!is.na(tcddry$WCRNmbr),]

length(unique(tcddry$WCRNmbr)) / length(unique(mad$WCRNmbr)) # percent of wells that go dry in whole data set based on tcd

tcdbygsp <- mad %>% group_by(GSP_Name, tcd_dry) %>% summarise(num_tcd = length(unique(WCRNmbr))) %>% filter(., tcd_dry == "Failing")
st_geometry(tcdbygsp) <- NULL
t <- as.data.frame(tcdbygsp)

totalnums <- mad %>% group_by(GSP_Name) %>% summarise(numdw = length(unique(WCRNmbr))) 
st_geometry(totalnums) <- NULL
tn <- as.data.frame(totalnums)


all <- left_join(tn, t, by="GSP_Name")

all$tcd_perc <- all$num_tcd / all$numdw



a <- gsp$GSP_Name
levels(a)[1] <- "180/400 Ft Aquifer - Salinas"
b <- tcddry %>% group_by(GSP_Name) %>% summarise(tcdd=length(unique(WCRNmbr)))
b <- b$tcdd
b <- c(2,0,0,1)
c <- dw_gsp %>% group_by(GSP_Name) %>% summarise(tcd=length(unique(WCRNmbr)))
c <- c(14, 6, 110, 181)


hi <- data.frame("name"=a,"tcd_dry"=b,"tcd"=c)
hi$tcd_pd <- 100*(hi$tcd_dry / hi$tcd)

output <- read_csv("Output/drywellnumbers.csv")
output$tcd_pd <- 100*output$`Percent Dry Wells - TCD`

df.long<-melt(output[,c(1:3)])
df.long$variable <- as.character(df.long$variable)
df.long$variable[df.long$variable == "tcd_pd"] <- "MT below TCD"
names(df.long)[1] <- "name"

df.long <- filter(df.long, value > 0)

ggplot(df.long, aes(reorder(name,value),value,fill=Region))+
  coord_flip()+
  scale_y_continuous() +
  geom_bar(stat="identity",position="dodge", width = .5) +
  xlab("GSP Name") +
  ylab("Number of domestic wells that are likey dry") +
  scale_fill_manual("Region", values = c("San Joaquin Valley" = "#70AB53", "Coastal" = "#535EAB", "Central" = "#AB6353")) +
  guides(fill=guide_legend(title="")) +
  theme_minimal(base_size = 16)


#### raster map ####
# rasterize dry wells
r <- raster(buffer_outline)
res(r) <- 1610*6 # township-level grids (1,609 meters in a mile)
r <- rasterize(buffer_outline, r)
plot(r, col="grey90", lwd=10)
quads <- as(r, "SpatialPolygons")


tcddry <- tcddry[!is.na(tcddry$WCRNmbr),]
pld <- as_Spatial(tcddry)
pld <- rasterize(coordinates(pld), r, fun="count")
pld <- as(pld, "SpatialPolygonsDataFrame")
pld <- st_as_sf(pld)
pld$dw_pl <- ifelse(pld$layer <=10, "0 - 2", 
                    ifelse(pld$layer > 10 & pld$layer <= 20, "11 - 20", 
                           ifelse(pld$layer > 20 & pld$layer <= 40, "21 - 40", 
                                  ifelse(pld$layer > 40, "41 - 60", "NA"))))

plot(pld$geometry)
plot(quads, add=T, col="NA", lwd=.1)
plot(gsp$geometry, add=T, col="NA")

ints <- st_intersection(pld, st_make_valid(gsp))
ints$dw_pl <- factor(ints$dw_pl, levels=c("0 - 2", "11 - 20", "21 - 40", "41 - 60"))
ggplot() +
  geom_sf(data=gsps, fill="grey90", lwd=.01) +
  geom_sf(data=gsp, fill = "grey100", col="black", lwd=0.5) + 
  geom_sf(data=ints, aes(fill=dw_pl))  +
  geom_sf(data=tcddry, cex=.5, col="#9C331B") +
  scale_fill_brewer(palette = "OrRd") + 
  theme_void(base_size = 16)+ 
  labs(fill = "Number of dry wells\nper 36 square miles")


#### pump box plot ####
mtminmax <- mad %>% group_by(GSP_Name) %>% summarise(Max_MT = max(mt_at_dw, na.rm=TRUE), Min_MT = min(mt_at_dw, na.rm=TRUE)) %>% st_drop_geometry()

# select max pump location = 90th percentile
# account for 90% of wells
top_90p <- mad %>% count(GSP_Name) %>% filter(n > 300) %>% pull(GSP_Name)
  # sanity check
nrow(filter(mad, GSP_Name %in% top_90p)) / nrow(mad) 
# only show GSPs with 90% of all wells
dw_90 <- filter(mad, GSP_Name %in% top_90p)

# only show wells that are dry in both CIs
dwdri <- dw_90[dw_90$tcd_dry == "Failing", ]
dw_90 <- dw_90 %>% filter(TtlCmpD < 500)
mmsubset <- filter(mtminmax, GSP_Name %in% top_90p)

ggplot() + 
  geom_jitter(data = dw_90, aes(x=GSP_Name, y=TtlCmpD), col="grey80") +
  geom_jitter(data = dwdri, aes(x=GSP_Name, y=TtlCmpD), col="#C85124") + 
  scale_y_reverse() + 
  geom_rect(data = mmsubset, aes(xmin=GSP_Name, xmax=GSP_Name,ymin=Min_MT,ymax=Max_MT), fill="#900C3F", alpha=.6) + 
  geom_errorbar(data=mmsubset, aes(x= GSP_Name, ymin = Min_MT, ymax = Max_MT), col="grey20") +
  labs(x = "GSP Name", y = "Depth to Pump Location",title = "Distribution of Pump Locations by GSP", subtitle = "Red = Dry Wells \nGray = Active Wells \nBlack = Range of MTs") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
#### Visualizations ####


ggplot() + 
  geom_sf(data=gsps$geometry, fill = "#cacaca") + 
  geom_sf(data = dw_gsp$geometry, col = "#806775", cex=0.1) +
  geom_sf(data = mad$geometry, col = "#2E9AA9", cex = .2) + 
  ggtitle("OSWCR Wells within MT Interpolation")

ggplot() + 
  geom_sf(data=gsps$geometry, fill = "#cacaca") + 
  geom_sf(data = mad$geometry, col = "grey20", cex=0.1) +
  geom_sf(data = tcddry$geometry, col = "#C85124", cex = .2) +   
  labs(title = "40% of Domestic Wells Go Dry", subtitle = "Red = Dry DWs \nBlack = Active DWs")
