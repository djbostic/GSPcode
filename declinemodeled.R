# calculate decline from 2019 raster

steps <- seq(0, 400, 10)
cgwl <- read_rds("InterpolationGWLevels/cgwl_raster.rds") # 2019 current groundwater level raster

# subtract
cgwl_x <- list()
for (i in 1:length(steps)){
  cgwl_x[[i]] <- cgwl + steps[i]
}

# well data
wells <- dwws 
dom <- dws
psw <- psww

wellclean <- function(y){
  cgwl_at_dw <- raster::extract(cgwl, y) # intersect to get value of current water level at well points
  cad_dw <- cbind(y, cgwl_at_dw) 
  cad_dw <- cad_dw[!is.na(cad_dw$cgwl_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
  cad_dw$tcd_dry <- ifelse(cad_dw$TCD <= cad_dw$cgwl_at_dw, "Failing", "Active")
  activedw <- cad_dw[cad_dw$tcd_dry == "Active", ]
  activedw <- activedw[!is.na(activedw$WCR), ]
  print(length(unique(activedw$WCR)) / length(unique(cad_dw$WCR))) # percent of wells whose TCD is below current groundwater levels (useable wells)
  return(activedw)
}

pw <- wellclean(psw)
dw <- wellclean(dom)

wf_pw <- list()
for (i in 1:length(steps)){
  wf_pw[[i]] <- wellanalysis(x=cgwl_x[[i]], y=pw)
}

wf_pw2 <- data.frame(lapply(wf_pw, function(x) unlist(x)))
wf_pw2 <- t(wf_pw2)
rownames(wf_pw2) <- NULL
wf_pw3 <- as.data.frame(wf_pw2)
wf_pw3$FtDecline <- steps
wf_pw3$Type <- "Public Supply"


wf_dw <- list()
for (i in 1:length(steps)){
  wf_dw[[i]] <- wellanalysis(x=cgwl_x[[i]], y=dw)
}

wf_dw2 <- data.frame(lapply(wf_dw, function(x) unlist(x)))
wf_dw2 <- t(wf_dw2)
rownames(wf_dw2) <- NULL
wf_dw3 <- as.data.frame(wf_dw2)
wf_dw3$FtDecline <- steps
wf_dw3$Type <- "Domestic"

wf3 <- rbind(wf_dw3, wf_pw3)

wfmap <- ggplot(wf3, 
       aes(x = FtDecline, y = FailureRate)) +
  geom_line(aes(color = Type), lwd=2) +
  geom_point(aes(color = Type), cex=2) +
  labs(color="Well Type")+

  geom_point(aes(x=56, y=.15, alpha=""), col="#37693E", cex=2.2)+
  scale_alpha_manual(values=10)+
  labs(alpha = "Public Supply Well \nFailure Percentage at MT")+


  geom_point(aes(x=73, y=.18, fill=""), col="#6c4a85", cex=2.2)+
  scale_fill_manual(values=1)+
  labs(fill="Domestic Well \nFailure Percentage at MT")+
  
  labs(x = "Groundwater level decline (ft)", y = "Failure percentage") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values=c("#B5A5C2", "#9BB49F"))+
  theme_minimal(base_size = 12)+
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x =element_text(size=12), 
        axis.text.y =element_text(size=12),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "right",
        legend.box="vertical")
wfmap  

ggsave(plot=wfmap, filename="Results/Images/failureest.png", bg = "white", type = "cairo", device = "png", height = 7, width = 14.5)


wellanalysis <- function(x, y){
  # TCD
  mt_at_dw <- raster::extract(x, y) # intersect to get value of current water level at well points
  mad <- cbind(y, mt_at_dw) 
  mad <- mad[!is.na(mad$mt_at_dw), ] # remove wells where there is no current groundwater level value (wells are likely outside of interpolation area)
  mad$tcddry <- ifelse(mad$TCD <= mad$mt_at_dw, "Failing", "Active")
  #print(table(mad$tcddry))
  print(table(mad$tcddry)/length(unique(mad$WCR)))
  
  tcdry <- filter(mad, tcddry=="Failing")
  # bottom of well screen
  bots <- filter(mad, is.na(mad$bot)==FALSE & mad$bot > 0 & mad$tcddry == "Active")
  bots$botdry <- ifelse(bots$bot <= bots$mt_at_dw, "Failing", "Active")
  #print(table(bots$botdry))
  #print(table(bots$botdry)/length(unique(bots$WCR)))
  
  # all 
  mad$top <- ifelse(is.na(mad$top)==TRUE, 0, mad$top)
  mad$bot <- ifelse(is.na(mad$bot)==TRUE, 0, mad$bot)
  mad$pump_loc <- ifelse(is.na(mad$pump_loc)==TRUE, 0, mad$pump_loc)
  mad$dry <- ifelse((mad$mt_at_dw >= mad$TCD & mad$TCD > 0), "TCDdry",
                    ifelse(mad$top > 0 & mad$mt_at_dw >= mad$top & mad$mt_at_dw < mad$TCD, "topdry",
                           ifelse(mad$pump_loc > 0 & mad$mt_at_dw >= mad$pump_loc & mad$mt_at_dw < mad$TCD & mad$mt_at_dw < mad$top, "pumpdry", "active")))
  print(table(mad$dry)/length(unique(mad$WCR)))
  
  output <- data.frame(FailureCount=length(unique(tcdry$WCR)),
                       FailureRate = length(unique(tcdry$WCR))/length(unique(mad$WCR)))
  return(output)
}


