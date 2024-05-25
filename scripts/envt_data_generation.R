###-- Script for generating bio_4 variable values for 2nd iteration of Bayesian models

#####-- Environmental Association Models #####
setwd("~/Documents/PhD/Miami_Project/data-analysis/climate_data/")
library(raster)

##-- Data downloaded from WORLDCLIM and PALEOCLIM websites (details in manuscript)
preclim2 <- raster("LH_v1_2_5m/bio_2.tif")
preclim2 <- rasterToPoints(preclim2)
preclim2 <- data.frame(preclim2)
preclim4 <- raster("LH_v1_2_5m/bio_4.tif")
preclim4 <- rasterToPoints(preclim4)
preclim4 <- data.frame(preclim4)
preclim7 <- raster("LH_v1_2_5m/bio_7.tif")
preclim7 <- rasterToPoints(preclim7)
preclim7 <- data.frame(preclim7)

sflclim2 <- raster("CHELSA_cur_V1_2B_r2_5m_79-13/2_5min/bio_2.tif")
sflclim2 <- rasterToPoints(sflclim2)
sflclim2 <- data.frame(sflclim2)
sflclim4 <- raster("CHELSA_cur_V1_2B_r2_5m_79-13/2_5min/bio_4.tif")
sflclim4 <- rasterToPoints(sflclim4)
sflclim4 <- data.frame(sflclim4)
sflclim7 <- raster("CHELSA_cur_V1_2B_r2_5m_79-13/2_5min/bio_7.tif")
sflclim7 <- rasterToPoints(sflclim7)
sflclim7 <- data.frame(sflclim7)

#see Paleoclim and Chelsa websites for scaling factors to interpret value

###-- trim down species ranges

## import occurrence points to delimit ranges; derived from GBIF datasets as detailed in manuscript
caropoints <- read.csv("geo_records/carolinensis_gbif_3-21-23.csv", sep="\t")
chloropoints <- read.csv("geo_records/chlorocyanus_gbif_3-21-23.csv", sep="\t")
cristapoints <- read.csv("geo_records/cristatellus_gbif_3-21-23.csv", sep="\t")
cybpoints <- read.csv("geo_records/cybotes_gbif_3-21-23.csv", sep="\t")
distpoints <- read.csv("geo_records/distichus_gbif_3-21-23.csv", sep="\t")
eqpoints <- read.csv("geo_records/equestris_gbif_3-21-23.csv", sep="\t")
sagpoints <- read.csv("geo_records/sagrei_gbif_3-21-23.csv", sep="\t")

caropoints <- caropoints[,c("decimalLongitude","decimalLatitude")]
chloropoints <- chloropoints[,c("decimalLongitude","decimalLatitude")]
cristapoints <- cristapoints[,c("decimalLongitude","decimalLatitude")]
cybpoints <- cybpoints[,c("decimalLongitude","decimalLatitude")]
distpoints <- distpoints[,c("decimalLongitude","decimalLatitude")]
eqpoints <- eqpoints[,c("decimalLongitude","decimalLatitude")]
sagpoints <- sagpoints[,c("decimalLongitude","decimalLatitude")]

caropoints <- caropoints[!is.na(caropoints$decimalLongitude),]
caropoints <- caropoints[!is.na(caropoints$decimalLatitude),]
chloropoints <- chloropoints[!is.na(chloropoints$decimalLongitude),]
chloropoints <- chloropoints[!is.na(chloropoints$decimalLatitude),]
cristapoints <- cristapoints[!is.na(cristapoints$decimalLongitude),]
cristapoints <- cristapoints[!is.na(cristapoints$decimalLatitude),]
cybpoints <- cybpoints[!is.na(cybpoints$decimalLongitude),]
cybpoints <- cybpoints[!is.na(cybpoints$decimalLatitude),]
distpoints <- distpoints[!is.na(distpoints$decimalLongitude),]
distpoints <- distpoints[!is.na(distpoints$decimalLatitude),]
eqpoints <- eqpoints[!is.na(eqpoints$decimalLongitude),]
eqpoints <- eqpoints[!is.na(eqpoints$decimalLatitude),]
sagpoints <- sagpoints[!is.na(sagpoints$decimalLongitude),]
sagpoints <- sagpoints[!is.na(sagpoints$decimalLatitude),]

#carolinensis
#introduced and native the same
carosquares <- subset(caropoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                        decimalLatitude<27 & decimalLatitude>25.1)
max(carosquares$decimalLongitude) #-80.00203
min(carosquares$decimalLongitude) #-82.38245
max(carosquares$decimalLatitude) #26.99368
min(carosquares$decimalLatitude) #25.12291

#presentclim clim
sflclim2car <- subset(sflclim2, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
sflclim4car <- subset(sflclim4, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
sflclim7car <- subset(sflclim7, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
sflclimcar <- data.frame(cbind(sflclim2car$bio_2, sflclim4car$bio_4,sflclim7car$bio_7))
names(sflclimcar) <- c("bio_2", "bio_4", "bio_7")
sflcaro_corrmat <- cor(sflclimcar)
ggcorrplot(sflcaro_corrmat)
boxplot(sflclimcar$bio_4)
sflpcacar <- princomp(sflclimcar)
summary(sflpcacar)
screeplot(sflpcacar)
fviz_contrib(sflpcacar, choice = "var", axes =1, top = 10) 
invcarmean <- mean(sflclim4car$bio_4) #3236.783
invcarvar <- var(sflclim4car$bio_4) #62408.28
#t.test between SFL and preclim for bio 4 reveal significantly different

#pre-invasion clim
#1773 at end, no NAs
preclim2caro <- subset(preclim2, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
preclim4caro <- subset(preclim4, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
preclim7caro <- subset(preclim7, x<(-80.00203) & x>-82.38245 & y<26.99368 & y>25.12291)
preclimscaro <- data.frame(cbind(preclim2caro$bio_2, preclim4caro$bio_4, preclim7caro$bio_7))
names(preclimscaro) <- c("bio_2", "bio_4", "bio_7")
precaro_corrmat <- cor(preclimscaro)
ggcorrplot(precaro_corrmat)
prepcacaro<- princomp(preclimscaro)
summary(prepcacaro) #basically everything is on PC1
screeplot(prepcacaro, addlabels = TRUE) #screeplot
fviz_contrib(prepcacaro, choice = "var", axes =1, top = 10) 
natcarmean <- mean(preclim4caro$bio_4) #3319.522
natcarvar <- var(preclim4caro$bio_4) #69857.59
  
  
#chlorocyanus
#native
natchloro <- subset(chloropoints, decimalLongitude<(-68.32) & decimalLongitude>-74.5 & 
                   decimalLatitude<20.1 & decimalLatitude>17.46)
max(natchloro$decimalLongitude) #-68.37216
min(natchloro$decimalLongitude) #-74.42285
max(natchloro$decimalLatitude) #20.04534
min(natchloro$decimalLatitude) #17.76715

# 4133 squares
natclim2chloro <- subset(preclim2, x<(-68.32) & x>-74.5 & y<20.1 & y>17.46)
natclim4chloro <- subset(preclim4, x<(-68.32) & x>-74.5 & y<20.1 & y>17.46)
natclim7chloro <- subset(preclim7, x<(-68.32) & x>-74.5 & y<20.1 & y>17.46)
natclimschloro <- data.frame(cbind(natclim2chloro$bio_2, natclim4chloro$bio_4, natclim7chloro$bio_7))
names(natclimschloro) <- c("bio_2", "bio_4", "bio_7")
natchloro_corrmat <- cor(natclimschloro)
ggcorrplot(natclimschloro)
natpcachloro<- princomp(natclimschloro)
summary(natpcachloro) #basically everything is on PC1
screeplot(natpcachloro, addlabels = TRUE) #screeplot
fviz_contrib(natpcachloro, choice = "var", axes =1, top = 10) 
natchloromean <- mean(natclim4chloro$bio_4) #1325.395
natchlorovar <- var(natclim4chloro$bio_4) #26279.88

#introduced
invchloro <- subset(chloropoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                      decimalLatitude<27 & decimalLatitude>25.1)
max(invchloro$decimalLongitude) #-80.06832
min(invchloro$decimalLongitude) #-80.39812
max(invchloro$decimalLatitude) #26.66959
min(invchloro$decimalLatitude) #26.01847

#128 squares
invclim2chloro <- subset(sflclim2, x<(-80.06832) & x>-80.39812 & y<26.66959 & y>26.01847)
invclim4chloro <- subset(sflclim4, x<(-80.06832) & x>-80.39812 & y<26.66959 & y>26.01847)
invclim7chloro <- subset(sflclim7, x<(-80.06832) & x>-80.39812 & y<26.66959 & y>26.01847)
invclimschloro <- data.frame(cbind(invclim2chloro$bio_2, invclim4chloro$bio_4, invclim7chloro$bio_7))
names(invclimschloro) <- c("bio_2", "bio_4", "bio_7")
invchloro_corrmat <- cor(invclimschloro)
ggcorrplot(invchloro_corrmat)
prepcachloro<- princomp(invclimschloro)
summary(prepcachloro) #basically everything is on PC1
screeplot(prepcachloro, addlabels = TRUE) #screeplot
fviz_contrib(prepcachloro, choice = "var", axes =1, top = 10) 
invchloromean <- mean(invclim4chloro$bio_4) #3059.977
invchlorovar <- var(invclim4chloro$bio_4) #4048.637
                           

#cristatellus
#native
natcrista <- subset(cristapoints, decimalLongitude<(-64.26) & decimalLongitude>(-67.27) & 
                  decimalLatitude<18.75 & decimalLatitude>17.92)
max(natcrista$decimalLongitude) #-64.27527
min(natcrista$decimalLongitude) #-67.26422
max(natcrista$decimalLatitude) #18.73
min(natcrista$decimalLatitude) #17.93191

# 614 squares
natclim2crista <- subset(preclim2, x<(-64.27527) & x>-67.26422 & y<18.73 & y>17.93191)
natclim4crista <- subset(preclim4, x<(-64.27527) & x>-67.26422 & y<18.73 & y>17.93191)
natclim7crista <- subset(preclim7, x<(-64.27527) & x>-67.26422 & y<18.73 & y>17.93191)
natclimscrista <- data.frame(cbind(natclim2crista$bio_2, natclim4crista$bio_4, natclim7crista$bio_7))
names(natclimscrista) <- c("bio_2", "bio_4", "bio_7")
natcrista_corrmat <- cor(natclimscrista)
ggcorrplot(natcrista_corrmat)
natpcacrista<- princomp(natclimscrista)
summary(natpcacrista) #basically everything is on PC1
screeplot(natpcacrista, addlabels = TRUE) #screeplot
fviz_contrib(natpcacrista, choice = "var", axes =1, top = 10) 
natcristamean <- mean(natclim4crista$bio_4) #1148.06
natcristavar <- var(natclim4crista$bio_4) #1064.602

#introduced
invcrista <- subset(cristapoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                      decimalLatitude<27 & decimalLatitude>25.1)
max(invcrista$decimalLongitude) #-80.0002
min(invcrista$decimalLongitude) #-82.19253
max(invcrista$decimalLatitude) #26.97939
min(invcrista$decimalLatitude) #25.12523
#1746 squares
invclim2crista <- subset(sflclim2, x<(-80.0002) & x>-82.19253 & y<26.97939 & y>25.12523)
invclim4crista <- subset(sflclim4, x<(-80.0002) & x>-82.19253 & y<26.97939 & y>25.12523)
invclim7crista <- subset(sflclim7, x<(-80.0002) & x>-82.19253 & y<26.97939 & y>25.12523)
invclimscrista <- data.frame(cbind(invclim2crista$bio_2, invclim4crista$bio_4, invclim7crista$bio_7))
names(invclimscrista) <- c("bio_2", "bio_4", "bio_7")
invcrista_corrmat <- cor(invclimscrista)
ggcorrplot(invcrista_corrmat)
prepcacrista<- princomp(invclimscrista)
summary(prepcacrista) #basically everything is on PC1
screeplot(prepcacrista, addlabels = TRUE) #screeplot
fviz_contrib(prepcacrista, choice = "var", axes =1, top = 10) 
invcristamean <- mean(invclim4crista$bio_4) #3229.372
invcristavar <- var(invclim4crista$bio_4) #5971.38


#cybotes
#native
natcyb <- subset(cybpoints, decimalLongitude<(-68.32) & decimalLongitude>-74.5 & 
                   decimalLatitude<20.1 & decimalLatitude>17.46)
max(natcyb$decimalLongitude) #-68.36886
min(natcyb$decimalLongitude) #-74.3986
max(natcyb$decimalLatitude) #20.04534
min(natcyb$decimalLatitude) #17.5733

#4098 squares
natclim2cyb <- subset(preclim2, x<(-68.36886) & x>-74.3986 & y<20.04534 & y>17.5733)
natclim4cyb <- subset(preclim4, x<(-68.36886) & x>-74.3986 & y<20.04534 & y>17.5733)
natclim7cyb <- subset(preclim7, x<(-68.36886) & x>-74.3986 & y<20.04534 & y>17.5733)
natclimscyb <- data.frame(cbind(natclim2cyb$bio_2, natclim4cyb$bio_4, natclim7cyb$bio_7))
names(natclimscyb) <- c("bio_2", "bio_4", "bio_7")
natcyb_corrmat <- cor(natclimscyb)
ggcorrplot(natcyb_corrmat)
natpcacyb<- princomp(natclimscyb)
summary(natpcacyb) #basically everything is on PC1
screeplot(natpcacyb, addlabels = TRUE) #screeplot
fviz_contrib(natpcacyb, choice = "var", axes =1, top = 10)#bio7 slightly more important here
natcybmean <- mean(natclim4cyb$bio_4) #1326.114
natcybvar <- var(natclim4cyb$bio_4) #26209.41


#introduced
invcyb <- subset(cybpoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                   decimalLatitude<27 & decimalLatitude>25.1)
max(invcyb$decimalLongitude) #-80.1772
min(invcyb$decimalLongitude) #-80.77662
max(invcyb$decimalLatitude) #26.9975
min(invcyb$decimalLatitude) #25.96269
#375 squares
invclim2cyb <- subset(sflclim2, x<(-80.1772) & x>-80.77662 & y<26.9975 & y>25.96269)
invclim4cyb <- subset(sflclim4, x<(-80.1772) & x>-80.77662 & y<26.9975 & y>25.96269)
invclim7cyb <- subset(sflclim7, x<(-80.1772) & x>-80.77662 & y<26.9975 & y>25.96269)
invclimscyb <- data.frame(cbind(invclim2cyb$bio_2, invclim4cyb$bio_4, invclim7cyb$bio_7))
names(invclimscyb) <- c("bio_2", "bio_4", "bio_7")
invcyb_corrmat <- cor(invclimscyb)
ggcorrplot(invcyb_corrmat)
prepcacyb<- princomp(invclimscyb)
summary(prepcacyb) #basically everything is on PC1
screeplot(prepcacyb, addlabels = TRUE) #screeplot
fviz_contrib(prepcacyb, choice = "var", axes =1, top = 10) 
invcybmean <- mean(invclim4cyb$bio_4) #3181.285
invcybvar <- var(invclim4cyb$bio_4) #11166.35

#distichus
#native
natdist <- subset(distpoints, decimalLongitude<(-68.32) & decimalLongitude>-74.5 & 
                 decimalLatitude<20.1 & decimalLatitude>17.46)
max(natdist$decimalLongitude) #-68.34635
min(natdist$decimalLongitude) #-74.42285
max(natdist$decimalLatitude) #20.04534
min(natdist$decimalLatitude) #17.5994

#4102 squares
natclim2dist <- subset(preclim2, x<(-68.34635) & x>-74.42285 & y<20.04534 & y>17.5994)
natclim4dist <- subset(preclim4, x<(-68.34635) & x>-74.42285 & y<20.04534 & y>17.5994)
natclim7dist <- subset(preclim7, x<(-68.34635) & x>-74.42285 & y<20.04534 & y>17.5994)
natclimsdist <- data.frame(cbind(natclim2dist$bio_2, natclim4dist$bio_4, natclim7dist$bio_7))
names(natclimsdist) <- c("bio_2", "bio_4", "bio_7")
natdist_corrmat <- cor(natclimsdist)
ggcorrplot(natdist_corrmat)
natpcadist<- princomp(natclimsdist)
summary(natpcadist) #basically everything is on PC1
screeplot(natpcadist, addlabels = TRUE) #screeplot
fviz_contrib(natpcadist, choice = "var", axes =1, top = 10) 
natdistmean <- mean(natclim4dist$bio_4) #1325.938
natdistvar <- var(natclim4dist$bio_4) #26215.46

#introduced -- double check
invdist <- subset(distpoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                    decimalLatitude<27 & decimalLatitude>25.1)
max(invdist$decimalLongitude) #-80.09112
min(invdist$decimalLongitude) #-80.45226
max(invdist$decimalLatitude) #26.6851
min(invdist$decimalLatitude) #25.31559
#257 squares
invclim2dist <- subset(sflclim2, x<(-80.09112) & x>-80.45226 & y<26.6851 & y>25.31559)
invclim4dist <- subset(sflclim4, x<(-80.09112) & x>-80.45226 & y<26.6851 & y>25.31559)
invclim7dist <- subset(sflclim7, x<(-80.09112) & x>-80.45226 & y<26.6851 & y>25.31559)
invclimsdist <- data.frame(cbind(invclim2dist$bio_2, invclim4dist$bio_4, invclim7dist$bio_7))
names(invclimsdist) <- c("bio_2", "bio_4", "bio_7")
invdist_corrmat <- cor(invclimsdist)
ggcorrplot(invdist_corrmat)
prepcadist<- princomp(invclimsdist)
summary(prepcadist) #basically everything is on PC1
screeplot(prepcadist, addlabels = TRUE) #screeplot
fviz_contrib(prepcadist, choice = "var", axes =1, top = 10) 
invdistmean <- mean(invclim4dist$bio_4) #2987.319
invdistvar <- var(invclim4dist$bio_4) #12642.18

#equestris
#native
nateq <- subset(eqpoints, decimalLongitude<(-74.12) & decimalLongitude>-85 & 
                  decimalLatitude<23.28 & decimalLatitude>19.82)
max(nateq$decimalLongitude) #-74.66794
min(nateq$decimalLongitude) #-83.7167
max(nateq$decimalLatitude) #23.19968
min(nateq$decimalLatitude) #19.9333

#6378 squares
natclim2eq <- subset(preclim2, x<(-74.66794) & x>-83.7167 & y<23.19968 & y>19.9333)
natclim4eq <- subset(preclim4, x<(-74.66794) & x>-83.7167 & y<23.19968 & y>19.9333)
natclim7eq <- subset(preclim7, x<(-74.66794) & x>-83.7167 & y<23.19968 & y>19.9333)
#second step of trimming to remove piece of bahamas in upper right
bahamabye <- subset(natclim4eq, x>-78 & y>21.9)
natclim2eq <- natclim2eq[!(row.names(natclim2eq) %in% rownames(bahamabye)),]
natclim4eq <- natclim4eq[!(row.names(natclim4eq) %in% rownames(bahamabye)),]
natclim7eq <- natclim7eq[!(row.names(natclim7eq) %in% rownames(bahamabye)),]
natclimseq <- data.frame(cbind(natclim2eq$bio_2, natclim4eq$bio_4, natclim7eq$bio_7))
names(natclimseq) <- c("bio_2", "bio_4", "bio_7")
nateq_corrmat <- cor(natclimseq)
ggcorrplot(nateq_corrmat)
natpcaeq<- princomp(natclimseq)
summary(natpcaeq) #basically everything is on PC1
screeplot(natpcaeq, addlabels = TRUE) #screeplot
fviz_contrib(natpcaeq, choice = "var", axes =1, top = 10) 
nateqmean <- mean(natclim4eq$bio_4) #1773.516
nateqvar <- var(natclim4eq$bio_4) #67058.43

#introduced
inveq <- subset(eqpoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                  decimalLatitude<27 & decimalLatitude>25.1)
max(inveq$decimalLongitude) #-80.00553
min(inveq$decimalLongitude) #-82.36154
max(inveq$decimalLatitude) #26.98279
min(inveq$decimalLatitude) #25.1074
#1773 squares
invclim2eq <- subset(sflclim2, x<(-80.00553) & x>-82.36154 & y<26.98279 & y>25.1074)
invclim4eq <- subset(sflclim4, x<(-80.00553) & x>-82.36154 & y<26.98279 & y>25.1074)
invclim7eq <- subset(sflclim7, x<(-80.00553) & x>-82.36154 & y<26.98279 & y>25.1074)
invclimseq <- data.frame(cbind(invclim2eq$bio_2, invclim4eq$bio_4, invclim7eq$bio_7))
names(invclimseq) <- c("bio_2", "bio_4", "bio_7")
inveq_corrmat <- cor(invclimseq)
ggcorrplot(inveq_corrmat)
prepcaeq<- princomp(invclimseq)
summary(prepcaeq) #basically everything is on PC1
screeplot(prepcaeq, addlabels = TRUE) #screeplot
fviz_contrib(prepcaeq, choice = "var", axes =1, top = 10) 
inveqmean <- mean(invclim4eq$bio_4) #3236.783
inveqvar <- var(invclim4eq$bio_4) #62408.28

ggplot(natclim4eq, aes(x = x, y = y, colour = bio_4)) +
  geom_point()

#sagrei
#native
natsag1 <- subset(sagpoints, decimalLongitude<(-74.12) & decimalLongitude>-85 & 
                    decimalLatitude<23.28 & decimalLatitude>19.82) #Cuba
natsag2 <- subset(sagpoints, decimalLongitude<(-73.5) & decimalLongitude>-79.31 & 
                    decimalLatitude<27.27 & decimalLatitude>22.12) #Bahamas
natsag <- rbind(natsag1, natsag2) #doesn't matter if some points are double counted bc this is just used to set bounds for climate polygons
max(natsag$decimalLongitude) #-74.16169
min(natsag$decimalLongitude) #-84.9508
max(natsag$decimalLatitude) #27.12381
min(natsag$decimalLatitude) #19.83928

#if I draw the initial polygon large enough to capture both then I should be able to just cut Florida and get it

natclim2sag <- subset(preclim2, x<(-74.16169) & x>-84.9508 & y<27.12381 & y>19.83928)
natclim4sag <- subset(preclim4, x<(-74.16169) & x>-84.9508 & y<27.12381 & y>19.83928)
natclim7sag <- subset(preclim7, x<(-74.16169) & x>-84.9508 & y<27.12381 & y>19.83928)
flcut <- subset(natclim4sag, x<(-79.5) & y>24.2)
natclim2sag <- natclim2sag[!(row.names(natclim2sag) %in% rownames(flcut)),]
natclim4sag <- natclim4sag[!(row.names(natclim4sag) %in% rownames(flcut)),]
natclim7sag <- natclim7sag[!(row.names(natclim7sag) %in% rownames(flcut)),]
natclimssag <- data.frame(cbind(natclim2sag$bio_2, natclim4sag$bio_4, natclim7sag$bio_7))
names(natclimssag) <- c("bio_2", "bio_4", "bio_7")
natsag_corrmat <- cor(natclimssag)
ggcorrplot(natsag_corrmat)
natpcasag<- princomp(natclimssag)
summary(natpcasag) #basically everything is on PC1
screeplot(natpcasag, addlabels = TRUE) #screeplot
fviz_contrib(natpcasag, choice = "var", axes =1, top = 10) 
natsagmean <- mean(natclim4sag$bio_4) #1862.801
natsagvar <- var(natclim4sag$bio_4) #108681.4


#introduced
invsag <- subset(sagpoints, decimalLongitude<(-80) & decimalLongitude>-82.45 & 
                   decimalLatitude<27 & decimalLatitude>25.1)
max(invsag$decimalLongitude) #-80.00339
min(invsag$decimalLongitude) #-82.40524
max(invsag$decimalLatitude) #26.99961
min(invsag$decimalLatitude) #25.10251
#1793 squares
invclim2sag <- subset(sflclim2, x<(-80.00339) & x>-82.40524 & y<26.99961 & y>25.10251)
invclim4sag <- subset(sflclim4, x<(-80.00339) & x>-82.40524 & y<26.99961 & y>25.10251)
invclim7sag <- subset(sflclim7, x<(-80.00339) & x>-82.40524 & y<26.99961 & y>25.10251)
invclimssag <- data.frame(cbind(invclim2sag$bio_2, invclim4sag$bio_4, invclim7sag$bio_7))
names(invclimssag) <- c("bio_2", "bio_4", "bio_7")
invsag_corrmat <- cor(invclimssag)
ggcorrplot(invsag_corrmat)
prepcasag<- princomp(invclimssag)
summary(prepcasag) #basically everything is on PC1
screeplot(prepcasag, addlabels = TRUE) #screeplot
fviz_contrib(prepcasag, choice = "var", axes =1, top = 10) 
invsagmean <- mean(invclim4sag$bio_4) #3232.857
invsagvar <- var(invclim4sag$bio_4) #64123.73

#process and export data
#meansand var for all spp pre and post - save
sppnames_alph <- c("carolinensis", "chlorocyanus","cristatellus","cybotes","distichus",
                   "equestris","sagrei")
premeans <- c(natcarmean, natchloromean, natcristamean, natcybmean, natdistmean,
              nateqmean, natsagmean)
invmeans <- c(invcarmean, invchloromean, invcristamean, invcybmean, invdistmean,
              inveqmean, invsagmean)
prevars <- c(natcarvar, natchlorovar, natcristavar, natcybvar, natdistvar,
             nateqvar, natsagvar)
invvars <- c(invcarvar, invchlorovar, invcristavar, invcybvar, invdistvar,
             inveqvar, invsagvar)
envts <- data.frame(cbind(sppnames_alph, premeans, prevars, invmeans, invvars))

save(envts, file="Miami_envtavgs.rda")
#write.csv(envts, file="Miami_envtavgs.csv")


#distributions for all spp pre and post - save
save(sflclim4car, natclim4chloro, natclim4crista, natclim4cyb,
              natclim4dist, natclim4eq, natclim4sag, preclim4caro, invclim4chloro,
     invclim4crista, invclim4cyb, invclim4dist, invclim4eq, invclim4sag, file="bio_4s.rda")

