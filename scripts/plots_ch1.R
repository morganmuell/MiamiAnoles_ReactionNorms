##### Plots for Miami Paper #####

setwd("~/Documents/PhD/Miami_Project/data-analysis")

library(RColorBrewer)
library(ggplot2)
library(emmeans)
library(MCMCglmm)
library(nadiv)
library(reshape)
library(ggridges)
library(viridis)
library(raster)
library(wesanderson)

rm(list=ls())

###-- Import data
#mializ <- read.csv(file="R_imports/miamidat_full.csv", header=TRUE, sep=",")

load("miami-lizards.rda") #***
load("bio_4s.rda") ###load for heat map fig
load("model2df.rda") #***
load("model1_results.rda") #***
load("inputdata_trimmed.rda") #***

sppnames <- c("cristatellus", "sagrei", "carolinensis", "chlorocyanus",
              "distichus", "cybotes", "equestris")
phenotypes <- c("WaterUptake", "Inc", "SVL", "SprintSpeed", "Endurance")

colsinter <- viridis(n=5)
cols <- palette(brewer.pal(n=7, name="Set2")) #intraspecific

### Heat Map Figure ###
### Make Heatmap Figure for Miami Paper ###

###-- Import data
preclim4 <- raster("LH_v1_2_5m/bio_4.tif")
sflclim4 <- raster("CHELSA_cur_V1_2B_r2_5m_79-13/2_5min/bio_4.tif")

#clip to: y>24.98, y<27.2, x>-82.58, x<-79.97
lims <- extent(-86, -65, 16.7, 27.2) #this works
presentclims <- crop(sflclim4, lims)
pastclims <- crop(preclim4, lims)

par(mfrow=c(2,1))
par(bty="n")
plot(presentclims, axes=F, bty="n", col=viridis(200), legend.args=list(cex=0.5))
par(bty="n")
plot(pastclims, axes=F, col=viridis(200))
dev.off()

#tiff(filename="presentclims.tiff", units="in", width=5, height=4, res=700)
svg(filename="presentclims.svg", width=5, height=4)
par(bty="n")
plot(presentclims, axes=F, col=viridis(200))
dev.off()

par(mfrow=c(2,1))
par(bty="n")
plot(presentclims, axes=F, bty="n", col=viridis(200), legend.args=list(cex=0.5))
par(bty="n")
plot(pastclims, axes=F, col=viridis(200))
dev.off()

#tiff(filename="presentclims.tiff", units="in", width=5, height=4, res=700)
svg(filename="pastclims.svg", width=5, height=4)
par(bty="n")
plot(pastclims, axes=F, col=viridis(200))
dev.off()


##### Interspecific figure NEW 3/19/24 ####
prior1 <- list(R = list(V = 1e-12, nu = -2), #<-- Non-informative improper: *marginal* posterior equal to REML estimate
               G = list(G1 = list(V = diag(2)*0.002, nu = 3, 
                                  alpha.mu = c(0,0), alpha.V = diag(2)*10000)))
#new plan - make posterior density plots with lines for 95% credible intervals for each species sensu fig. 4 Gross et al. 2024 https://onlinelibrary-wiley-com.spot.lib.auburn.edu/doi/epdf/10.1111/brv.13025
#pin in this for now; good idea but is log-scaled and want to confirm whether that will change decision making

# Water Uptake
#would be more efficient to load it in
set.seed(330)
bwumodel1 <- MCMCglmm(logWU ~ 1 + EggMass + nTreatment + PercDev + Species
                      + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = wu, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = TRUE, prior = prior1)

#I think, Species is y horizontals and effect is x, so, Species is Class in example
#then I need to make 5 of these plots, one for each phenotype, and line up side-by-side in Inkscape

# Pull marginals posterior modes for overall model effect
wutr.out <- emtrends(object = bwumodel1, pairwise ~ 1, var="nTreatment", data = wu) #for BOTTOMMOST row; all species + overall
spp.wu.1 <- emtrends(object = bwumodel1, pairwise ~ Species, var="nTreatment", data=wu)
#wutr.df.spp <- data.frame(hpd.summary(spp.wu.1, point.est = posterior.mode))
                          
# Bind columns of posterior estimates, melt into proper format for plotting
spp.wu.mc <- as.data.frame(as.mcmc.emmGrid(spp.wu.1$emtrends))
wutr.df.mc <- as.data.frame(as.mcmc.emmGrid(wutr.out$emtrends))
colnames(spp.wu.mc)<-gsub("Species","",colnames(spp.wu.mc))
colnames(wutr.df.mc)<-gsub("1 overall","Overall",colnames(wutr.df.mc))
wu.em <- cbind(spp.wu.mc, wutr.df.mc)
melt.wu.em <- melt(wu.em)

## Generate plot
svg(file = "R_exports/m1_trt_eff/wu_trteff.svg")
ggplot(melt.wu.em, aes(x = value, y = variable)) +
  coord_cartesian(ylim = c(1.4, 8.2)) +
  xlim(-1.1, 1.1) +
  scale_y_discrete(limits = unique(rev(melt.wu.em$variable))) +
  geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
             linewidth = 1.75, color = "grey10") +
  geom_density_ridges2(rel_min_height = 5e-11, scale = 0.95,
                       stat = "density_ridges_HPDCrI", 
                       quantile_lines = TRUE,
                       quantiles = c(0.95),
                       quantile_fun = HPD_fun,
                       fill = colsinter[1],
                       alpha = 0.6) +  #<-- transparency (0 more -- 1 less/solid)
  theme_classic(base_size=15) +
  ylab("") + xlab("Treatment Effect Size (95% Credible Interval)")
dev.off()

#incubation period
set.seed(330)
bincmodel1 <- MCMCglmm(logInc ~ 1 + Species + EggMass + nTreatment + Species*nTreatment, 
                       random = ~us(1 + nTreatment):Cage, pr = TRUE,
                       data = incp, nitt = 275000, thin = 50, burnin = 25000, 
                       verbose = TRUE, prior = prior1)

# Pull marginals posterior modes for overall model and Species fixed effect
inctr.out <- emtrends(object = bincmodel1, pairwise ~ 1, var="nTreatment", data = incp) #for BOTTOMMOST row; all species + overall
spp.inc.1 <- emtrends(object = bincmodel1, pairwise ~ Species, var="nTreatment", data=incp)

# Bind columns of posterior estimates, melt into proper format for plotting
spp.inc.mc <- as.data.frame(as.mcmc.emmGrid(spp.inc.1$emtrends))
inctr.df.mc <- as.data.frame(as.mcmc.emmGrid(inctr.out$emtrends))
colnames(spp.inc.mc)<-gsub("Species","",colnames(spp.inc.mc))
colnames(inctr.df.mc)<-gsub("1 overall","Overall",colnames(inctr.df.mc))
inc.em <- cbind(spp.inc.mc, inctr.df.mc)
melt.inc.em <- melt(inc.em)

#melt.tmp <- melt(spp.inc.mc)

## Generate plot
svg(file = "R_exports/m1_trt_eff/inc_trteff.svg")
ggplot(melt.inc.em, aes(x = value, y = variable)) +
  coord_cartesian(ylim = c(1, 8.2)) +
  xlim(-0.5,0.5) +
  scale_y_discrete(limits = unique(rev(melt.inc.em$variable))) +
  geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
             linewidth = 1.75, color = "grey10") +
  geom_density_ridges2(rel_min_height = 5e-11, scale = 0.95,
                       stat = "density_ridges_HPDCrI", 
                       quantile_lines = TRUE,
                       quantiles = c(0.95),
                       quantile_fun = HPD_fun,
                       fill = colsinter[2],
                       alpha = 0.8) +  #<-- transparency (0 more -- 1 less/solid)
  theme_classic(base_size=15) +
  ylab("") + xlab("Treatment Effect Size (95% Credible Interval)")
dev.off()


#svl
set.seed(330)
bSVLmodel1 <- MCMCglmm(logSVL ~ 1 + Species + EggMass + nTreatment  + Species*nTreatment, 
                       random = ~us(1 + nTreatment):Cage, pr = TRUE,
                       data = hatchmorph, nitt = 275000, thin = 50, burnin = 25000,
                       verbose = TRUE, prior = prior1)

# Pull marginals posterior modes for overall model and Species fixed effect
svltr.out <- emtrends(object = bSVLmodel1, pairwise ~ 1, var="nTreatment", data = hatchmorph) #for BOTTOMMOST row; all species + overall
spp.svl.1 <- emtrends(object = bSVLmodel1, pairwise ~ Species, var="nTreatment", data=hatchmorph)

# Bind columns of posterior estimates, melt into proper format for plotting
spp.svl.mc <- as.data.frame(as.mcmc.emmGrid(spp.svl.1$emtrends))
svltr.df.mc <- as.data.frame(as.mcmc.emmGrid(svltr.out$emtrends))
colnames(spp.svl.mc)<-gsub("Species","",colnames(spp.svl.mc))
colnames(svltr.df.mc)<-gsub("1 overall","Overall",colnames(svltr.df.mc))
svl.em <- cbind(spp.svl.mc, svltr.df.mc)
melt.svl.em <- melt(svl.em)

## Generate plot
svg(file = "R_exports/m1_trt_eff/svl_trteff.svg")
ggplot(melt.svl.em, aes(x = value, y = variable)) +
  coord_cartesian(ylim = c(1.4, 8.2)) +
  xlim(-0.5, 0.5) +
  scale_y_discrete(limits = unique(rev(melt.svl.em$variable))) +
  geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
             linewidth = 1.75, color = "grey10") +
  geom_density_ridges2(rel_min_height = 5e-11, scale = 0.95,
                       stat = "density_ridges_HPDCrI", 
                       quantile_lines = TRUE,
                       quantiles = c(0.95),
                       quantile_fun = HPD_fun,
                       fill = colsinter[3],
                       alpha = 0.8) +  #<-- transparency (0 more -- 1 less/solid)
  theme_classic(base_size=15) +
  ylab("") + xlab("Treatment Effect Size (95% Credible Interval)")
dev.off()

#sprint speed
set.seed(330)
bssmodel1 <- MCMCglmm(logSS ~ 1 + Species + SVL + Stops + Fails +
                        SprintTemp + AgeSprint + nTreatment + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = ss, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = TRUE, prior = prior1)

# Pull marginals posterior modes for overall model and Species fixed effect
sstr.out <- emtrends(object = bssmodel1, pairwise ~ 1, var="nTreatment", data = ss) #for BOTTOMMOST row; all species + overall
spp.ss.1 <- emtrends(object = bssmodel1, pairwise ~ Species, var="nTreatment", data=ss)

# Bind columns of posterior estimates, melt into proper format for plotting
spp.ss.mc <- as.data.frame(as.mcmc.emmGrid(spp.ss.1$emtrends))
sstr.df.mc <- as.data.frame(as.mcmc.emmGrid(sstr.out$emtrends))
colnames(spp.ss.mc)<-gsub("Species","",colnames(spp.ss.mc))
colnames(sstr.df.mc)<-gsub("1 overall","Overall",colnames(sstr.df.mc))
ss.em <- cbind(spp.ss.mc, sstr.df.mc)
melt.ss.em <- melt(ss.em)

## Generate plot
svg(file = "R_exports/m1_trt_eff/ss_trteff.svg")
ggplot(melt.ss.em, aes(x = value, y = variable)) +
  coord_cartesian(ylim = c(1.4, 8.2)) +
  xlim(-0.7, 0.7) +
  scale_y_discrete(limits = unique(rev(melt.ss.em$variable))) +
  geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
             linewidth = 1.75, color = "grey10") +
  geom_density_ridges2(rel_min_height = 5e-11, scale = 0.95,
                       stat = "density_ridges_HPDCrI", 
                       quantile_lines = TRUE,
                       quantiles = c(0.95),
                       quantile_fun = HPD_fun,
                       fill = colsinter[4],
                       alpha = 0.8) +  #<-- transparency (0 more -- 1 less/solid)
  theme_classic(base_size=15) +
  ylab("") + xlab("Treatment Effect Size (95% Credible Interval)")
dev.off()

#endurance
set.seed(330)
bendmodel1 <- MCMCglmm(logEnd ~ 1 + Species + nTreatment + SVL + Velocity_cm.s +
                         AgeEndurance + Species*nTreatment, 
                       random = ~us(1 + nTreatment):Cage, pr = TRUE,
                       data = endur, nitt = 275000, thin = 50, burnin = 25000,
                       verbose = TRUE, prior = prior1)

# Pull marginals posterior modes for overall model and Species fixed effect
endurtr.out <- emtrends(object = bendmodel1, pairwise ~ 1, var="nTreatment", data = endur) #for BOTTOMMOST row; all species + overall
spp.endur.1 <- emtrends(object = bendmodel1, pairwise ~ Species, var="nTreatment", data=endur)

# Bind columns of posterior estimates, melt into proper format for plotting
spp.endur.mc <- as.data.frame(as.mcmc.emmGrid(spp.endur.1$emtrends))
endurtr.df.mc <- as.data.frame(as.mcmc.emmGrid(endurtr.out$emtrends))
colnames(spp.endur.mc)<-gsub("Species","",colnames(spp.endur.mc))
colnames(endurtr.df.mc)<-gsub("1 overall","Overall",colnames(endurtr.df.mc))
endur.em <- cbind(spp.endur.mc, endurtr.df.mc)
melt.endur.em <- melt(endur.em)

## Generate plot
svg(file = "R_exports/m1_trt_eff/endur_trteff.svg")
ggplot(melt.endur.em, aes(x = value, y = variable)) +
  coord_cartesian(ylim = c(1.4, 8.2)) +
  xlim(-0.7, 0.7) +
  scale_y_discrete(limits = unique(rev(melt.endur.em$variable))) +
  geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
             linewidth = 1.75, color = "grey10") +
  geom_density_ridges2(rel_min_height = 5e-11, scale = 0.95,
                       stat = "density_ridges_HPDCrI", 
                       quantile_lines = TRUE,
                       quantiles = c(0.95),
                       quantile_fun = HPD_fun,
                       fill = colsinter[5],
                       alpha = 0.8) +  #<-- transparency (0 more -- 1 less/solid)
  theme_classic(base_size=15) +
  ylab("") + xlab("")
dev.off()

##### Paper figures - Intraspecific####

#Have to calculate predicted effects and then convert them back from log to normal scale before plotting
X <- bwumodel1$X
Z <- bwumodel1$Z
unique(X[, "nTreatment"])  #<-- original levels
# Get NewCage names:
Zclnms <- strsplit(Z@Dimnames[[2L]], split = ".", fixed = TRUE)
nms <- unique(unlist(lapply(Zclnms, FUN = "[", i = 3)))
nNewCage <- length(nms)
Spp <- c("carolinensis", "chlorocyanus", "cristatellus", "cybotes", "distichus",
          "equestris", "sagrei") #this has to be in alphabetical order for the apply line to work
nSpp <- length(Spp)
# Get Cage in order of data (i.e., order in X)
## first get subset of Z for just random intercepts (first half)
Zint <- Z[, 1:nNewCage]
# replace column names with simple cage name
Zint@Dimnames[[2L]] <- unlist(lapply(strsplit(Zint@Dimnames[[2L]],
                                              split = ".", fixed = TRUE),
                                     FUN = "[", i = 3))
# re-create data.frame from X design matrix (with only useful variables)
datf <- as.data.frame(as.matrix(X))[, c("EggMass", "nTreatment", "PercDev")]
# use Z design matrix for cage intercepts to assign cages to each observation
datf$Cage <- as.character(Zint@Dimnames[[2L]][(Zint %*% matrix(seq(70), ncol = 1))@x])
# Use X design matrix to assign populations to each observation
datf$Species <- Spp[1 + apply(as.matrix(X[, paste0("Species", Spp)[-1]]),
                                  MARGIN = 1,
                                  FUN = function(r){ ifelse(sum(r) > 0, which(r == 1), 0) })]
# Need placeholder response for MCMCglmm to find and be happy
## Also do predictions below marginal to other fixed effects
sppMeanCovs <- aggregate(cbind(EggMass, PercDev) ~ Species,
                         data = datf, FUN = mean)
### Map these back into datf
meanMap <- match(datf$Species, sppMeanCovs$Species)

datf <- within(datf, {
  logWU <- 1  #<-- placeholder for predict.MCMCglmm
  # within-population average values of covariates  
  ## add tiny bit of random noise or else MCMCglmm complains about priors not
  ### strong enough to estimate fixed effects and drops these
  EggMass <- sppMeanCovs[meanMap, "EggMass"] + rnorm(nrow(datf), 0, 0.005)
  PercDev <- sppMeanCovs[meanMap, "PercDev"] + rnorm(nrow(datf), 0, 0.005)
})
datf$predRfx <- predict(bwumodel1, newdata = datf,
                        marginal = NULL,    # so includes random effects in predictions
                        posterior = "all")
## Now get prediction for mean regression
### (make prediction marginal to all random effects = Default in predict.MCMCglmm)
datf$predMean <- predict(bwumodel1, newdata = datf,
                         posterior = "all")
wudatf <- datf

#svg(file = "R_exports/intraspecific_fig/log_scaled/wulog.svg")
#par(mfrow=c(1,7))
for(p in 1:nSpp){
  #print(paste("wulog_",sppnames[p],".svg"), sep="")
  svg(file = paste("R_exports/intraspecific_fig/log_scaled/wu/wulog_",sppnames[p],".svg", sep=""))
  plot(predMean ~ nTreatment, data = wudatf,
       type = "n", 
       axes = FALSE,  
       xlim = c(-0.55, 0.55), 
       ylim = c(-2.5,0.5), #for log scale, c(-2.5, .5)
       xlab = NA, ylab = NA)
  axis(1, at=c(-0.55,0.55),labels=NA, lwd=10)
  axis(2, at=c(-2.5,0.5), labels=NA, lwd=10)  #<-- so it ends up close to plotting area for "intecerpt line"
  
  # Subset data for population
  sub_datf <- subset(wudatf, Species == sppnames[p])

  ### Go through each Cage to plot cage reaction norms:
  uCageID <- unique(as.character(sub_datf$Cage))
  for(i in uCageID){
    i_int <- match(i, uCageID)  #<-- because group IDs are character/factors

    cagesub_datf <- subset(sub_datf, Cage == i)
    cagesub_datf <- cagesub_datf[order(cagesub_datf$nTreatment), ]

    ### Plot Cage model prediction
    if(nrow(cagesub_datf)>1){
      lines(x = cagesub_datf$nTreatment,
            y= cagesub_datf$predRfx,
            col = cols[p], lwd = 6,
            type="o") #<-- could make color ramp within each population
    } else{
      print(c(p,i,cagesub_datf$nTreatment))
    }
  }
  osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
        y = osub_datf$predMean,
         col = "black", lwd = 10, lty=2, type="l")
  dev.off()
} 


### Incubation Period
X <- bincmodel1$X
Z <- bincmodel1$Z
unique(X[, "nTreatment"])  #<-- original levels
# Get NewCage names:
Zclnms <- strsplit(Z@Dimnames[[2L]], split = ".", fixed = TRUE)
nms <- unique(unlist(lapply(Zclnms, FUN = "[", i = 3)))
nNewCage <- length(nms)
Spp <- c("carolinensis", "chlorocyanus", "cristatellus", "cybotes", "distichus",
         "equestris", "sagrei") #this has to be in alphabetical order for the apply line to work
nSpp <- length(Spp)
# Get Cage in order of data (i.e., order in X)
## first get subset of Z for just random intercepts (first half)
Zint <- Z[, 1:nNewCage]
# replace column names with simple cage name
Zint@Dimnames[[2L]] <- unlist(lapply(strsplit(Zint@Dimnames[[2L]],
                                              split = ".", fixed = TRUE),
                                     FUN = "[", i = 3))
# re-create data.frame from X design matrix (with only useful variables)
datf <- as.data.frame(as.matrix(X))[, c("EggMass", "nTreatment")]
# use Z design matrix for cage intercepts to assign cages to each observation
datf$Cage <- as.character(Zint@Dimnames[[2L]][(Zint %*% matrix(seq(70), ncol = 1))@x])
# Use X design matrix to assign populations to each observation
datf$Species <- Spp[1 + apply(as.matrix(X[, paste0("Species", Spp)[-1]]),
                              MARGIN = 1,
                              FUN = function(r){ ifelse(sum(r) > 0, which(r == 1), 0) })]
# Need placeholder response for MCMCglmm to find and be happy
## Also do predictions below marginal to other fixed effects
sppMeanCovs <- aggregate(cbind(EggMass) ~ Species,
                         data = datf, FUN = mean)
### Map these back into datf
meanMap <- match(datf$Species, sppMeanCovs$Species)

datf <- within(datf, {
  logInc <- 1  #<-- placeholder for predict.MCMCglmm
  # within-population average values of covariates  
  ## add tiny bit of random noise or else MCMCglmm complains about priors not
  ### strong enough to estimate fixed effects and drops these
  EggMass <- sppMeanCovs[meanMap, "EggMass"] + rnorm(nrow(datf), 0, 0.005)
})
datf$predRfx <- predict(bincmodel1, newdata = datf,
                        marginal = NULL,    # so includes random effects in predictions
                        posterior = "all")
## Now get prediction for mean regression
### (make prediction marginal to all random effects = Default in predict.MCMCglmm)
datf$predMean <- predict(bincmodel1, newdata = datf,
                         posterior = "all")
incdatf <- datf

for(p in 1:nSpp){
  
  #png(filename=paste(popnames[p], "_max_intra_std.png", sep=""), width=5, height=4, res=300, units="in")
  svg(file = paste("R_exports/intraspecific_fig/log_scaled/inc/inclog_",sppnames[p],".svg", sep=""))
  plot(predMean ~ nTreatment, data = incdatf,
       type = "n", axes = FALSE,  
       xlim = c(-0.55, 0.55),
       ylim = c(3.2,4.4),
       xlab = NA, ylab = NA)
  axis(1, at = c(-0.55,0.55), lwd=10, labels=NA)
  axis(2, at=c(3.2,4.4), lwd=10, labels=NA)  #<-- so it ends up close to plotting area for "intecerpt line"
  
  # Subset data for population
  sub_datf <- subset(incdatf, Species == sppnames[p])

  
  ### Go through each Cage to plot cage reaction norms:
  uCageID <- unique(as.character(sub_datf$Cage))
  for(i in uCageID){
    # First give an integer index value from `i` which is a character
    i_int <- match(i, uCageID)  #<-- because group IDs are character/factors
    ### Subset data for i-th Cage
    cagesub_datf <- subset(sub_datf, Cage == i)
    cagesub_datf <- cagesub_datf[order(cagesub_datf$nTreatment), ]
    ### Plot Cage model prediction
    if(nrow(cagesub_datf)>1){
      lines(x = cagesub_datf$nTreatment,
            y = cagesub_datf$predRfx,
            col = cols[p], lwd = 6,
            type="o") #<-- could make color ramp within each population
    } else{
      print(c(p,i,cagesub_datf$nTreatment))
    }
  }
  osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
        y = osub_datf$predMean, #<-- what does model predict for each x?
        col = "black", lwd = 10, lty=2)
  dev.off()
} 



### SVL
X <- bSVLmodel1$X
Z <- bSVLmodel1$Z
unique(X[, "nTreatment"])  #<-- original levels
# Get NewCage names:
Zclnms <- strsplit(Z@Dimnames[[2L]], split = ".", fixed = TRUE)
nms <- unique(unlist(lapply(Zclnms, FUN = "[", i = 3)))
nNewCage <- length(nms)
Spp <- c("carolinensis", "chlorocyanus", "cristatellus", "cybotes", "distichus",
         "equestris", "sagrei") #this has to be in alphabetical order for the apply line to work
nSpp <- length(Spp)
# Get Cage in order of data (i.e., order in X)
## first get subset of Z for just random intercepts (first half)
Zint <- Z[, 1:nNewCage]
# replace column names with simple cage name
Zint@Dimnames[[2L]] <- unlist(lapply(strsplit(Zint@Dimnames[[2L]],
                                              split = ".", fixed = TRUE),
                                     FUN = "[", i = 3))
# re-create data.frame from X design matrix (with only useful variables)
datf <- as.data.frame(as.matrix(X))[, c("EggMass", "nTreatment")]
# use Z design matrix for cage intercepts to assign cages to each observation
datf$Cage <- as.character(Zint@Dimnames[[2L]][(Zint %*% matrix(seq(70), ncol = 1))@x])
# Use X design matrix to assign populations to each observation
datf$Species <- Spp[1 + apply(as.matrix(X[, paste0("Species", Spp)[-1]]),
                              MARGIN = 1,
                              FUN = function(r){ ifelse(sum(r) > 0, which(r == 1), 0) })]
# Need placeholder response for MCMCglmm to find and be happy
## Also do predictions below marginal to other fixed effects
sppMeanCovs <- aggregate(cbind(EggMass) ~ Species,
                         data = datf, FUN = mean)
### Map these back into datf
meanMap <- match(datf$Species, sppMeanCovs$Species)

datf <- within(datf, {
  logSVL <- 1  #<-- placeholder for predict.MCMCglmm
  # within-population average values of covariates  
  ## add tiny bit of random noise or else MCMCglmm complains about priors not
  ### strong enough to estimate fixed effects and drops these
  EggMass <- sppMeanCovs[meanMap, "EggMass"] + rnorm(nrow(datf), 0, 0.005)
})
datf$predRfx <- predict(bSVLmodel1, newdata = datf,
                        marginal = NULL,    # so includes random effects in predictions
                        posterior = "all")
## Now get prediction for mean regression
### (make prediction marginal to all random effects = Default in predict.MCMCglmm)
datf$predMean <- predict(bSVLmodel1, newdata = datf,
                         posterior = "all")
svldatf <- datf

for(p in 1:nSpp){
  
  svg(file = paste("R_exports/intraspecific_fig/log_scaled/svl/svllog_",sppnames[p],".svg", sep=""))
  plot(predMean ~ nTreatment, data = svldatf,
       type = "n", axes = FALSE,  
       xlim = c(-0.55, 0.55),
       ylim = c(2,4),
       xlab = NA, ylab = NA)
  axis(1, at = c(-0.55,0.55), lwd=10, labels=NA)
  axis(2, at = c(2,4), lwd=10, labels=NA)  #<-- so it ends up close to plotting area for "intecerpt line"
  
  # Subset data for population
  sub_datf <- subset(svldatf, Species == sppnames[p])

  ### Go through each Cage to plot cage reaction norms:
  uCageID <- unique(as.character(sub_datf$Cage))
  for(i in uCageID){
    # First give an integer index value from `i` which is a character
    i_int <- match(i, uCageID)  #<-- because group IDs are character/factors
    ### Subset data for i-th Cage
    cagesub_datf <- subset(sub_datf, Cage == i)
    #### Check/Order subsetted data by increasing nTreatment
    ##### Needed to plot straight lines
    cagesub_datf <- cagesub_datf[order(cagesub_datf$nTreatment), ]

    ### Plot Cage model prediction
    if(nrow(cagesub_datf)>1){
      lines(cagesub_datf$predRfx ~ cagesub_datf$nTreatment,
            col = cols[p], lwd = 6,
            type="o") #<-- could make color ramp within each population
    } else{
      print(c(p,i,cagesub_datf$nTreatment))
    }
  }
  osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
        y = osub_datf$predMean, #<-- what does model predict for each x?
        col = "black", lwd = 10, lty=2)
  dev.off()
} 


### Sprint Speed
X <- bssmodel1$X
Z <- bssmodel1$Z
unique(X[, "nTreatment"])  #<-- original levels
# Get NewCage names:
Zclnms <- strsplit(Z@Dimnames[[2L]], split = ".", fixed = TRUE)
nms <- unique(unlist(lapply(Zclnms, FUN = "[", i = 3)))
nNewCage <- length(nms)
Spp <- c("carolinensis", "chlorocyanus", "cristatellus", "cybotes", "distichus",
         "equestris", "sagrei") #this has to be in alphabetical order for the apply line to work
nSpp <- length(Spp)
# Get Cage in order of data (i.e., order in X)
## first get subset of Z for just random intercepts (first half)
Zint <- Z[, 1:nNewCage]
# replace column names with simple cage name
Zint@Dimnames[[2L]] <- unlist(lapply(strsplit(Zint@Dimnames[[2L]],
                                              split = ".", fixed = TRUE),
                                     FUN = "[", i = 3))
# re-create data.frame from X design matrix (with only useful variables)
datf <- as.data.frame(as.matrix(X))[, c("SVL", "nTreatment", "Stops", "Fails", "SprintTemp", "AgeSprint")]
# use Z design matrix for cage intercepts to assign cages to each observation
datf$Cage <- as.character(Zint@Dimnames[[2L]][(Zint %*% matrix(seq(68), ncol = 1))@x]) #number in "seq" must equal number of cages
# Use X design matrix to assign populations to each observation
datf$Species <- Spp[1 + apply(as.matrix(X[, paste0("Species", Spp)[-1]]),
                              MARGIN = 1,
                              FUN = function(r){ ifelse(sum(r) > 0, which(r == 1), 0) })]
# Need placeholder response for MCMCglmm to find and be happy
## Also do predictions below marginal to other fixed effects
sppMeanCovs <- aggregate(cbind(SVL, Stops, Fails, SprintTemp, AgeSprint) ~ Species,
                         data = datf, FUN = mean)
### Map these back into datf
meanMap <- match(datf$Species, sppMeanCovs$Species)

datf <- within(datf, {
  logSS <- 1  #<-- placeholder for predict.MCMCglmm
  # within-population average values of covariates  
  ## add tiny bit of random noise or else MCMCglmm complains about priors not
  ### strong enough to estimate fixed effects and drops these
  SVL <- sppMeanCovs[meanMap, "SVL"] + rnorm(nrow(datf), 0, 0.005)
  Stops <- sppMeanCovs[meanMap, "Stops"] + rnorm(nrow(datf), 0, 0.005)
  Fails <- sppMeanCovs[meanMap, "Fails"] + rnorm(nrow(datf), 0, 0.005)
  SprintTemp <- sppMeanCovs[meanMap, "SprintTemp"] + rnorm(nrow(datf), 0, 0.005)
  AgeSprint <- sppMeanCovs[meanMap, "AgeSprint"] + rnorm(nrow(datf), 0, 0.005)
})
datf$predRfx <- predict(bssmodel1, newdata = datf,
                        marginal = NULL,    # so includes random effects in predictions
                        posterior = "all")
## Now get prediction for mean regression
### (make prediction marginal to all random effects = Default in predict.MCMCglmm)
datf$predMean <- predict(bssmodel1, newdata = datf,
                         posterior = "all")
ssdatf <- datf

for(p in 1:nSpp){
  
  svg(file = paste("R_exports/intraspecific_fig/log_scaled/ss/sslog_",sppnames[p],".svg", sep=""))
  plot(predRfx ~ nTreatment, data = ssdatf,
       type = "n", axes = FALSE,  
       xlim = c(-0.55, 0.55),
       ylim = c(-2.8,-2),
       xlab = NA, ylab = NA)
  axis(1, at = c(-0.55,0.55), lwd=10, labels=NA)
  axis(2, at = c(-2.8,-2), lwd=10, labels=NA)  #<-- so it ends up close to plotting area for "intecerpt line"
  
  # Subset data for population
  sub_datf <- subset(ssdatf, Species == sppnames[p])

  ### Go through each Cage to plot cage reaction norms:
  uCageID <- unique(as.character(sub_datf$Cage))
  for(i in uCageID){
    # First give an integer index value from `i` which is a character
    i_int <- match(i, uCageID)  #<-- because group IDs are character/factors
    ### Subset data for i-th Cage
    cagesub_datf <- subset(sub_datf, Cage == i)
    #### Check/Order subsetted data by increasing nTreatment
    ##### Needed to plot straight lines
    cagesub_datf <- cagesub_datf[order(cagesub_datf$nTreatment), ]

    ### Plot Cage model prediction
    if(nrow(cagesub_datf)>1){
      lines(cagesub_datf$predRfx ~ cagesub_datf$nTreatment,
            col = cols[p], lwd = 6,
            type="o") #<-- could make color ramp within each population
    } else{
      print(c(p,i,cagesub_datf$nTreatment))
    }
  }
  osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
        y = osub_datf$predMean, #<-- what does model predict for each x?
        col = "black", lwd = 10, lty=2)
  dev.off()
} 


#endurance
X <- bendmodel1$X
Z <- bendmodel1$Z
unique(X[, "nTreatment"])  #<-- original levels
# Get Cage names:
Zclnms <- strsplit(Z@Dimnames[[2L]], split = ".", fixed = TRUE)
nms <- unique(unlist(lapply(Zclnms, FUN = "[", i = 3)))
nNewCage <- length(nms)
Spp <- c("carolinensis", "chlorocyanus", "cristatellus", "cybotes", "distichus",
         "equestris", "sagrei") #this has to be in alphabetical order for the apply line to work
nSpp <- length(Spp)
# Get Cage in order of data (i.e., order in X)
## first get subset of Z for just random intercepts (first half)
Zint <- Z[, 1:nNewCage]
# replace column names with simple cage name
Zint@Dimnames[[2L]] <- unlist(lapply(strsplit(Zint@Dimnames[[2L]],
                                              split = ".", fixed = TRUE),
                                     FUN = "[", i = 3))
# re-create data.frame from X design matrix (with only useful variables)
datf <- as.data.frame(as.matrix(X))[, c("SVL", "nTreatment", "Velocity_cm.s", "AgeEndurance")]
# use Z design matrix for cage intercepts to assign cages to each observation
datf$Cage <- as.character(Zint@Dimnames[[2L]][(Zint %*% matrix(seq(69), ncol = 1))@x]) #number in "seq" must equal number of cages
# Use X design matrix to assign populations to each observation
datf$Species <- Spp[1 + apply(as.matrix(X[, paste0("Species", Spp)[-1]]),
                              MARGIN = 1,
                              FUN = function(r){ ifelse(sum(r) > 0, which(r == 1), 0) })]
# Need placeholder response for MCMCglmm to find and be happy
## Also do predictions below marginal to other fixed effects
sppMeanCovs <- aggregate(cbind(SVL, Velocity_cm.s, AgeEndurance) ~ Species,
                         data = datf, FUN = mean)
### Map these back into datf
meanMap <- match(datf$Species, sppMeanCovs$Species)

datf <- within(datf, {
  logEnd <- 1  #<-- placeholder for predict.MCMCglmm
  # within-population average values of covariates  
  ## add tiny bit of random noise or else MCMCglmm complains about priors not
  ### strong enough to estimate fixed effects and drops these
  SVL <- sppMeanCovs[meanMap, "SVL"] + rnorm(nrow(datf), 0, 0.005)
  Velocity_cm.s <- sppMeanCovs[meanMap, "Velocity_cm.s"] + rnorm(nrow(datf), 0, 0.005)
  AgeEndurance <- sppMeanCovs[meanMap, "AgeEndurance"] + rnorm(nrow(datf), 0, 0.005)
})
datf$predRfx <- predict(bendmodel1, newdata = datf,
                        marginal = NULL,    # so includes random effects in predictions
                        posterior = "all")
## Now get prediction for mean regression
### (make prediction marginal to all random effects = Default in predict.MCMCglmm)
datf$predMean <- predict(bendmodel1, newdata = datf,
                         posterior = "all")
endurdatf <- datf

#currently exactly the same plot as sprint speed, there's an error somewhere
for(p in 1:nSpp){
  
  svg(file = paste("R_exports/intraspecific_fig/log_scaled/endur/endurlog_",sppnames[p],".svg", sep=""))
  plot(predRfx ~ nTreatment, data = endurdatf,
       type = "n", axes = FALSE,  
       xlim = c(-0.55, 0.55),
       ylim = c(4.5,6),
       xlab = NA, ylab = NA)
  axis(1, at = c(-0.55,0.55), lwd=10, labels=NA)
  axis(2, at = c(4.5,6), lwd=10, labels=NA)  #<-- so it ends up close to plotting area for "intecerpt line"
  
  # Subset data for population
  sub_datf <- subset(endurdatf, Species == sppnames[p])
  
  ### Go through each Cage to plot cage reaction norms:
  uCageID <- unique(as.character(sub_datf$Cage))
  for(i in uCageID){
    # First give an integer index value from `i` which is a character
    i_int <- match(i, uCageID)  #<-- because group IDs are character/factors
    ### Subset data for i-th Cage
    cagesub_datf <- subset(sub_datf, Cage == i)
    #### Check/Order subsetted data by increasing nTreatment
    ##### Needed to plot straight lines
    cagesub_datf <- cagesub_datf[order(cagesub_datf$nTreatment), ]
    
    ### Plot Cage model prediction
    if(nrow(cagesub_datf)>1){
      lines(cagesub_datf$predRfx ~ cagesub_datf$nTreatment,
            col = cols[p], lwd = 6,
            type="o") #<-- could make color ramp within each population
    } else{
      print(c(p,i,cagesub_datf$nTreatment))
    }
  }
  osub_datf <- sub_datf[order(sub_datf$nTreatment), ]

  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
        y = osub_datf$predMean, #<-- what does model predict for each x?
        col = "black", lwd = 10, lty=2, pch=NA)
  dev.off()
} 


##### Paper figures - Plot of Model 1 Cage Posteriors ####

prior1 <- list(R = list(V = 1e-12, nu = -2), #<-- Non-informative improper: *marginal* posterior equal to REML estimate
               G = list(G1 = list(V = diag(2)*0.002, nu = 3, 
                                  alpha.mu = c(0,0), alpha.V = diag(2)*10000)))

clP16 <- wes_palette("GrandBudapest2", n = 16, type = "continuous") 
clPurp <- clP16[c(1:6, 15:16)]


#par(mar = c(5.5, 6, 5.5, 2.3))
par(mar = c(3.5, 3, 1, 3))
png(filename="slope_postvars_model1.png", height=2, width=7, units="in", res=300)
par(mar = c(5.5, 6, 5.5, 2.3), mfrow=c(1,5))
## Water Uptake
postPlot(bwumodel1$VCV[,4], plotHist = TRUE,
         xlim = c(0, 0.085), ylim = c(0,55),
         xlab = "Slope Variance", ylab = "Density",
         denscol = clPurp[4])
# Add prior
abline(h = 0.5, col = "grey50", lwd = 4)
mtext(text = expression(bolditalic("Water Uptake")), side=3, cex=1.3,
      adj = 0
      #line = lettLine, 
      )

## Incubation Period
postPlot(bincmodel1$VCV[,4], plotHist = TRUE,
         xlim = c(0, 0.0004), ylim = c(0,30000),
         xlab = "Slope Variance", ylab = "Density",
         denscol = clPurp[4])
# Add prior
abline(h = 0.5, col = "grey50", lwd = 4)
mtext(text = expression(bolditalic("Incubation Period")), side=3, cex=1.3,
      adj = 0
      #line = lettLine, 
)

## SVL
postPlot(bSVLmodel1$VCV[,4], plotHist = TRUE,
         xlim = c(0, 0.0004), ylim = c(0,30000),
         xlab = "Slope Variance", ylab = "Density",
         denscol = clPurp[4])
# Add prior
abline(h = 0.5, col = "grey50", lwd = 4)
mtext(text = expression(bolditalic("Snout-vent Length")), side=3, cex=1.3,
      adj = 0
      #line = lettLine, 
)

## Sprint Speed
postPlot(bssmodel1$VCV[,4], plotHist = TRUE,
         xlim = c(0, 0.04), ylim = c(0,400),
         xlab = "Slope Variance", ylab = "Density",
         denscol = clPurp[4])
# Add prior
abline(h = 0.5, col = "grey50", lwd = 4)
mtext(text = expression(bolditalic("Sprint Speed")), side=3, cex=1.3,
      adj = 0
      #line = lettLine, 
)

## Endurance
postPlot(bendmodel1$VCV[,4], plotHist = TRUE,
         xlim = c(0, 0.08), ylim = c(0,175),
         xlab = "Slope Variance", ylab = "Density",
         denscol = clPurp[4])
# Add prior
abline(h = 0.5, col = "grey50", lwd = 4)
mtext(text = expression(bolditalic("Endurance")), side=3, cex=1.3,
      adj = 0
      #line = lettLine, 
)
dev.off()

