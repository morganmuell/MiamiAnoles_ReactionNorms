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


library(wesanderson)  #<-- contains color palette to use

rm(list=ls())

###-- Import data
#mializ <- read.csv(file="R_imports/miamidat_full.csv", header=TRUE, sep=",")

load("miami-lizards.rda") #***
load("wu_reaction-norms.rda")
load("inc_reaction-norms.rda")
load("SVL_reaction-norms.rda")
load("ss_reaction-norms.rda")
load("end_reaction-norms.rda")
load("ctmax_reaction-norms.rda")
load("Miami_envtavgs.rda")
load("bio_4s.rda") ###load for heat map fig
load("species-average_reaction-norms.rda") ##don't use for phylosig, based on REML results
load("wufittedrns_log.rda")
load("incfittedrns_log.rda")
load("svlfittedrns_log.rda")
load("ssfittedrns_log.rda")
load("endfittedrns_log.rda")
load("model2df.rda") #***
load("model1_results.rda") #***
load("inputdata_trimmed.rda") #***

sppnames <- c("cristatellus", "sagrei", "carolinensis", "chlorocyanus",
              "distichus", "cybotes", "equestris")
phenotypes <- c("WaterUptake", "Inc", "SVL", "SprintSpeed", "Endurance")

colsinter <- viridis(n=5)
cols <- palette(brewer.pal(n=7, name="Set2")) #intraspecific

#mods <- list(bwumodel1, bincmodel, bSVLmodel1, bssmodel1, bendmodel1)
#save(bwumodel1, bincmodel1, bSVLmodel1, bssmodel1, bendmodel1, file="model1_results.rda")
#save(wu, incp, hatchmorph, ss, endur, file="inputdata_trimmed.rda")


####-- Divide df for datasets -- include response and all potential covariates ######### Plotting Data #####
##Water Uptake
#fixed effects: egg mass, % developed, Treatment, incubator, lay day
#random effects: Cage, Species
wutmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                   "EggMass","scEggMass","PercDev","WaterUptake","logWU")]
wu <- wutmp[which(!is.na(wutmp$logWU) & !is.na(wutmp$PercDev)),] #do instead, alter for all variables
rm(wutmp)

##Incubation period
#fixed effects: egg mass, incubator, Treatment, lay day
#random effects: Cage, species
incptmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                     "EggMass","scEggMass","logInc","Inc")]
incp <- incptmp[which(!is.na(incptmp$logInc) & !is.na(incptmp$EggMass)),] #do instead, alter for all variables

rm(incptmp)

##Hatchling morphology (SVL, TL, mass)
#fixed effects: egg mass, treatment, other morphology, incubator, Treatment, lay day
#random effects: Cage, Species
tmphatchmorph <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                           "EggMass","scEggMass","logSVL", "SVL","TL", "Mass", "Sex")]
hatchmorph <- tmphatchmorph[which(!is.na(tmphatchmorph$logSVL) & !is.na(tmphatchmorph$EggMass)),]
rm(tmphatchmorph)

##Sprint Speed
#fixed effects: SVL, SprintTemp, Treatment, Age, Stops, Fails, Sex
#random effects: Cage, Species
sstmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "LayDay",
                   "SVL", "Sex", "logSS", "Stops", "Fails", "SprintTemp",
                   "AgeSprint", "FastestSprint")]
ss <- sstmp[which(!is.na(sstmp$logSS) & !is.na(sstmp$SVL) & !is.na(sstmp$Stops) &
                    !is.na(sstmp$Fails) & !is.na(sstmp$SprintTemp) & !is.na(sstmp$AgeSprint)),]
rm(sstmp)

##Endurance
#fixed effects: SVL, Fastest sprint?, Treatment, Age, include both velocities, Sex
#random effects: Cage, Species
endurtmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "LayDay",
                      "SVL", "Sex", "FastestSprint", "Velocity_m.s", "Velocity_cm.s",
                      "logEnd", "EnduranceTime", "AgeEndurance", "DistanceRun_cm")]
endur <- endurtmp[which(!is.na(endurtmp$SVL) & !is.na(endurtmp$FastestSprint) & !is.na(endurtmp$Velocity_cm.s) &
                          !is.na(endurtmp$logEnd) & !is.na(endurtmp$EnduranceTime) & !is.na(endurtmp$AgeEndurance)),]
rm(endurtmp)

##### Interspecific fig OLD #####
##-- Water Uptake
#wutrim <- wu[!(wu$Cage %in% c("1","4","18","19","21","25","32","36", 
#                              "40","47","48","51","64","65","75",
#                              "52","70","76")),]

#calculate averages
wunetrns <- matrix(NA, nrow=7, ncol=7)
rownames(wunetrns) <- sppnames
colnames(wunetrns) <- c("coolmean", "warmmean", "diff", "cupper", "clower", "wupper", "wlower")

for(i in sppnames){
    spp <- subset(wu, Species==i)
    tmpcool <- subset(spp, nTreatment==-1)
    wunetrns[i,1] <- mean(tmpcool$WaterUptake)
    tmpwarm <- subset(spp, nTreatment==1)
    wunetrns[i,2] <- mean(tmpwarm$WaterUptake)
    wunetrns[i,3] <- round(mean(tmpwarm$WaterUptake) - mean(tmpcool$WaterUptake), digits=4)
    cerr <- std.error(tmpcool$WaterUptake)
    wunetrns[i,4] <- (mean(tmpcool$WaterUptake)+(1.96*cerr))
    wunetrns[i,5] <- (mean(tmpcool$WaterUptake)-(1.96*cerr))
    werr <- std.error(tmpwarm$WaterUptake)
    wunetrns[i,6] <- (mean(tmpwarm$WaterUptake)+(1.96*werr))
    wunetrns[i,7] <- (mean(tmpwarm$WaterUptake)-(1.96*werr))
}

wunetrns #yes

#line plot with equestris
svg(filename="R_exports/interspecific_fig/water-uptake/wu_line-eq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(0, 1.6))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(0,0.5,1,1.5), labels=c(0,0.5,1,1.5))
legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Mass Difference (g)")
while (n < 8){
  cool <- wunetrns[n,1]
  warm <- wunetrns[n,2]
#  points(1, cool, pch=16, col=n)
#  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#inset without equestris
svg(filename="R_exports/interspecific_fig/water-uptake/wu_line-noeq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(0, 0.4))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(0,0.2, 0.4), labels=c(0,0.2, 0.4))
#legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Mass Difference (g)")
while (n < 8){
  cool <- wunetrns[n,1]
  warm <- wunetrns[n,2]
  #points(1, cool, pch=16, col=n)
  #points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#reduced box plot
svg(filename="R_exports/interspecific_fig/water-uptake/wu_box.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,7.5), ylim=c(0, 2))
axis(1, at=c(0,8), labels=NA)
axis(2, at=c(0,0.5,1.0,1.5,2.0), labels=c(0,0.5,1,1.5,2))
title(main=NA, xlab="Species", ylab="Mass Difference (g)")
mtext(sppnames, side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125,7.125), cex=0.8)
while (n < 8){
  cool <- wunetrns[n,1]
  segments(cspac, wunetrns[n,5], cspac, wunetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, wunetrns[n,5], cspac+0.05, wunetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, wunetrns[n,4], cspac+0.05, wunetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, wunetrns[n,7], wspac, wunetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, wunetrns[n,7], wspac+0.05, wunetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, wunetrns[n,6], wspac+0.05, wunetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, wunetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, wunetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,wunetrns[n,1],wspac,wunetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}

dev.off()

#full box plot
#boxplot
svg(file="R_exports/interspecific_fig/water-uptake/wu_box_gg.svg")
modwu <- wu
colsfix <- c("#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494",
             "#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494")
modwu$nTreatment <- as.factor(modwu$nTreatment)
modwu$Species <- factor(modwu$Species, levels=c("cristatellus","distichus","sagrei","cybotes","carolinensis","chlorocyanus","equestris"))
wuggplot <- ggplot(modwu, aes(x=Species, y=WaterUptake, fill=interaction(Species,nTreatment))) + geom_boxplot(width=0.75, show.legend = FALSE) + theme_classic() + scale_fill_manual(values=colsfix)
wuggplot + ylab("Mass Difference (g)")
dev.off()

##--Incubation Period
#calculate averages
incnetrns <- matrix(NA, nrow=7, ncol=7)
rownames(incnetrns) <- sppnames
colnames(incnetrns) <- c("coolmean", "warmmean", "diff", "cupper", "clower", "wupper", "wlower")

for(i in sppnames){
  spp <- subset(incp, Species==i)
  tmpcool <- subset(spp, nTreatment==-1)
  incnetrns[i,1] <- mean(tmpcool$Inc)
  tmpwarm <- subset(spp, nTreatment==1)
  incnetrns[i,2] <- mean(tmpwarm$Inc)
  incnetrns[i,3] <- round(mean(tmpwarm$Inc) - mean(tmpcool$Inc), digits=4)
  cerr <- std.error(tmpcool$Inc)
  incnetrns[i,4] <- (mean(tmpcool$Inc)+(1.96*cerr))
  incnetrns[i,5] <- (mean(tmpcool$Inc)-(1.96*cerr))
  werr <- std.error(tmpwarm$Inc)
  incnetrns[i,6] <- (mean(tmpwarm$Inc)+(1.96*werr))
  incnetrns[i,7] <- (mean(tmpwarm$Inc)-(1.96*werr))
}

incnetrns #yes

#line plot with equestris
svg(filename="R_exports/interspecific_fig/inc/inc_line-eq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(0, 71))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(0,10,20,30,40,50,60,70), labels=c(0,10,20,30,40,50,60,70))
legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Number of Days to Hatch")
while (n < 8){
  cool <- incnetrns[n,1]
  warm <- incnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#line plot without equestris
svg(filename="R_exports/interspecific_fig/inc/inc_line-noeq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(25, 55))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(25,35,45,55), labels=c(25,35,45,55))
#legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Number of Days to Hatch")
while (n < 7){
  cool <- incnetrns[n,1]
  warm <- incnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#reduced box plot
svg(filename="R_exports/interspecific_fig/inc/inc_box.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,7.5), ylim=c(20, 80))
axis(1, at=c(0,8), labels=NA)
axis(2, at=c(20,40,60,80), labels=c(20,40,60,80))
title(main=NA, xlab="Species", ylab="Number of Days to Hatch")
mtext(sppnames, side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125,7.125), cex=0.8)
while (n < 8){
  cool <- incnetrns[n,1]
  segments(cspac, incnetrns[n,5], cspac, incnetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, incnetrns[n,5], cspac+0.05, incnetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, incnetrns[n,4], cspac+0.05, incnetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, incnetrns[n,7], wspac, incnetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, incnetrns[n,7], wspac+0.05, incnetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, incnetrns[n,6], wspac+0.05, incnetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, incnetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, incnetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,incnetrns[n,1],wspac,incnetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}
dev.off()

#boxplot
svg(filename = "R_exports/interspecific_fig/inc/inc_box_gg.svg")
modincp <- incp
colsfix <- c("#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494",
             "#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494")
modincp$nTreatment <- as.factor(modincp$nTreatment)
modincp$Species <- factor(modincp$Species, levels=c("cristatellus","distichus","sagrei","cybotes","carolinensis","chlorocyanus","equestris"))
incggplot <- ggplot(modincp, aes(x=Species, y=Inc, fill=interaction(Species,nTreatment))) + geom_boxplot(width=0.75, show.legend = FALSE) + theme_classic() + scale_fill_manual(values=colsfix)
incggplot + ylab("Time to hatch (days)") #reorder colors
dev.off()

##--SVL
#calculate averages
svlnetrns <- matrix(NA, nrow=7, ncol=7)
rownames(svlnetrns) <- sppnames
colnames(svlnetrns) <- c("coolmean", "warmmean", "diff", "cupper", "clower", "wupper", "wlower")

for(i in sppnames){
  spp <- subset(hatchmorph, Species==i)
  tmpcool <- subset(spp, nTreatment==-1)
  svlnetrns[i,1] <- mean(tmpcool$SVL)
  tmpwarm <- subset(spp, nTreatment==1)
  svlnetrns[i,2] <- mean(tmpwarm$SVL)
  svlnetrns[i,3] <- round(mean(tmpwarm$SVL) - mean(tmpcool$SVL), digits=4)
  cerr <- std.error(tmpcool$SVL)
  svlnetrns[i,4] <- (mean(tmpcool$SVL)+(1.96*cerr))
  svlnetrns[i,5] <- (mean(tmpcool$SVL)-(1.96*cerr))
  werr <- std.error(tmpwarm$SVL)
  svlnetrns[i,6] <- (mean(tmpwarm$SVL)+(1.96*werr))
  svlnetrns[i,7] <- (mean(tmpwarm$SVL)-(1.96*werr))
}

svlnetrns #yes

#line plot with equestris
svg(filename="R_exports/interspecific_fig/svl/svl_line-eq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(0, 45))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(0,15,30,45), labels=c(0,15,30,45))
legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Snout-vent length (mm)")
while (n < 8){
  cool <- svlnetrns[n,1]
  warm <- svlnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#line plot without equestris
svg(filename="R_exports/interspecific_fig/svl/svl_line-noeq.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(15, 25))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(15,17.5,20,22.5,25), labels=c(15,17.5,20,22.5,25))
title(main=NA, xlab="Incubation Treatment", ylab="Snout-vent length (mm)")
while (n < 7){
  cool <- svlnetrns[n,1]
  warm <- svlnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()



#reduced box plot full
svg(filename="R_exports/interspecific_fig/svl/svl_box.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,7.5), ylim=c(15, 45))
axis(1, at=c(0,7.5), labels=NA)
axis(2, at=c(15,30,45), labels=c(15,30,45))
title(main=NA, xlab="Species", ylab="Snout-vent length (mm)")
mtext(sppnames, side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125,7.125), cex=0.8)
while (n < 8){
  cool <- svlnetrns[n,1]
  segments(cspac, svlnetrns[n,5], cspac, svlnetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, svlnetrns[n,5], cspac+0.05, svlnetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, svlnetrns[n,4], cspac+0.05, svlnetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, svlnetrns[n,7], wspac, svlnetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, svlnetrns[n,7], wspac+0.05, svlnetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, svlnetrns[n,6], wspac+0.05, svlnetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, svlnetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, svlnetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,svlnetrns[n,1],wspac,svlnetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}
dev.off()

#reduced box plot no eq + eq separate
svg(filename="R_exports/interspecific_fig/svl/svl_box_noeq1.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,6.5), ylim=c(17, 25))
axis(1, at=c(0,7.5), labels=NA)
axis(2, at=c(17.5,20,22.5,25), labels=c(17.5,20,22.5,25))
title(main=NA, xlab="Species", ylab="Snout-vent length (mm)")
mtext(sppnames[c(1:6)], side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125), cex=0.8)
while (n < 7){
  cool <- svlnetrns[n,1]
  segments(cspac, svlnetrns[n,5], cspac, svlnetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, svlnetrns[n,5], cspac+0.05, svlnetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, svlnetrns[n,4], cspac+0.05, svlnetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, svlnetrns[n,7], wspac, svlnetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, svlnetrns[n,7], wspac+0.05, svlnetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, svlnetrns[n,6], wspac+0.05, svlnetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, svlnetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, svlnetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,svlnetrns[n,1],wspac,svlnetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}
dev.off()

#boxplot
svg(filename = "R_exports/interspecific_fig/svl/svl_box_gg.svg")
modhatchmorph <- hatchmorph
colsfix <- c("#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494",
             "#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494")
modhatchmorph$nTreatment <- as.factor(modhatchmorph$nTreatment)
modhatchmorph$Species <- factor(modhatchmorph$Species, levels=c("cristatellus","distichus","sagrei","cybotes","carolinensis","chlorocyanus","equestris"))
svlgg <- ggplot(modhatchmorph, aes(x=Species, y=SVL, fill=interaction(Species,nTreatment))) + geom_boxplot(width=0.75, show.legend = FALSE) + theme_classic() + scale_fill_manual(values=colsfix)
svlgg + ylab("Snout-vent length (mm)") #reorder colors
dev.off()

##--Sprint Speed
#calculate averages
ssnetrns <- matrix(NA, nrow=7, ncol=7)
rownames(ssnetrns) <- sppnames
colnames(ssnetrns) <- c("coolmean", "warmmean", "diff", "cupper", "clower", "wupper", "wlower")

for(i in sppnames){
  spp <- subset(ss, Species==i)
  tmpcool <- subset(spp, nTreatment==-1)
  ssnetrns[i,1] <- mean(tmpcool$FastestSprint)
  tmpwarm <- subset(spp, nTreatment==1)
  ssnetrns[i,2] <- mean(tmpwarm$FastestSprint)
  ssnetrns[i,3] <- round(mean(tmpwarm$FastestSprint) - mean(tmpcool$FastestSprint), digits=4)
  cerr <- std.error(tmpcool$FastestSprint)
  ssnetrns[i,4] <- (mean(tmpcool$FastestSprint)+(1.96*cerr))
  ssnetrns[i,5] <- (mean(tmpcool$FastestSprint)-(1.96*cerr))
  werr <- std.error(tmpwarm$FastestSprint)
  ssnetrns[i,6] <- (mean(tmpwarm$FastestSprint)+(1.96*werr))
  ssnetrns[i,7] <- (mean(tmpwarm$FastestSprint)-(1.96*werr))
}

ssnetrns #yes

#line plot with equestris
svg(filename="R_exports/interspecific_fig/sprint-speed/ss_line.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(0, 2))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(0,0.5,1,1.5,2), labels=c(0,0.5,1,1.5,2))
legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
title(main=NA, xlab="Incubation Treatment", ylab="Fastest Sprint (m/s)")
while (n < 8){
  cool <- ssnetrns[n,1]
  warm <- ssnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#reduced box plot
svg(filename="R_exports/interspecific_fig/sprint-speed/ss_box_2.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,7.5), ylim=c(0, 2))
axis(1, at=c(0,7.5), labels=NA)
axis(2, at=c(0,0.5,1,1.5,2), labels=c(0,0.5,1,1.5,2))
title(main=NA, xlab="Species", ylab="Fastest Sprint (m/s")
mtext(sppnames, side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125,7.125), cex=0.8)
while (n < 8){
  cool <- ssnetrns[n,1]
  segments(cspac, ssnetrns[n,5], cspac, ssnetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, ssnetrns[n,5], cspac+0.05, ssnetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, ssnetrns[n,4], cspac+0.05, ssnetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, ssnetrns[n,7], wspac, ssnetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, ssnetrns[n,7], wspac+0.05, ssnetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, ssnetrns[n,6], wspac+0.05, ssnetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, ssnetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, ssnetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,ssnetrns[n,1],wspac,ssnetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}
dev.off()

#box plot
svg(filename = "R_exports/interspecific_fig/sprint-speed/ss_box_gg.svg")
modss <- ss
colsfix <- c("#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494",
             "#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494")
modss$nTreatment <- as.factor(modss$nTreatment)
modss$Species <- factor(modss$Species, levels=c("cristatellus","distichus","sagrei","cybotes","carolinensis","chlorocyanus","equestris"))
ssgg <- ggplot(modss, aes(x=Species, y=FastestSprint, fill=interaction(Species,nTreatment))) + geom_boxplot(width=0.75, show.legend = FALSE) + theme_classic() + scale_fill_manual(values=colsfix)
ssgg + ylab("Fastest Sprint (m/s)")  #reorder colors
dev.off()


##--Endurance
#calculate averages
endnetrns <- matrix(NA, nrow=7, ncol=7)
rownames(endnetrns) <- sppnames
colnames(endnetrns) <- c("coolmean", "warmmean", "diff", "cupper", "clower", "wupper", "wlower")

for(i in sppnames){
  spp <- subset(endur, Species==i)
  tmpcool <- subset(spp, nTreatment==-1)
  endnetrns[i,1] <- mean(tmpcool$DistanceRun_cm)
  tmpwarm <- subset(spp, nTreatment==1)
  endnetrns[i,2] <- mean(tmpwarm$DistanceRun_cm)
  endnetrns[i,3] <- round(mean(tmpwarm$DistanceRun_cm) - mean(tmpcool$DistanceRun_cm), digits=4)
  cerr <- std.error(tmpcool$DistanceRun_cm)
  endnetrns[i,4] <- (mean(tmpcool$DistanceRun_cm)+(1.96*cerr))
  endnetrns[i,5] <- (mean(tmpcool$DistanceRun_cm)-(1.96*cerr))
  werr <- std.error(tmpwarm$DistanceRun_cm)
  endnetrns[i,6] <- (mean(tmpwarm$DistanceRun_cm)+(1.96*werr))
  endnetrns[i,7] <- (mean(tmpwarm$DistanceRun_cm)-(1.96*werr))
}

endnetrns #yes

#line plot with equestris
svg(filename="R_exports/interspecific_fig/endurance/endur_line.svg")
n <- 1
plot.new()
plot.window(xlim=c(0,3), ylim=c(140, 280))
axis(1, at=c(1,2), labels=c("cool", "warm"))
axis(1, at=c(0,3), labels=NA)
axis(2, at=c(140,170,200,230,260,290), labels=c(140,170,200,230,260,290))
title(main=NA, xlab="Incubation Treatment", ylab="Distance Run (cm)")
legend("topright", legend=sppnames, fill=cols, cex=0.7, bty='n')
while (n < 8){
  cool <- endnetrns[n,1]
  warm <- endnetrns[n,2]
  #  points(1, cool, pch=16, col=n)
  #  points(2, warm, pch=16, col=n)
  lines(c(cool,warm) ~ c(1,2), col=n, lwd=4)
  n <- n+1
}
dev.off()

#reduced box plot
svg(filename="R_exports/interspecific_fig/endurance/endur_box_2.svg")
n <- 1
cspac <- 1
wspac <- 1.25
plot.new()
plot.window(xlim=c(0.5,7.5), ylim=c(110,320))
axis(1, at=c(0,7.5), labels=NA)
axis(2, at=c(110,140,170,200,230,260,290,320), labels=c(110,140,170,200,230,260,290,320))
title(main=NA, xlab="Species", ylab="Distance run (cm)")
mtext(sppnames, side=1, at=c(1.125,2.125,3.125,4.125,5.125,6.125,7.125), cex=0.8)
while (n < 8){
  cool <- endnetrns[n,1]
  segments(cspac, endnetrns[n,5], cspac, endnetrns[n,4], col=n, lwd=4) #vertical line
  segments(cspac-0.05, endnetrns[n,5], cspac+0.05, endnetrns[n,5], col=n, lwd=4) #lower horizontal line
  segments(cspac-0.05, endnetrns[n,4], cspac+0.05, endnetrns[n,4], col=n, lwd=4) #upper horizontal line
  warm <- wunetrns[n,2]
  segments(wspac, endnetrns[n,7], wspac, endnetrns[n,6], col=n, lwd=4)
  segments(wspac-0.05, endnetrns[n,7], wspac+0.05, endnetrns[n,7], col=n, lwd=4) #lower horizontal line
  segments(wspac-0.05, endnetrns[n,6], wspac+0.05, endnetrns[n,6], col=n, lwd=4) #upper horizontal line
  points(cspac, endnetrns[n,1], col=n, pch=18, cex=1.5) #cool
  points(wspac, endnetrns[n,2], col=n, pch=18, cex=1.5) #warm
  segments(cspac,endnetrns[n,1],wspac,endnetrns[n,2], col=n, lwd=4) #draw reaction norm
  n <- n+1
  cspac <- cspac+1
  wspac <- wspac+1
}
dev.off()

svg(filename = "R_exports/interspecific_fig/endurance/end_box_gg.svg")
modendur <- endur
colsfix <- c("#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494",
             "#66C2A5","#A6D854","#FC8D62","#FFD92F","#8DA0CB","#E78AC3","#E5C494")
modendur$nTreatment <- as.factor(modendur$nTreatment)
modendur$Species <- factor(modendur$Species, levels=c("cristatellus","distichus","sagrei","cybotes","carolinensis","chlorocyanus","equestris"))
endurgg <- ggplot(modendur, aes(x=Species, y=DistanceRun_cm, fill=interaction(Species,nTreatment))) + geom_boxplot(width=0.75, show.legend = FALSE) + theme_classic() + scale_fill_manual(values=colsfix)
endurgg + ylab("Distance Run (cm)") + theme(axis.text=element_text(size=12))
dev.off()


## Intraspecific
#carolinensis
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_car.svg")
dumb <- ggplot(subset(wu, Species=="carolinensis"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3, colour=cols[1])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
dumb + theme(axis.text.x=element_blank(), #remove x axis labels
             axis.ticks.x=element_blank(), #remove x axis ticks
             axis.text.y=element_blank(),  #remove y axis labels
             axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_car.svg")
inccar <- ggplot(subset(incp, Species=="carolinensis"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[1])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_car.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="carolinensis"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[1])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_car.svg")
sscar <- ggplot(subset(ss, Species=="carolinensis"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[1])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_car.svg")
endurcar <- ggplot(subset(endur, Species=="carolinensis"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[1])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#chlorocyanus
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_chloro.svg")
chlorodumb <- ggplot(subset(wu, Species=="chlorocyanus"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[2])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
chlorodumb + theme(axis.text.x=element_blank(), #remove x axis labels
                   axis.ticks.x=element_blank(), #remove x axis ticks
                   axis.text.y=element_blank(),  #remove y axis labels
                   axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_chloro.svg")
inccar <- ggplot(subset(incp, Species=="chlorocyanus"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[2])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_chloro.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="chlorocyanus"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[2])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_chloro.svg")
sscar <- ggplot(subset(ss, Species=="chlorocyanus"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[2])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_chloro.svg")
endurcar <- ggplot(subset(endur, Species=="chlorocyanus"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[2])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#cristatellus
svg(filename = "R_exports/intraspecific_fig/water-uptake/wu_crista.svg")
cristadumb <- ggplot(subset(wu, Species=="cristatellus"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[3])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
cristadumb + theme(axis.text.x=element_blank(), #remove x axis labels
                   axis.ticks.x=element_blank(), #remove x axis ticks
                   axis.text.y=element_blank(),  #remove y axis labels
                   axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_crista.svg")
inccar <- ggplot(subset(incp, Species=="cristatellus"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[3])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_crista.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="cristatellus"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[3])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_crista.svg")
sscar <- ggplot(subset(ss, Species=="cristatellus"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[3])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_crista.svg")
endurcar <- ggplot(subset(endur, Species=="cristatellus"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[3])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#cybotes
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_cyb.svg")
cybwu<- ggplot(subset(wu, Species=="cybotes"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[4])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
cybwu + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_cyb.svg")
inccar <- ggplot(subset(incp, Species=="cybotes"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[4])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_cyb.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="cybotes"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[4])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_cyb.svg")
sscar <- ggplot(subset(ss, Species=="cybotes"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[4])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_cyb.svg")
endurcar <- ggplot(subset(endur, Species=="cybotes"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[4])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#distichus
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_dist.svg")
distwu<- ggplot(subset(wu, Species=="distichus"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[5])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
distwu + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_dist.svg")
inccar <- ggplot(subset(incp, Species=="distichus"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[5])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_dist.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="distichus"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[5])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_dist.svg")
sscar <- ggplot(subset(ss, Species=="distichus"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[5])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_dist.svg")
endurcar <- ggplot(subset(endur, Species=="distichus"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[5])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#equestris
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_eq.svg")
wueq <- ggplot(subset(wu, Species=="equestris"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[6])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
wueq + theme(axis.text.x=element_blank(), #remove x axis labels
             axis.ticks.x=element_blank(), #remove x axis ticks
             axis.text.y=element_blank(),  #remove y axis labels
             axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_eq.svg")
inccar <- ggplot(subset(incp, Species=="equestris"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[6])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_eq.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="equestris"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[6])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_eq.svg")
sscar <- ggplot(subset(ss, Species=="equestris"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[6])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_eq.svg")
endurcar <- ggplot(subset(endur, Species=="equestris"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[6])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()
#sagrei
svg(filename="R_exports/intraspecific_fig/water-uptake/wu_sag.svg")
wusag <- ggplot(subset(wu, Species=="sagrei"),aes(x=nTreatment,y=WaterUptake,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[7])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic()
wusag + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/inc/inc_sag.svg")
inccar <- ggplot(subset(incp, Species=="sagrei"),aes(x=nTreatment,y=Inc,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[7])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
inccar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/svl/SVL_sag.svg")
SVLcar <- ggplot(subset(hatchmorph, Species=="sagrei"),aes(x=nTreatment,y=SVL,group=Cage))+ geom_line(alpha=1,lwd=3, colour=cols[7])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
SVLcar + theme(axis.text.x=element_blank(), #remove x axis labels
               axis.ticks.x=element_blank(), #remove x axis ticks
               axis.text.y=element_blank(),  #remove y axis labels
               axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/ss/ss_sag.svg")
sscar <- ggplot(subset(ss, Species=="sagrei"),aes(x=nTreatment,y=FastestSprint,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[7])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
sscar + theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank())
dev.off()
svg(filename="R_exports/intraspecific_fig/endur/endur_sag.svg")
endurcar <- ggplot(subset(endur, Species=="sagrei"),aes(x=nTreatment,y=DistanceRun_cm,group=Cage))+ geom_line(alpha=1, lwd=3,colour=cols[7])+ scale_x_continuous(breaks=c(-1,0,1))+ theme_classic() 
endurcar + theme(axis.text.x=element_blank(), #remove x axis labels
                 axis.ticks.x=element_blank(), #remove x axis ticks
                 axis.text.y=element_blank(),  #remove y axis labels
                 axis.ticks.y=element_blank())
dev.off()

wuplot <- ggplot(wurns, aes(x=nTreatment, y=WaterUptake, color=Cage)) + geom_line(show.legend=FALSE)
wuplot + facet_wrap( ~ wuspp, ncol=3)

noneqwu <- subset(wurns, Species!="equestris")
wuplot2 <- ggplot(noneqwu, aes(x=nTreatment, y=WaterUptake, color=Cage)) + geom_line(show.legend=FALSE)
wuplot2 + facet_wrap( ~ wuspp, ncol=3)


##--Incubation Period
incrns

#cage 2 - int: 39.38498, Slope: -6.873540
-6.87354*-1 + 39.38498 #46.25852
-6.87354 + 39.38498 #32.51144
test <- data.frame(cbind(c(-1,1),c(46.25852,32.51144)))
colnames(test) <- c("Treatment","Response")

incplotvals <- matrix(NA, nrow=nrow(incrns), ncol=3)
plot(Response ~ Treatment, data=test, xlim=c(-2,2), ylim=c(28,50))

for(i in incrns$Cage){
  int <- incrns[i,2]
  slope <- incrns[i,3]
  val1 <- slope*-1 + int
  val2 <- slope + int
  
  
}


incptrim <- incp[!(incp$Cage %in% c("1","21", "40","47","48",
                                    "52","76")),] #this seems to work

noneqinc <- subset(incptrim, Species!="equestris")

ggplot(cristainc, aes(x=nTreatment, y=Inc, color=Cage)) + 
  geom_line(show.legend = FALSE) + geom_point(show.legend=FALSE)

incplot <- ggplot(incptrim, aes(x=nTreatment, y=Inc, color=Cage)) + geom_line(show.legend=FALSE)
incplot + facet_wrap( ~ Species, ncol=3)

incplot2 <- ggplot(noneqinc, aes(x=nTreatment, y=Inc, color=Cage)) + geom_line(show.legend=FALSE)
incplot2 + facet_wrap( ~ Species, ncol=3)

##--SVL
SVLrns
hmtrim <- hatchmorph[!(hatchmorph$Cage %in% c("1","18","21", "40","47","48",
                                              "52","76")),]

noneqhm <- subset(hmtrim, Species!="equestris")

svlplot <- ggplot(hmtrim, aes(x=nTreatment, y=SVL, color=Cage)) + geom_line(show.legend=FALSE)
svlplot + facet_wrap( ~ Species, ncol=3)

svlplot2 <- ggplot(noneqhm, aes(x=nTreatment, y=SVL, color=Cage)) + geom_line(show.legend=FALSE)
svlplot2 + facet_wrap( ~ Species, ncol=3)



##--Sprint Speed
ssrns
sstrim <- ss[!(ss$Cage %in% c("1","18","21","25","32",
                              "36","40","47","48","51",
                              "64","65","72","75","76")),]
ssplot <- ggplot(sstrim, aes(x=nTreatment, y=FastestSprint, color=Cage)) + geom_line(show.legend=FALSE)
ssplot + facet_wrap( ~ Species, ncol=3)

##--Endurance
endrns
endtrim <- endur[!(endur$Cage %in% c("1","18","21","25","32",
                                     "36","40","46","47","48","49","51","61",
                                     "65","72","74","75","76")),]
endplot <- ggplot(endtrim, aes(x=nTreatment, y=DistanceRun_cm, color=Cage)) + geom_line(show.legend=FALSE)
endplot + facet_wrap( ~ Species, ncol=3)

##--CTmax
ctrns
cttrim <- ctmaxdat[!(ctmaxdat$Cage %in% c("1","4","6","7","11","13","14","19","20","21","25","27","31","32",
                                          "33","35","36","42","45","47","48","49","51","57","58","59",
                                          "65","66","67","72","73","74","75","76")),]
ctplot <- ggplot(cttrim, aes(x=nTreatment, y=Ctmax, color=Cage)) + geom_line(show.legend=FALSE)
ctplot + facet_wrap( ~ Species, ncol=3)


##### Paper figures - Interspecific

##### Interspecific fig NEW 3/19/24

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
#bro I cannot do this anymore
#3/19/24 okay this one - the issue before was that data were scaled and centered so line directions didn't make sense
#LOG SCALE should make sense direction-wise, keeps units consistent for comparison among species, so if the only goal is to look at intraspecific variation
#   and compare that directly among species, I should keep it on the log scale and remove units on axes
#   so I think this figure code is perfectly fine and ready to roll, just needs to regenerate with new model objects and new log scaled response values

#I do not have the posterior object for model 1 but I could recreate it
#I really should change the script to spit out the model 1 object as an .Rda every time
#using run 1 random seed which is 330

# I want to quit so bad

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
### ^^^I barely know wtf is going on up there but I know it works - produces predicted response data

#back transform to raw data, not log transformed
##  by the way, doing this exponentiation caused all the plotting problems for some reason
##  when I returned to using the predicted data for plots the problem totally vanished
##  one error implied it might have to do with some dumbass structure, I bet an issue with embedded vectors
#wudatf$predwu <- exp(wudatf$predMean)
#wudatf$rawrfx <- exp(wudatf$predRfx)

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
  #print(sub_datf)
  #### Add Population mean regression line
  # Need to create new data.frame that is ordered according to the covariate
  ## Do this to get a straight line, otherwise lines will just connect-the-dots
  ### in whatever order they appear and gives an awful line
  # osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  #  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
  #        y = osub_datf$predEAHT, #<-- what does model predict for each x?
  #        col = "black", lwd = 5)
  # Add the points now
  #points(TODO)
  
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

### Water Uptake
bwumodel1 #pull from bayesian analyses - make sure it's the right one!

#pull random slopes from model 1
tmpMod <- bwumodel1
SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
CageNms <- SolColNms[87:156]
CageNumsUn <- gsub("nTreatment.","",CageNums) #OHMYGOD FUCKING GSUB COMMAND WTF SO NICE WTF

#construct data frame
wufittedrns <- matrix(NA, nrow=70, ncol=2)
colnames(wufittedrns) <- c("AvgIntercept", "AvgSlope")
rownames(wufittedrns) <- CageNumsUn #fuck yeah brother

#intercepts
wubints <- bwumodel1$Sol[,17:86]
i <- 1
while (i<71){
  vec <- wubints[,i]
  wufittedrns[i,1] <- mean(vec)
  i <- i+1
}

#slopes
wubslopes <- bwumodel1$Sol[,87:156]
i <- 1
while (i<71){
  vec <- wubslopes[,i]
  wufittedrns[i,2] <- mean(vec)
  i <- i+1
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
### ^^^I barely know wtf is going on up there but I know it works - produces predicted response data

#back transform to raw data, not log transformed
#incdatf$predinc <- exp(incdatf$predMean)
#incdatf$rawrfx <- exp(incdatf$predRfx)

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
  #print(sub_datf)
  #### Add Population mean regression line
  # Need to create new data.frame that is ordered according to the covariate
  ## Do this to get a straight line, otherwise lines will just connect-the-dots
  ### in whatever order they appear and gives an awful line
  # osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  #  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
  #        y = osub_datf$predEAHT, #<-- what does model predict for each x?
  #        col = "black", lwd = 5)
  # Add the points now
  #points(TODO)
  
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


bincmodel1

#pull random slopes from model 1
tmpMod <- bincmodel1
SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
CageNums <- SolColNms[86:155]
CageNumsUn <- gsub("nTreatment.","",CageNums) #OHMYGOD FUCKING GSUB COMMAND WTF SO NICE WTF

#construct data frame
incfittedrns <- matrix(NA, nrow=70, ncol=2)
colnames(incfittedrns) <- c("AvgIntercept", "AvgSlope")
rownames(incfittedrns) <- CageNumsUn

#intercepts
incbints <- bincmodel1$Sol[,16:85]
i <- 1
while (i<71){
  vec <- incbints[,i]
  incfittedrns[i,1] <- mean(vec)
  i <- i+1
}

#slopes
incbslopes <- bincmodel1$Sol[,86:155]
i <- 1
while (i<71){
  vec <- incbslopes[,i]
  incfittedrns[i,2] <- mean(vec)
  i <- i+1
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
### ^^^I barely know wtf is going on up there but I know it works - produces predicted response data

#back transform to raw data, not log transformed
#svldatf$predsvl <- exp(svldatf$predMean)
#svldatf$rawrfx <- exp(svldatf$predRfx)

for(p in 1:nSpp){
  
  #png(filename=paste(popnames[p], "_max_intra_std.png", sep=""), width=5, height=4, res=300, units="in")
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

#scale is too different to visualize differences if all species are on same scale with raw values

bSVLmodel1

#pull random slopes from model 1
tmpMod <- bSVLmodel1
SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
CageNums <- SolColNms[86:155]
CageNumsUn <- gsub("nTreatment.","",CageNums)

#construct data frame
svlfittedrns <- matrix(NA, nrow=70, ncol=2)
colnames(svlfittedrns) <- c("AvgIntercept", "AvgSlope")
rownames(svlfittedrns) <- CageNumsUn

#intercepts
svlbints <- bSVLmodel1$Sol[,16:85]
i <- 1
while (i<71){
  vec <- svlbints[,i]
  svlfittedrns[i,1] <- mean(vec)
  i <- i+1
}

#slopes
svlbslopes <- bSVLmodel1$Sol[,86:155]
i <- 1
while (i<71){
  vec <- svlbslopes[,i]
  svlfittedrns[i,2] <- mean(vec)
  i <- i+1
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
### ^^^I barely know wtf is going on up there but I know it works - produces predicted response data

#back transform to raw data, not log transformed
#ssdatf$predss <- exp(ssdatf$predMean)
#ssdatf$rawrfx <- exp(ssdatf$predRfx)

for(p in 1:nSpp){
  
  #png(filename=paste(popnames[p], "_max_intra_std.png", sep=""), width=5, height=4, res=300, units="in")
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
  #print(sub_datf)
  #### Add Population mean regression line
  # Need to create new data.frame that is ordered according to the covariate
  ## Do this to get a straight line, otherwise lines will just connect-the-dots
  ### in whatever order they appear and gives an awful line
  # osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  #  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
  #        y = osub_datf$predEAHT, #<-- what does model predict for each x?
  #        col = "black", lwd = 5)
  # Add the points now
  #points(TODO)
  
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

#well. do I use predicted slopes. Or raw slopes. Or modeled slope effects. BEcause they're different.

bssmodel1

#endurance
bendmodel1

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
### ^^^I barely know wtf is going on up there but I know it works - produces predicted response data

#back transform to raw data, not log transformed
#ssdatf$predss <- exp(ssdatf$predMean)
#ssdatf$rawrfx <- exp(ssdatf$predRfx)

#currently exactly the same plot as sprint speed, there's an error somewhere
for(p in 1:nSpp){
  
  #png(filename=paste(popnames[p], "_max_intra_std.png", sep=""), width=5, height=4, res=300, units="in")
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
  #print(sub_datf)
  #### Add Population mean regression line
  # Need to create new data.frame that is ordered according to the covariate
  ## Do this to get a straight line, otherwise lines will just connect-the-dots
  ### in whatever order they appear and gives an awful line
  # osub_datf <- sub_datf[order(sub_datf$nTreatment), ]
  #  lines(x = osub_datf$nTreatment,   #<-- the entire range of x-axis
  #        y = osub_datf$predEAHT, #<-- what does model predict for each x?
  #        col = "black", lwd = 5)
  # Add the points now
  #points(TODO)
  
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

### OLD ###
#pull random slopes from model 1
tmpMod <- bssmodel1
SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
CageNums <- SolColNms[89:157]
CageNumsUn <- gsub("nTreatment.","",CageNums)

#construct data frame
ssfittedrns <- matrix(NA, nrow=69, ncol=2)
colnames(ssfittedrns) <- c("AvgIntercept", "AvgSlope")
rownames(ssfittedrns) <- CageNumsUn

#intercepts
ssbints <- bssmodel1$Sol[,20:88]
i <- 1
while (i<70){
  vec <- ssbints[,i]
  ssfittedrns[i,1] <- mean(vec)
  i <- i+1
}

#slopes
ssbslopes <- bssmodel1$Sol[,89:157]
i <- 1
while (i<70){
  vec <- ssbslopes[,i]
  ssfittedrns[i,2] <- mean(vec)
  i <- i+1
}

### Endurance
bendmodel1

#pull random slopes from model 1
tmpMod <- bendmodel1
SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
CageNums <- SolColNms[87:155]
CageNumsUn <- gsub("nTreatment.","",CageNums)

#construct data frame
endfittedrns <- matrix(NA, nrow=69, ncol=2)
colnames(endfittedrns) <- c("AvgIntercept", "AvgSlope")
rownames(endfittedrns) <- CageNumsUn

#intercepts
endbints <- bendmodel1$Sol[,18:86]
i <- 1
while (i<70){
  vec <- endbints[,i]
  endfittedrns[i,1] <- mean(vec)
  i <- i+1
}

#slopes
endbslopes <- bendmodel1$Sol[,87:155]
i <- 1
while (i<70){
  vec <- endbslopes[,i]
  endfittedrns[i,2] <- mean(vec)
  i <- i+1
}

#save output
save(wufittedrns, file="wufittedrns_log.rda")
save(incfittedrns, file="incfittedrns_log.rda")
save(svlfittedrns, file="svlfittedrns_log.rda")
save(ssfittedrns, file="ssfittedrns_log.rda")
save(endfittedrns, file="endfittedrns_log.rda")

## Make intraspecific panels
# cristatellus #
wucristafitted <- wufittedrns[c("Cage.1", "Cage.2", "Cage.3", "Cage.4", "Cage.6", "Cage.7", "Cage.8", "Cage.9", "Cage.10"),]
inccristafitted <- incfittedrns[c("Cage.1", "Cage.2", "Cage.3", "Cage.4", "Cage.6", "Cage.7", "Cage.8", "Cage.9", "Cage.10"),]
svlcristafitted <- svlfittedrns[c("Cage.1", "Cage.2", "Cage.3", "Cage.4", "Cage.6", "Cage.7", "Cage.8", "Cage.9", "Cage.10"),]
sscristafitted <- ssfittedrns[c("Cage.1", "Cage.2", "Cage.3", "Cage.4", "Cage.6", "Cage.7", "Cage.8", "Cage.9", "Cage.10"),]
endcristafitted <- endfittedrns[c("Cage.1", "Cage.2", "Cage.3", "Cage.4", "Cage.6", "Cage.7", "Cage.8", "Cage.9", "Cage.10"),]

#water uptake
svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_crista.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5))
axis(1, at=c(1,2), labels=NA, lwd=4)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=4)
n <- 1
while (n < 10){
  tmpint <- wucristafitted[n,1]
  tmpslp <- wucristafitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[1])
  points(2, wtmp, pch=16, col=cols[1])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[1], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_crista.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- inccristafitted[n,1]
  tmpslp <- inccristafitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[1])
  points(2, wtmp, pch=16, col=cols[1])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[1], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_crista.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- svlcristafitted[n,1]
  tmpslp <- svlcristafitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[1])
  points(2, wtmp, pch=16, col=cols[1])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[1], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_crista.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2,0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- sscristafitted[n,1]
  tmpslp <- sscristafitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[1])
  points(2, wtmp, pch=16, col=cols[1])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[1], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_crista.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- endcristafitted[n,1]
  tmpslp <- endcristafitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[1])
  points(2, wtmp, pch=16, col=cols[1])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[1], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# sagrei #
wusagfitted <- wufittedrns[c("Cage.11", "Cage.12", "Cage.13", "Cage.14", "Cage.15", "Cage.16", "Cage.17", "Cage.18","Cage.19", "Cage.20"),]
incsagfitted <- incfittedrns[c("Cage.11", "Cage.12", "Cage.13", "Cage.14", "Cage.15", "Cage.16", "Cage.17", "Cage.18","Cage.19", "Cage.20"),]
svlsagfitted <- svlfittedrns[c("Cage.11", "Cage.12", "Cage.13", "Cage.14", "Cage.15", "Cage.16", "Cage.17", "Cage.18","Cage.19", "Cage.20"),]
sssagfitted <- ssfittedrns[c("Cage.11", "Cage.12", "Cage.13", "Cage.14", "Cage.15", "Cage.16", "Cage.17", "Cage.18","Cage.19", "Cage.20"),]
endsagfitted <- endfittedrns[c("Cage.11", "Cage.12", "Cage.13", "Cage.14", "Cage.15", "Cage.16", "Cage.17", "Cage.18","Cage.19", "Cage.20"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_sag.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- wusagfitted[n,1]
  tmpslp <- wusagfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[2])
  points(2, wtmp, pch=16, col=cols[2])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[2], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_sag.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- incsagfitted[n,1]
  tmpslp <- incsagfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[2])
  points(2, wtmp, pch=16, col=cols[2])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[2], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_sag.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- svlsagfitted[n,1]
  tmpslp <- svlsagfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[2])
  points(2, wtmp, pch=16, col=cols[2])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[2], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_sag.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- sssagfitted[n,1]
  tmpslp <- sssagfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[2])
  points(2, wtmp, pch=16, col=cols[2])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[2], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_sag.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- endsagfitted[n,1]
  tmpslp <- endsagfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[2])
  points(2, wtmp, pch=16, col=cols[2])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[2], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# carolinensis #
wucarfitted <- wufittedrns[c("Cage.21", "Cage.22", "Cage.23", "Cage.24", "Cage.25", "Cage.26", "Cage.27", "Cage.28","Cage.29", "Cage.30"),]
inccarfitted <- incfittedrns[c("Cage.21", "Cage.22", "Cage.23", "Cage.24", "Cage.25", "Cage.26", "Cage.27", "Cage.28","Cage.29", "Cage.30"),]
svlcarfitted <- svlfittedrns[c("Cage.21", "Cage.22", "Cage.23", "Cage.24", "Cage.25", "Cage.26", "Cage.27", "Cage.28","Cage.29", "Cage.30"),]
sscarfitted <- ssfittedrns[c("Cage.21", "Cage.22", "Cage.23", "Cage.24", "Cage.25", "Cage.26", "Cage.27", "Cage.28","Cage.29", "Cage.30"),]
endcarfitted <- endfittedrns[c("Cage.21", "Cage.22", "Cage.23", "Cage.24", "Cage.25", "Cage.26", "Cage.27", "Cage.28","Cage.29", "Cage.30"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_car.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- wucarfitted[n,1]
  tmpslp <- wucarfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[3])
  points(2, wtmp, pch=16, col=cols[3])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[3], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_car.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- inccarfitted[n,1]
  tmpslp <- inccarfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[3])
  points(2, wtmp, pch=16, col=cols[3])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[3], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_car.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- svlcarfitted[n,1]
  tmpslp <- svlcarfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[3])
  points(2, wtmp, pch=16, col=cols[3])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[3], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_car.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- sscarfitted[n,1]
  tmpslp <- sscarfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[3])
  points(2, wtmp, pch=16, col=cols[3])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[3], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_car.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- endcarfitted[n,1]
  tmpslp <- endcarfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[3])
  points(2, wtmp, pch=16, col=cols[3])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[3], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# chlorocyanus #

wuchlfitted <- wufittedrns[c("Cage.31", "Cage.32", "Cage.33", "Cage.34", "Cage.35", "Cage.36", "Cage.37", "Cage.38","Cage.40", "Cage.42"),]
incchlfitted <- incfittedrns[c("Cage.31", "Cage.32", "Cage.33", "Cage.34", "Cage.35", "Cage.36", "Cage.37", "Cage.38","Cage.40", "Cage.42"),]
svlchlfitted <- svlfittedrns[c("Cage.31", "Cage.32", "Cage.33", "Cage.34", "Cage.35", "Cage.36", "Cage.37", "Cage.38","Cage.40", "Cage.42"),]
sschlfitted <- ssfittedrns[c("Cage.31", "Cage.32", "Cage.33", "Cage.34", "Cage.35", "Cage.36", "Cage.37", "Cage.38","Cage.40", "Cage.42"),]
endchlfitted <- endfittedrns[c("Cage.31", "Cage.32", "Cage.33", "Cage.34", "Cage.35", "Cage.36", "Cage.37", "Cage.38","Cage.40", "Cage.42"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_chl.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- wuchlfitted[n,1]
  tmpslp <- wuchlfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[4])
  points(2, wtmp, pch=16, col=cols[4])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[4], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_chl.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- incchlfitted[n,1]
  tmpslp <- incchlfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[4])
  points(2, wtmp, pch=16, col=cols[4])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[4], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_chl.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- svlchlfitted[n,1]
  tmpslp <- svlchlfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[4])
  points(2, wtmp, pch=16, col=cols[4])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[4], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_chl.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- sschlfitted[n,1]
  tmpslp <- sschlfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[4])
  points(2, wtmp, pch=16, col=cols[4])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[4], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_chl.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- endchlfitted[n,1]
  tmpslp <- endchlfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[4])
  points(2, wtmp, pch=16, col=cols[4])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[4], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# distichus #

wudistfitted <- wufittedrns[c("Cage.43", "Cage.44", "Cage.45", "Cage.46", "Cage.47", "Cage.48", "Cage.49", "Cage.50","Cage.51", "Cage.52"),]
incdistfitted <- incfittedrns[c("Cage.43", "Cage.44", "Cage.45", "Cage.46", "Cage.47", "Cage.48", "Cage.49", "Cage.50","Cage.51", "Cage.52"),]
svldistfitted <- svlfittedrns[c("Cage.43", "Cage.44", "Cage.45", "Cage.46", "Cage.47", "Cage.48", "Cage.49", "Cage.50","Cage.51", "Cage.52"),]
ssdistfitted <- ssfittedrns[c("Cage.43", "Cage.44", "Cage.45", "Cage.46", "Cage.47", "Cage.48", "Cage.49", "Cage.50","Cage.51"),]
enddistfitted <- endfittedrns[c("Cage.43", "Cage.44", "Cage.45", "Cage.46", "Cage.47", "Cage.48", "Cage.49", "Cage.50","Cage.51"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_dist.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- wudistfitted[n,1]
  tmpslp <- wudistfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[5])
  points(2, wtmp, pch=16, col=cols[5])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[5], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_dist.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- incdistfitted[n,1]
  tmpslp <- incdistfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[5])
  points(2, wtmp, pch=16, col=cols[5])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[5], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_dist.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- svldistfitted[n,1]
  tmpslp <- svldistfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[5])
  points(2, wtmp, pch=16, col=cols[5])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[5], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_dist.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- ssdistfitted[n,1]
  tmpslp <- ssdistfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[5])
  points(2, wtmp, pch=16, col=cols[5])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[5], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_dist.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 10){
  tmpint <- enddistfitted[n,1]
  tmpslp <- enddistfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[5])
  points(2, wtmp, pch=16, col=cols[5])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[5], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# cybotes #

wucybfitted <- wufittedrns[c("Cage.57", "Cage.58", "Cage.59", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
inccybfitted <- incfittedrns[c("Cage.57", "Cage.58", "Cage.59", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
svlcybfitted <- svlfittedrns[c("Cage.57", "Cage.58", "Cage.59", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
sscybfitted <- ssfittedrns[c("Cage.57", "Cage.58", "Cage.59", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
endcybfitted <- endfittedrns[c("Cage.57", "Cage.58", "Cage.59", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_cyb.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.65, 0.5)) #-0.16, 0.16
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 12){
  tmpint <- wucybfitted[n,1]
  tmpslp <- wucybfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[6])
  points(2, wtmp, pch=16, col=cols[6])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[6], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_cyb.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2, 1.7))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 12){
  tmpint <- inccybfitted[n,1]
  tmpslp <- inccybfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[6])
  points(2, wtmp, pch=16, col=cols[6])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[6], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_cyb.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-2.3, 1.5))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 12){
  tmpint <- svlcybfitted[n,1]
  tmpslp <- svlcybfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[6])
  points(2, wtmp, pch=16, col=cols[6])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[6], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_cyb.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 12){
  tmpint <- sscybfitted[n,1]
  tmpslp <- sscybfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[6])
  points(2, wtmp, pch=16, col=cols[6])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[6], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_cyb.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-3, 4))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 12){
  tmpint <- endcybfitted[n,1]
  tmpslp <- endcybfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[6])
  points(2, wtmp, pch=16, col=cols[6])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[6], lwd=10)
 # lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

# equestris #

wueqfitted <- wufittedrns[c("Cage.68", "Cage.69", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
inceqfitted <- incfittedrns[c("Cage.68", "Cage.69", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
svleqfitted <- svlfittedrns[c("Cage.68", "Cage.69", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
sseqfitted <- ssfittedrns[c("Cage.68", "Cage.69", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]
endeqfitted <- endfittedrns[c("Cage.68", "Cage.69", "Cage.60", "Cage.61", "Cage.62", "Cage.63", "Cage.64","Cage.65", "Cage.66", "Cage.67"),]

svg(filename="R_exports/intraspecific_fig/log_scaled/wu/wulog_eq.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.2))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.65,0.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- wueqfitted[n,1]
  tmpslp <- wueqfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[7])
  points(2, wtmp, pch=16, col=cols[7])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[7], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#incubation period

svg(filename="R_exports/intraspecific_fig/log_scaled/inc/inclog_eq.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.2))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2,1.7), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- inceqfitted[n,1]
  tmpslp <- inceqfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[7])
  points(2, wtmp, pch=16, col=cols[7])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[7], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#svl
svg(filename="R_exports/intraspecific_fig/log_scaled/SVL/svllog_eq.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.2))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-2.3,1.5), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- svleqfitted[n,1]
  tmpslp <- svleqfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[7])
  points(2, wtmp, pch=16, col=cols[7])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[7], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#sprint speed
svg(filename="R_exports/intraspecific_fig/log_scaled/sprint-speed/sslog_eq.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.2))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-0.2,0.4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- sseqfitted[n,1]
  tmpslp <- sseqfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[7])
  points(2, wtmp, pch=16, col=cols[7])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[7], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()

#endurance
svg(filename="R_exports/intraspecific_fig/log_scaled/endurance/endlog_eq.svg")
plot.new()
plot.window(xlim=c(0.95,2.05), ylim=c(-0.2, 0.2))
axis(1, at=c(1,2), labels=NA, lwd=10)
axis(2, at=c(-3,4), labels=NA, lwd=10)
n <- 1
while (n < 11){
  tmpint <- endeqfitted[n,1]
  tmpslp <- endeqfitted[n,2]
  ctmp <- (tmpint - tmpslp)
  wtmp <- (tmpint + tmpslp)
  points(1, ctmp, pch=16, col=cols[7])
  points(2, wtmp, pch=16, col=cols[7])
  lines(c(ctmp,wtmp) ~ c(1,2), col=cols[7], lwd=10)
  #lines(c(0,0) ~ c(1,2), col="black", lwd=10, lty="dashed")
  n <- n+1
}
dev.off()




## Make interspecific panels
#Water Uptake


#Incubation Period


#SVL


#Sprint Speed


#Endurance


## Assemble full multipanel plot


### Heatmap of island environments




# (done to enable variable data trimming and transformation based on phenotype)

##Water Uptake
#fixed effects: egg mass, % developed, Treatment, incubator, lay day
#random effects: Cage, Species
wutmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                   "EggMass","scEggMass","PercDev","WaterUptake","scWU")]
#wu <- wutmp[!is.na(wutmp$WaterUptake),]
wu <- wutmp[which(!is.na(wutmp$WaterUptake) & !is.na(wutmp$PercDev)),] #do instead, alter for all variables
rm(wutmp)

##Incubation period
#fixed effects: egg mass, incubator, Treatment, lay day
#random effects: Cage, species
incptmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                     "EggMass","scEggMass","Inc")]
incp <- incptmp[which(!is.na(incptmp$Inc) & !is.na(incptmp$EggMass)),] #do instead, alter for all variables

rm(incptmp)

##Hatchling morphology (SVL, TL, mass)
#fixed effects: egg mass, treatment, other morphology, incubator, Treatment, lay day
#random effects: Cage, Species
tmphatchmorph <- mializ[,c("ID", "Species", "Cage", "nTreatment", "Incubator", "LayDay",
                           "EggMass","scEggMass","SVL", "TL", "Mass", "Sex")]
hatchmorph <- tmphatchmorph[which(!is.na(tmphatchmorph$SVL) & !is.na(tmphatchmorph$EggMass)),]
rm(tmphatchmorph)

##Sprint Speed
#fixed effects: SVL, SprintTemp, Treatment, Age, Stops, Fails, Sex
#random effects: Cage, Species
sstmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "LayDay",
                   "SVL", "Sex", "FastestSprint", "Stops", "Fails", "SprintTemp",
                   "AgeSprint")]
ss <- sstmp[which(!is.na(sstmp$FastestSprint) & !is.na(sstmp$SVL) & !is.na(sstmp$Stops) &
                    !is.na(sstmp$Fails) & !is.na(sstmp$SprintTemp) & !is.na(sstmp$AgeSprint)),]
rm(sstmp)

##Endurance
#fixed effects: SVL, Fastest sprint?, Treatment, Age, include both velocities, Sex
#random effects: Cage, Species
endurtmp <- mializ[,c("ID", "Species", "Cage", "nTreatment", "LayDay",
                      "SVL", "Sex", "FastestSprint", "Velocity_m.s", "Velocity_cm.s",
                      "DistanceRun_cm", "EnduranceTime", "AgeEndurance")]
endur <- endurtmp[which(!is.na(endurtmp$SVL) & !is.na(endurtmp$FastestSprint) & !is.na(endurtmp$Velocity_cm.s) &
                          !is.na(endurtmp$DistanceRun_cm) & !is.na(endurtmp$EnduranceTime) & !is.na(endurtmp$AgeEndurance)),]
rm(endurtmp)





#blah blah 

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

