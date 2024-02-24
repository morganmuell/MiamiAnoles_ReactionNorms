##--Load libraries
library(phytools)
library(MCMCglmm)
library(nadiv)

##--Prep Data and Priors
#load(insertdata.Rdata)
load("model2df.rda")

tree <- read.tree(file="R_imports/pruned_anole_tree.treefile")

Ainv <- inverseA(tree)$Ainv

prior1 <- list(R = list(V = 1e-12, nu = -2), #<-- Non-informative improper: *marginal* posterior equal to REML estimate
               G = list(G1 = list(V = diag(2)*0.002, nu = 3, 
                                  alpha.mu = c(0,0), alpha.V = diag(2)*10000)))
pr2A <- list(R = list(V = 1, nu = 0.002), 
             G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10), 
                      G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10))) 

pr2B <- list(R = list(V = 1, nu = 0.002), 
             G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10)))

##-- Run Model 1
set.seed(330)
bwumodel1 <- MCMCglmm(logWU ~ 1 + EggMass + nTreatment + PercDev + Species
                      + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = wu, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = TRUE, prior = prior1)

##--Housekeeping prior to looop
postsize <- nrow(bwumodel1$VCV) 

#output dataframes for model effect sizes
preenveff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion time: iteration, post.mean, two slots for confidence interval, eff.size, p-value
colnames(preenveff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
invenveff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion environment
colnames(invenveff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
invtimeeff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion time
colnames(invtimeeff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
inteff <- matrix(NA, nrow=postsize, ncol=6) #effects for intercept
colnames(inteff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")

phyleff <- matrix(NA, nrow=postsize, ncol=5) #effects for species random effect: iteration, postmean, two CI slots, effsamp
colnames(phyleff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp")
sppeff <- matrix(NA, nrow=postsize, ncol=5) #effects for species random effect: iteration, postmean, two CI slots, effsamp
colnames(sppeff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp")
unitseff <- matrix(NA, nrow=postsize, ncol=5) #effects for residuals: iteration, postmean, two CI slots, effsamp
colnames(unitseff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp")

polypreenveff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion time: iteration, post.mean, two slots for confidence interval, eff.size, p-value
colnames(polypreenveff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
polyinvenveff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion environment
colnames(polyinvenveff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
polyinvtimeeff <- matrix(NA, nrow=postsize, ncol=6) #effects for invasion time
colnames(polyinvtimeeff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
polyinteff <- matrix(NA, nrow=postsize, ncol=6) #effects for intercept
colnames(polyinteff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp", "pMCMC")
polysppeff <- matrix(NA, nrow=postsize, ncol=5) #effects for species random effect: iteration, postmean, two CI slots, effsamp
colnames(polysppeff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp")
polyunitseff <- matrix(NA, nrow=postsize, ncol=5) #effects for residuals: iteration, postmean, two CI slots, effsamp
colnames(polyunitseff) <- c("Iteration", "post.mean", "lower.95%.CI", "upper.95%.CI", "eff.samp")

#posterior modes
postfmodes <- matrix(NA, nrow=postsize, ncol=5) #posterior modes for fixed effects: iteration, intercept, preclim, invclim, invtime
colnames(postfmodes) <- c("Iteration", "Intercept", "Preclim", "InvClim", "InvTime")
postsppmodes <- matrix(NA, nrow=postsize, ncol=4) #posterior modes for random effect: iteration, species, units
colnames(postsppmodes) <- c("Iteration", "Species","AbSpec", "Units")

polypostfmodes <- matrix(NA, nrow=postsize, ncol=5) #posterior modes for fixed effects: iteration, intercept, preclim, invclim, invtime
colnames(polypostfmodes) <- c("Iteration", "Intercept", "Preclim", "InvClim", "InvTime")
polypostsppmodes <- matrix(NA, nrow=postsize, ncol=3) #posterior modes for random effect: iteration, species, units
colnames(polypostsppmodes) <- c("Iteration", "Species", "Units")

#and convergence diagnostics
heidelfout <- matrix(NA, nrow=postsize, ncol=17) #p-values for stationarity, fixed effects: iteration, passint, pint, passpre, ppre, passinv, pinv, inthalfpass,, prepasshalf, invpasshalf
colnames(heidelfout) <- c("Iteration", "Pass.Stat.Int", "Pvalue.Int", "Pass.Preclim", "Pvalue.Preclim",
                          "Pass.InvClim", "Pvalue.InvClim", "Pass.InvTime", "Pvalue.InvTime", "Pass.Halfwidth.Int", "Ratio.Int.",
                          "Pass.Halfwidth.Preclim", "Ratio.Preclim", "Pass.Halfwidth.InvClim",
                          "Ratio.InvClim", "Pass.Halfwidth.InvTime", "Ratio.InvTime")
heidelrout <- matrix(NA, nrow=postsize, ncol=13) #p-values for stationarity, random effect: iteration, spppass, sppp, unitspass, unitsp, spphalfpass, unitshalfpass
colnames(heidelrout) <- c("Iteration", "Pass.Stat.Spp", "Pvalue.Spp","Pass.Stat.AbSpec","Pvalue.AbSpec", "Pass.Stat.Resid", 
                          "Pvalue.Resid", "Pass.Halfwidth.Spp", "Ratio.Spp", "Pass.Halfwidth.AbSpec","Ratio.AbSpec","Pass.Halfwidth.Resid",
                          "Ratio.Resid")

polyheidelfout <- matrix(NA, nrow=postsize, ncol=17) #p-values for stationarity, fixed effects: iteration, passint, pint, passpre, ppre, passinv, pinv, inthalfpass,, prepasshalf, invpasshalf
colnames(polyheidelfout) <- c("Iteration", "Pass.Stat.Int", "Pvalue.Int", "Pass.Preclim", "Pvalue.Preclim",
                              "Pass.InvClim", "Pvalue.InvClim", "Pass.InvTime", "Pvalue.InvTime", "Pass.Halfwidth.Int", "Ratio.Int.",
                              "Pass.Halfwidth.Preclim", "Ratio.Preclim", "Pass.Halfwidth.InvClim",
                              "Ratio.InvClim", "Pass.Halfwidth.InvTime", "Ratio.InvTime")
polyheidelrout <- matrix(NA, nrow=postsize, ncol=9) #p-values for stationarity, random effect: iteration, spppass, sppp, unitspass, unitsp, spphalfpass, unitshalfpass
colnames(polyheidelrout) <- c("Iteration", "Pass.Stat.Spp", "Pvalue.Spp", "Pass.Stat.Resid", 
                              "Pvalue.Resid", "Pass.Halfwidth.Spp", "Ratio.Spp", "Pass.Halfwidth.Resid",
                              "Ratio.Resid")
autofix <- matrix(NA, nrow=postsize, ncol=5) #mean autocorr for fixed effects
colnames(autofix) <- c("Iteration", "Intercept", "scPreClim", "scInvClim","InvTime")
autorand <- matrix(NA, nrow=postsize, ncol=4) #mean autocorr for random effects
colnames(autorand) <- c("Iteration", "Species", "AbSpec", "units")

polyautofix <- matrix(NA, nrow=postsize, ncol=5) #mean autocorr for fixed effects
colnames(polyautofix) <- c("Iteration", "Intercept", "scPreClim", "scInvClim","InvTime")
polyautorand <- matrix(NA, nrow=postsize, ncol=3) #mean autocorr for random effects
colnames(polyautorand) <- c("Iteration", "AbSpec", "units")

dicscore <- matrix(NA, nrow=postsize, ncol=15)
colnames(dicscore) <- c("Iteration","PhyloDIC","PolyDIC", "Diff")

hpd <- matrix(NA, nrow=postsize, ncol=15)
colnames(hpd) <- c("Iteration", "lowerInt" "upperInt", "lPreClim", "uPreClim", "lInvClim", "uInvClim","lInvTime","uInvTime",
						"lSpecies", "uSpecies", "lAbSpec", "uAbSpec", "lUnits", "uUnits")
polyhpd <- matrix(NA, nrow=postsize, ncol=13)
colnames(polyhpd) <- c("Iteration", "lowerInt" "upperInt", "lPreClim", "uPreClim", "lInvClim", "uInvClim","lInvTime","uInvTime",
						"lAbSpec", "uAbSpec", "lUnits", "uUnits")


#loop management
i <- 1
j <- 1 #use i to index run in model 1, use j to index output files

##--NOW GO!
while(i<postsize){
  print(j)
  
  ## Pull slopes from model iteration
  tmpMod <- bwumodel1
  SolColNms <- dimnames(tmpMod$Sol)[[2L]] #pulls random effects only
  slpColInd <- which(startsWith(SolColNms, "nTreatment.Cage.")) #cuts off inconvenient names just to grab the column index number
  devs <- tmpMod$Sol[i, slpColInd] #these are deviations of cage from species x treatment average
  lstNmWrds <- strsplit(x = names(devs), split = ".", fixed = TRUE) 
  cageNumbs <- sapply(lstNmWrds, FUN = "[", i = 3) #this returns vector of only cage numbers, in order
  islps <- matrix(NA, nrow=length(devs), ncol=1)
  slopeframe <- data.frame(cageNumbs, devs, islps)
  
  #convert cage deviations to deviations from model average -- CageDev + SpecAverage + RefSpp
  cg <- 1 #index row number for loop
  while(cg < length(devs)+1){
    CID <- slopeframe[cg,"cageNumbs"] #make variable for cage ID
    if(CID<11){ #cristatellus
      slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Speciescristatellus"]
    } else { #if the species is a green anole (range 11-20)
      if(CID>10 | CID<21){ #sagrei
        slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Speciessagrei"]
      } else {
        if(CID>20 | CID<31){ #carolinensis
          slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"]
        } else {
          if(CID>30 | CID<43){ #chlorocyanus
            slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Specieschlorocyanus"]
          } else {
            if(CID>42 | CID<53){ #distichus
              slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Speciesdistichus"]
            } else {
              if(CID>56| CID<68){ #cybotes
                slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Speciescybotes"]
              } else {
                if(CID>67| CID<78){ #equestris
                  slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + bwumodel1$Sol[i,"nTreatment"] + bwumodel1$Sol[i,"nTreatment:Speciesequestris"]
                } else {slopeframe[cg,"islps"] <- "ERROR: CAGE NUMBER NOT FOUND IN DATASET"}
              }
            }
          }
        }
      }
    }
    cg <- cg + 1
  }
  
  ordslp <- slopeframe[order(as.integer(slopeframe$cageNumbs)),]
  iframe <- model2df[(model2df$Cage %in% cageNumbs),]
  iframe$islps <- ordslp$islps[match(iframe$Cage, ordslp$cageNumbs)]
  iframe$scInvTime[is.na(iframe$scInvTime)] <- 0 #<-- trick MCMCglmm that no missing values in InvTime
  
  ## Run Model 2 with phylogeny 
  set.seed(330)
  bwumodel2 <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                        random=~Species + AbSpec,
                        ginverse=list(Species=Ainv), 
                        data=iframe,
                        prior=pr2A,
                        family = rep("gaussian", 1),
                        nitt=500000, burnin=50000, thin=250, verbose=FALSE)
  
  
  
  ## Run Model 2 without phylogeny
  set.seed(330)
  bwumodel2poly <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                            random=~AbSpec,
                            #ginverse=list(Species=Ainv), 
                            data=iframe,
                            prior=pr2B,
                            family = rep("gaussian", 1),
                            nitt=500000, burnin=50000, thin=250, verbose=FALSE)
  
  
  #grab relevant output
  tmpsum <- summary(bwumodel2)
  tmpheifx <- heidel.diag(bwumodel2$Sol)
  tmpheiran <- heidel.diag(bwumodel2$VCV)
  tmpsumpoly <- summary(bwumodel2poly)
  tmpheifxpoly <- heidel.diag(bwumodel2poly$Sol)
  tmpheiranpoly <- heidel.diag(bwumodel2poly$VCV)
  tmpacrand <- autocorr.diag(bwumodel2$VCV)[2,]
  tmpacrandpoly <- autocorr.diag(bwumodel2poly$VCV)[2,]
  tmpacfix <- autocorr.diag(bwumodel2$Sol)[2,]
  tmpacfixpoly <- autocorr.diag(bwumodel2poly$Sol)[2,]
  tmphpd <- rbind(HPDinterval(bwumodel2$Sol), HPDinterval(bwumodel2$VCV))
  tmphpdpoly <- rbind(HPDinterval(bwumodel2poly$Sol), HPDinterval(bwumodel2poly$VCV))
  
  ##--pull effects and convergence diagnostics and store 
  hpd[j,1] <- i
  hpd[j,2] <- tmphpd[1,"lower"]
  hpd[j,3] <- tmphpd[1,"upper"]
  hpd[j,4] <- tmphpd[2,"lower"]
  hpd[j,5] <- tmphpd[2,"upper"]
  hpd[j,6] <- tmphpd[3,"lower"]
  hpd[j,7] <- tmphpd[3,"upper"]
  hpd[j,8] <- tmphpd[4,"lower"]
  hpd[j,9] <- tmphpd[4,"upper"]
  hpd[j,10] <- tmphpd["Species","lower"]
  hpd[j,11] <- tmphpd["Species", "upper"]
  hpd[j,12] <- tmphpd["AbSpec","lower"]
  hpd[j,13] <- tmphpd["AbSpec","upper"]
  hpd[j,14] <- tmphpd["units","lower"]
  hpd[j,15] <- tmphpd["units","upper"]
  
  polyhpd[j,1] <- i
  polyhpd[j,2] <- tmphpdpoly[1,"lower"]
  polyhpd[j,3] <- tmphpdpoly[1,"upper"]
  polyhpd[j,4] <- tmphpdpoly[2,"lower"]
  polyhpd[j,5] <- tmphpdpoly[2,"upper"]
  polyhpd[j,6] <- tmphpdpoly[3,"lower"]
  polyhpd[j,7] <- tmphpdpoly[3,"upper"]
  polyhpd[j,8] <- tmphpdpoly[4,"lower"]
  polyhpd[j,9] <- tmphpdpoly[4,"upper"]
  polyhpd[j,10] <- tmphpdpoly[5,"lower"]
  polyhpd[j,11] <- tmphpdpoly[5,"upper"]
  polyhpd[j,12] <- tmphpdpoly[6,"lower"]
  polyhpd[j,13] <- <- tmphpdpoly[6,"upper"]
  
  #DIC scores
  dicscore[j,1] <- i
  dicscore[j,2] <- tmpsum$DIC
  dicscore[j,3] <- tmpsumpoly$DIC
  dicscore[j,4] <- tmpsum$DIC - tmpsumpoly$DIC
  
  #autocorrelation
  autofix[j,1] <- i
  autofix[j,c(2:5)] <- tmpacfix[c(1:4)]
  autorand[j,1] <- i
  autorand[j,c(2:4)] <- tmpacrand[c(1:3)]
  
  polyautofix[j,1] <- i
  polyautofix[j,c(2:5)] <- tmpacfixpoly[c(1:4)]
  polyautorand[j,1] <- i
  polyautorand[j,c(2:3)] <- tmpacrandpoly[c(1:2)]
  
  #effects
  phyleff[j,1] <- i
  phyleff[j,c(2:5)] <- tmpsum$Gcovariances[1,]
  sppeff[j,1] <- i
  sppeff[j,c(2:5)] <- tmpsum$Gcovariances[2,]
  unitseff[j,1] <- i
  unitseff[j,c(2:5)] <- tmpsum$Rcovariances
  inteff[j,1] <- i
  inteff[j,c(2:6)] <- tmpsum$solutions[1,]
  preenveff[j,1] <- i
  preenveff[j,c(2:6)] <- tmpsum$solutions[2,]
  invenveff[j,1] <- i
  invenveff[j,c(2:6)] <- tmpsum$solutions[3,]
  invtimeeff[j,1] <- i
  invtimeeff[j,c(2:6)] <- tmpsum$solutions[4,]
  
  polysppeff[j,1] <- i
  polysppeff[j,c(2:5)] <- tmpsumpoly$Gcovariances
  polyunitseff[j,1] <- i
  polyunitseff[j,c(2:5)] <- tmpsumpoly$Rcovariances
  polyinteff[j,1] <- i
  polyinteff[j,c(2:6)] <- tmpsumpoly$solutions[1,]
  polypreenveff[j,1] <- i
  polypreenveff[j,c(2:6)] <- tmpsumpoly$solutions[2,]
  polyinvenveff[j,1] <- i
  polyinvenveff[j,c(2:6)] <- tmpsumpoly$solutions[3,]
  polyinvtimeeff[j,1] <- i
  polyinvtimeeff[j,c(2:6)] <- tmpsumpoly$solutions[4,]
  
  #posterior modes
  postfmodes[j,1] <- i
  postfmodes[j,c(2:5)] <- posterior.mode(bwumodel2$Sol)
  postsppmodes[j,1] <- i
  postsppmodes[j,c(2:4)] <- posterior.mode(bwumodel2$VCV)
  
  polypostfmodes[j,1] <- i
  polypostfmodes[j,c(2:5)] <- posterior.mode(bwumodel2poly$Sol)
  polypostsppmodes[j,1] <- i
  polypostsppmodes[j,c(2:3)] <- posterior.mode(bwumodel2poly$VCV)
  
  #convergence diagnostics
  ##fixed effects
  heidelfout[j,1] <- i
  heidelfout[j,2] <- tmpheifx[1,1] #intercept pass/fail, stationarity
  heidelfout[j,3] <- tmpheifx[1,3] #intercept p-value, stationarity
  heidelfout[j,4] <- tmpheifx[2,1] #preclim pass/fail, stationarity
  heidelfout[j,5] <- tmpheifx[2,3] #preclim p-value, stationarity
  heidelfout[j,6] <- tmpheifx[3,1] #postclim pass/fail, stationarity
  heidelfout[j,7] <- tmpheifx[3,3] #postclim p-value, stationarity
  heidelfout[j,8] <- tmpheifx[4,1] #invTime pass/fail, stationarity
  heidelfout[j,9] <- tmpheifx[4,3] #invtime p-value, stationarity
  heidelfout[j,10] <- tmpheifx[1,4] #intercept pass/fail, halfwidth
  heidelfout[j,11] <- (tmpheifx[1,6] / tmpheifx[1,5]) #intercept ratio; value must be below 0.1 to pass
  heidelfout[j,12] <- tmpheifx[2,4] #preclim pass/fail, halfwidth
  heidelfout[j,13] <- (tmpheifx[2,6] / tmpheifx[2,5]) #preclim ratio
  heidelfout[j,14] <- tmpheifx[3,4] #postclim pass/fail
  heidelfout[j,15] <- (tmpheifx[3,6] / tmpheifx[3,5]) #postclim ratio 
  heidelfout[j,16] <- tmpheifx[4,4] #invtime pass/fail
  heidelfout[j,17] <- (tmpheifx[4,6] / tmpheifx[4,5]) #invtime ratio
  
  polyheidelfout[j,1] <- i
  polyheidelfout[j,2] <- tmpheifxpoly[1,1] #intercept pass/fail, stationarity
  polyheidelfout[j,3] <- tmpheifxpoly[1,3] #intercept p-value, stationarity
  polyheidelfout[j,4] <- tmpheifxpoly[2,1] #preclim pass/fail, stationarity
  polyheidelfout[j,5] <- tmpheifxpoly[2,3] #preclim p-value, stationarity
  polyheidelfout[j,6] <- tmpheifxpoly[3,1] #postclim pass/fail, stationarity
  polyheidelfout[j,7] <- tmpheifxpoly[3,3] #postclim p-value, stationarity
  polyheidelfout[j,8] <- tmpheifxpoly[4,1] #invTime pass/fail, stationarity
  polyheidelfout[j,9] <- tmpheifxpoly[4,3] #invtime p-value, stationarity
  polyheidelfout[j,10] <- tmpheifxpoly[1,4] #intercept pass/fail, halfwidth
  polyheidelfout[j,11] <- (tmpheifxpoly[1,6] / tmpheifxpoly[1,5]) #intercept ratio; value must be below 0.1 to pass
  polyheidelfout[j,12] <- tmpheifxpoly[2,4] #preclim pass/fail, halfwidth
  polyheidelfout[j,13] <- (tmpheifxpoly[2,6] / tmpheifxpoly[2,5]) #preclim ratio
  polyheidelfout[j,14] <- tmpheifxpoly[3,4] #postclim pass/fail
  polyheidelfout[j,15] <- (tmpheifxpoly[3,6] / tmpheifxpoly[3,5]) #postclim ratio 
  polyheidelfout[j,16] <- tmpheifxpoly[4,4] #invtime pass/fail
  polyheidelfout[j,17] <- (tmpheifxpoly[4,6] / tmpheifxpoly[4,5]) #invtime ratio
  
  ##random effect and residuals
  heidelrout[j,1] <- i
  heidelrout[j,2] <- tmpheiran[1,1] #species pass/fail, stationarity
  heidelrout[j,3] <- tmpheiran[1,3] #species p-value, stationarity
  heidelrout[j,4] <- tmpheiran[2,1] #absec pass/fail, stationarity
  heidelrout[j,5] <- tmpheiran[2,3] #abspec p-value, stationarity
  heidelrout[j,6] <- tmpheiran[3,1] #residuals pass/fail, stationarity
  heidelrout[j,7] <- tmpheiran[3,3] #residuals p-value, stationarity
  heidelrout[j,8] <- tmpheiran[1,4] #species pass/fail, halfwidth
  heidelrout[j,9] <- (tmpheiran[1,6] / tmpheiran[1,5]) #species ratio
  heidelrout[j,10] <- tmpheiran[2,4] #absec pass/fail, halfwidth
  heidelrout[j,11] <- (tmpheiran[2,6] / tmpheiran[2,5]) #abspec ratio)
  heidelrout[j,12] <- tmpheiran[3,4] #residuals pass/fail, halfwidth
  heidelrout[j,13] <- (tmpheiran[3,6] / tmpheiran[3,5]) #residuals ratio
  
  polyheidelrout[j,1] <- i
  polyheidelrout[j,2] <- tmpheiranpoly[1,1] #species pass/fail, stationarity
  polyheidelrout[j,3] <- tmpheiranpoly[1,3] #species p-value, stationarity
  polyheidelrout[j,4] <- tmpheiranpoly[2,1] #residuals pass/fail, stationarity
  polyheidelrout[j,5] <- tmpheiranpoly[2,3] #residuals p-value, stationarity
  polyheidelrout[j,6] <- tmpheiranpoly[1,4] #species pass/fail, halfwidth
  polyheidelrout[j,7] <- (tmpheiranpoly[1,6] / tmpheiranpoly[1,5]) #species ratio
  polyheidelrout[j,8] <- tmpheiranpoly[2,4] #residuals pass/fail, halfwidth
  polyheidelrout[j,9] <- (tmpheiranpoly[2,6] / tmpheiranpoly[2,5]) #residuals ratio
  
  #last but not least
  i<- i+50
  j<- j+1
  
}

##--export results
#Model 1 Output
wumodel1out <- summary(bwumodel1)
sink("pheno_summ_model1.txt")
print(wumodel1out)
sink()
wumode1 <- posterior.mode(bwumodel1$VCV)
write.csv(wumode1, file="pheno_mode_1.csv")
wuautocor1 <- autocorr.diag(bwumodel1$VCV)
write.csv(wuautocor1, file="pheno_autocorr_VCV1.csv")
wuheidelrout1 <- heidel.diag(bwumodel1$VCV)
write.csv(wuheidelrout1, file="pheno_heidelrout_1.csv")
wuheidelfout1 <- heidel.diag(bwumodel1$Sol)
write.csv(wuheidelfout1, file="pheno_heidelfout_1.csv")
wuHPD1 <- rbind(HPDinterval(bwumodel1$VCV), HPDinterval(bwumodel1$Sol))
write.csv(wuHPD1, file="pheno_HPD_1.csv")


#Model 2 Output
write.csv(sppeff, file="pheno_sppeff_longrun.csv")
write.csv(unitseff, file="pheno_unitseff_longrun.csv")
write.csv(inteff, file="pheno_inteff_longrun.csv")
write.csv(preenveff, file="pheno_preenveff_longrun.csv")
write.csv(invenveff, file="pheno_invenveff_longrun.csv")
write.csv(invtimeeff, file="pheno_invtimeeff_longrun.csv")
write.csv(postfmodes, file="pheno_postfmodes_longrun.csv")
write.csv(postsppmodes, file="pheno_postsppmodes_longrun.csv")
write.csv(heidelfout, file="pheno_heidelfout_longrun.csv")
write.csv(heidelrout, file="pheno_heidelrout_longrun.csv")

write.csv(polysppeff, file="poly_pheno_sppeff_longrun.csv")
write.csv(polyunitseff, file="poly_pheno_unitseff_longrun.csv")
write.csv(polyinteff, file="poly_pheno_inteff_longrun.csv")
write.csv(polypreenveff, file="poly_pheno_preenveff_longrun.csv")
write.csv(polyinvenveff, file="poly_pheno_invenveff_longrun.csv")
write.csv(polyinvtimeeff, file="poly_pheno_invtimeeff_longrun.csv")
write.csv(polypostfmodes, file="poly_pheno_postfmodes_longrun.csv")
write.csv(polypostsppmodes, file="poly_pheno_postsppmodes_longrun.csv")
write.csv(polyheidelfout, file="poly_pheno_heidelfout_longrun.csv")
write.csv(polyheidelrout, file="poly_pheno_heidelrout_longrun.csv")


