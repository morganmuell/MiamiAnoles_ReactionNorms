##--Load libraries
library(MCMCglmm)
library(nadiv)
library(phytools)

##--Prep Data and Priors
load("inctrim.Rdata")
load("../model2df.rda")

tree <- read.tree(file="../pruned_anole_tree.treefile")

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
bincmodel1 <- MCMCglmm(logInc ~ 1 + EggMass + nTreatment + Species
                      + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = incp, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = FALSE, prior = prior1)

print("Model 1 Complete")

##--Housekeeping prior to looop
postsize <- 100 

# Output mcmc object for all model names
incmod2.mcmc <- as.data.frame(matrix(NA, nrow=200000, ncol=12))
colnames(incmod2.mcmc) <- c("Intercept","PreClim","InvClim","InvTime","Phylo","Spp","units","h2phylo","s2","CVphylo","CVspp","CVresid")
incmod2poly.mcmc <- as.data.frame(matrix(NA, nrow=200000, ncol=9))
colnames(incmod2poly.mcmc) <- c("Intercept","PreClim","InvClim","InvTime","Spp","units","s2","CVspp","CVresid")

# Create new model summary frame
#full model
m2sum <- as.data.frame(matrix(NA, nrow=12, ncol=5))
rownames(m2sum) <- c("Intercept","PreClim","InvClim","InvTime","Phylo","Spp","units","h2phylo","s2","CVphylo","CVspp","CVresid")
colnames(m2sum) <- c("post.mean", "lower.HPD.95", "upper.HPD.95","post.mode","pMCMC")

#reduced model
polysum <- as.data.frame(matrix(NA, nrow=9, ncol=5))
rownames(polysum) <- c("Intercept","PreClim","InvClim","InvTime","Spp","units","s2","CVspp","CVresid")
colnames(polysum) <- c("post.mean", "lower.HPD.95", "upper.HPD.95","post.mode","pMCMC")

#write functions for h2phylo, s2, and coefficients of variation - compute generation by generation
CVfun <- function(var, mean) {
  sqrt(var)/mean
}

s2fun <- function(spp, resid){
  spp/(spp + resid)
}

h2fun <- function(phylo, spp, resid){
  phylo/(phylo + spp + resid)
}

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

dicscore <- matrix(NA, nrow=postsize, ncol=4)
colnames(dicscore) <- c("Iteration","PhyloDIC","PolyDIC", "Diff")

hpd <- matrix(NA, nrow=postsize, ncol=15)
colnames(hpd) <- c("Iteration", "lowerInt", "upperInt", "lPreClim", "uPreClim", "lInvClim", "uInvClim","lInvTime","uInvTime",
						"lSpecies", "uSpecies", "lAbSpec", "uAbSpec", "lUnits", "uUnits")
polyhpd <- matrix(NA, nrow=postsize, ncol=13)
colnames(polyhpd) <- c("Iteration", "lowerInt", "upperInt", "lPreClim", "uPreClim", "lInvClim", "uInvClim","lInvTime","uInvTime",
						"lAbSpec", "uAbSpec", "lUnits", "uUnits")

sampcor <- matrix(NA, nrow=100, ncol=1)

print("Output Dataframes Constructed")


#loop management
i <- 1
j <- 1 #use i to index run in model 1, use j to index output files
rlow <- 1 #use g to index net posterior positions
rhigh <- 2000

##--NOW GO!
while(j<postsize+1){
  print(j)
  
  ## Pull slopes from model iteration
  tmpMod <- bincmodel1
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
      slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciescristatellus"]
    } else { #if the species is a green anole (range 11-20)
      if(CID>10 | CID<21){ #sagrei
        slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciessagrei"]
      } else {
        if(CID>20 | CID<31){ #carolinensis
          slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"]
        } else {
          if(CID>30 | CID<43){ #chlorocyanus
            slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Specieschlorocyanus"]
          } else {
            if(CID>42 | CID<53){ #distichus
              slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciesdistichus"]
            } else {
              if(CID>56| CID<68){ #cybotes
                slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciescybotes"]
              } else {
                if(CID>67| CID<78){ #equestris
                  slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciesequestris"]
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
  bincmodel2 <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                        random=~Species + AbSpec,
                        ginverse=list(Species=Ainv), 
                        data=iframe,
                        prior=pr2A,
                        family = rep("gaussian", 1),
                        nitt=550000, burnin=50000, thin=250, verbose=FALSE)
  
  
  
  ## Run Model 2 without phylogeny
  set.seed(330)
  bincmodel2poly <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                            random=~AbSpec,
                            #ginverse=list(Species=Ainv), 
                            data=iframe,
                            prior=pr2B,
                            family = rep("gaussian", 1),
                            nitt=550000, burnin=50000, thin=250, verbose=FALSE)
  
  ## save posteriors
  tmpbig <- as.data.frame(cbind(bincmodel2$Sol, bincmodel2$VCV))
  tmpbigpoly <- as.data.frame(cbind(bincmodel2poly$Sol,bincmodel2poly$VCV))  

  #calculate h2, s2, and coefficients of variation
  tmpbig$h2phylo <- h2fun(tmpbig$Species, tmpbig$AbSpec, tmpbig$units)
  tmpbig$s2 <- h2fun(tmpbig$AbSpec, tmpbig$Species, tmpbig$units)
  tmpbigpoly$s2 <- s2fun(tmpbig$AbSpec, tmpbig$units)
  
  tmpbig$CVphylo <- CVfun(tmpbig$Species, tmpbig$`(Intercept)`)
  tmpbig$CVspp <- CVfun(tmpbig$AbSpec, tmpbig$`(Intercept)`)
  tmpbig$CVresid <- CVfun(tmpbig$units, tmpbig$`(Intercept)`)
  tmpbigpoly$CVspp <- CVfun(tmpbig$AbSpec, tmpbig$`(Intercept)`)
  tmpbigpoly$CVresid <- CVfun(tmpbig$units, tmpbig$`(Intercept)`)
  
  #add summary stats to growing net posterior distribution
  incmod2.mcmc[rlow:rhigh,] <- as.data.frame(tmpbig)
  incmod2poly.mcmc[rlow:rhigh,] <- as.data.frame(tmpbigpoly)

  #grab relevant output
  tmpsum <- summary(bincmodel2)
  tmpheifx <- heidel.diag(bincmodel2$Sol)
  tmpheiran <- heidel.diag(bincmodel2$VCV)
  tmpsumpoly <- summary(bincmodel2poly)
  tmpheifxpoly <- heidel.diag(bincmodel2poly$Sol)
  tmpheiranpoly <- heidel.diag(bincmodel2poly$VCV)
  tmpacrand <- autocorr.diag(bincmodel2$VCV)[2,]
  tmpacrandpoly <- autocorr.diag(bincmodel2poly$VCV)[2,]
  tmpacfix <- autocorr.diag(bincmodel2$Sol)[2,]
  tmpacfixpoly <- autocorr.diag(bincmodel2poly$Sol)[2,]
  tmphpd <- rbind(HPDinterval(bincmodel2$Sol), HPDinterval(bincmodel2$VCV))
  tmphpdpoly <- rbind(HPDinterval(bincmodel2poly$Sol), HPDinterval(bincmodel2poly$VCV))
  
  #random effect sampling correlation for full model
  sampcor[j,1] <- cor(bincmodel2$VCV[,1], bincmodel2$VCV[,2])
  
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
  rlow <- rlow + 2000
  rhigh <- rhigh + 2000
}

##-- Calculate Summary Statistics

write.csv(incmod2.mcmc, file="results/run1/inc_fullmod_mcmc_out.csv")
write.csv(incmod2poly.mcmc, file="results/run1/inc_reducemod_mcmc_out.csv")

#perform mean calculation while still data frames
for(n in rownames(m2sum)){
  m2sum[n,"post.mean"] <- mean(incmod2.mcmc[,n])
}

for(n in rownames(polysum)){
  polysum[n,"post.mean"] <- mean(incmod2poly.mcmc[,n])
}

#convert to mcmc objects to calculate HPD intervals and posterior modes
incmod2.mcmc <- as.mcmc(incmod2.mcmc)
incmod2poly.mcmc <- as.mcmc(incmod2poly.mcmc)

#calculate HPD and modes
for(n in rownames(m2sum)){
  m2sum[n,"lower.HPD.95"] <- HPDinterval(incmod2.mcmc[,n])[1]
  m2sum[n,"upper.HPD.95"] <- HPDinterval(incmod2.mcmc[,n])[2]
  m2sum[n,"post.mode"] <- posterior.mode(incmod2.mcmc[,n])
}

for(n in rownames(polysum)){
  polysum[n,"lower.HPD.95"] <- HPDinterval(incmod2poly.mcmc[,n])[1]
  polysum[n,"upper.HPD.95"] <- HPDinterval(incmod2poly.mcmc[,n])[2]
  polysum[n,"post.mode"] <- posterior.mode(incmod2poly.mcmc[,n])
}

#calculate pMCMC
#calculate pMCMC
for (n in rownames(m2sum)){
  m2sum[n,"pMCMC"] <- 2 * pmax(0.5/dim(incmod2.mcmc)[1], pmin(colSums(incmod2.mcmc[,n, drop = FALSE] > 0)/dim(incmod2.mcmc)[1],
 1 - colSums(incmod2.mcmc[, n, drop = FALSE] > 0)/dim(incmod2.mcmc)[1]))
}

for (n in rownames(polysum)){
  polysum[n,"pMCMC"] <- 2 * pmax(0.5/dim(incmod2poly.mcmc)[1], pmin(colSums(incmod2poly.mcmc[,n, drop = FALSE] > 0)/dim(incmod2poly.mcmc)[1],
 1 - colSums(incmod2poly.mcmc[, n, drop = FALSE] > 0)/dim(incmod2poly.mcmc)[1]))
}

##--export results
#Model 1 Output
incmodel1out <- summary(bincmodel1)
sink("results/run1/inc_summ_model1.txt")
print(incmodel1out)
sink()
incmode1 <- posterior.mode(bincmodel1$VCV)
write.csv(incmode1, file="results/run1/inc_mode_1.csv")
incautocor1 <- autocorr.diag(bincmodel1$VCV)
write.csv(incautocor1, file="results/run1/inc_autocorr_VCV1.csv")
incheidelrout1 <- heidel.diag(bincmodel1$VCV)
write.csv(incheidelrout1, file="results/run1/inc_heidelrout_1.csv")
incheidelfout1 <- heidel.diag(bincmodel1$Sol)
write.csv(incheidelfout1, file="results/run1/inc_heidelfout_1.csv")
write.csv(rbind(HPDinterval(bincmodel1$Sol), HPDinterval(bincmodel1$VCV)), file="results/run1/inc_hpds_model1")

#Model 2 Output
write.csv(m2sum, file="results/run1/inc_fullmod_summary.csv")
write.csv(polysum, file="results/run1/inc_reducemod_summary.csv")

write.csv(dicscore, file="results/run1inc_dics.csv")

write.csv(heidelfout, file="results/run1/inc_heidelfout_longrun.csv")
write.csv(heidelrout, file="results/run1/inc_heidelrout_longrun.csv")
write.csv(autofix, file="results/run1/inc_autocorr_fix_longrun.csv")
write.csv(autorand, file="results/run1/inc_autocorr_rand_longrun.csv")

write.csv(polyheidelfout, file="results/run1/poly_inc_heidelfout_longrun.csv")
write.csv(polyheidelrout, file="results/run1/poly_inc_heidelrout_longrun.csv")
write.csv(polyautofix, file="results/run1/poly_inc_autocorr_fix_longrun.csv")
write.csv(polyautorand, file="results/run1/poly_inc_autocorr_rand_longrun.csv")

write.csv(sampcor, file="results/run1/rf_sampling_correlations_longrun.csv")
