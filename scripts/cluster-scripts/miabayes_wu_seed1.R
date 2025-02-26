setwd("/Users/morganmuell/Documents/PhD/Miami_Project/data-analysis/github/MiamiAnoles_ReactionNorms/datasets/phenotype_data")

##--Load libraries
library(MCMCglmm)
library(nadiv)
library(phytools)
library(gremlin)

##--Prep Data and Priors
load("wutrim.Rdata")
load("model2df.rda")

tree <- read.tree(file="../../pruned_anole_tree.treefile")

Ainv <- inverseA(tree)$Ainv

prior1 <- list(R = list(V = 1e-12, nu = -2), #<-- Non-informative improper: *marginal* posterior equal to REML estimate
               G = list(G1 = list(V = diag(2)*0.002, nu = 3, 
                                  alpha.mu = c(0,0), alpha.V = diag(2)*10000)))
#pr2A <- list(R = list(V = 1, nu = 0.002), 
#             G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10), 
#                      G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10))) 

prior2 <- list(R = list(V = 1, nu = 0.002), 
             G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*10)))

##-- Run Model 1
set.seed(330)
bwumodel1 <- MCMCglmm(logWU ~ 1 + EggMass + nTreatment + PercDev + Species
                      + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = wu, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = FALSE, prior = prior1)

print("Model 1 Complete")

##--Housekeeping prior to looop
postsize <- 100

# Output mcmc object for all model names
wumod2.mcmc <- as.data.frame(matrix(NA, nrow=200000, ncol=12))
colnames(wumod2.mcmc) <- c("Intercept","PreClim","InvClim","InvTime","Phylo","Spp","units","h2phylo","s2","CVphylo","CVspp","CVresid")
wumod2poly.mcmc <- as.data.frame(matrix(NA, nrow=200000, ncol=9))
colnames(wumod2poly.mcmc) <- c("Intercept","PreClim","InvClim","InvTime","Spp","units","s2","CVspp","CVresid")

# Create new model summary frame
#full model
m2sum <- as.data.frame(matrix(NA, nrow=12, ncol=5))
rownames(m2sum) <- c("Intercept","PreClim","InvClim","InvTime","Phylo","Spp","units","h2phylo","s2","CVphylo","CVspp","CVresid")
colnames(m2sum) <- c("post.mean", "lower.HPD.95", "upper.HPD.95","post.mode","pMCMC")

#reduced model
polysum <- as.data.frame(matrix(NA, nrow=9, ncol=5))
rownames(polysum) <- c("Intercept","PreClim","InvClim","InvTime","Spp","units","s2","CVspp","CVresid")
colnames(polysum) <- c("post.mean", "lower.HPD.95", "upper.HPD.95","post.mode","pMCMC")

#gremlin summary
gremsum <- as.data.frame(matrix(NA, nrow=9, ncol=5))
rownames(gremsum) <- c("Tree", "SpeciesCat", "Residual", "Intercept", "scPreClim", "scInvClim","InvTime","h2phylo","s2")
colnames(gremsum) <- c("post.mean", "lower.HPD.95", "upper.HPD.95","post.mode","pMCMC")

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
heidelrout <- matrix(NA, nrow=postsize, ncol=9) #p-values for stationarity, random effect: iteration, spppass, sppp, unitspass, unitsp, spphalfpass, unitshalfpass
colnames(heidelrout) <- c("Iteration", "Pass.Stat.Spp", "Pvalue.Spp", "Pass.Stat.Resid", 
                          "Pvalue.Resid", "Pass.Halfwidth.Spp", "Ratio.Spp", "Pass.Halfwidth.Resid",
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
autorand <- matrix(NA, nrow=postsize, ncol=3) #mean autocorr for random effects
colnames(autorand) <- c("Iteration", "Species", "units")

polyautofix <- matrix(NA, nrow=postsize, ncol=5) #mean autocorr for fixed effects
colnames(polyautofix) <- c("Iteration", "Intercept", "scPreClim", "scInvClim","InvTime")
polyautorand <- matrix(NA, nrow=postsize, ncol=3) #mean autocorr for random effects
colnames(polyautorand) <- c("Iteration", "AbSpec", "units")

dicscore <- matrix(NA, nrow=postsize, ncol=4)
colnames(dicscore) <- c("Iteration","PhyloDIC","PolyDIC", "Diff")

gremfout <- as.data.frame(matrix(NA, nrow=100, ncol=14))
colnames(gremfout) <- c("Intercept","Intercept.StdError","Intercept.z.value", 
                        "scPreClim", "scPreClim.StdError", "scPreClim.z.value",
                        "scInvClim","scInvClim.StdError","scInvClim.z.value",
                        "InvTime","InvTime.StdError","InvTime.z.value","h2phylo","s2")

gremrout <- as.data.frame(matrix(NA, nrow=100, ncol=6))
colnames(gremrout) <- c("Tree","Tree.StdError", "SpeciesCat", 
                        "SpeciesCat.StdError","Residual","Residual.StdError")

sampcor <- matrix(NA, nrow=100, ncol=1)
gremcor <- matrix(NA, nrow=100, ncol=1) #correlations btw tree and species random effects

print("Output Dataframes Constructed")


#loop management
i <- 1 #use i to index run in model 1
j <- 1 # use j to index output files
rlow <- 1 #use g to index net posterior positions
rhigh <- 2000

##--NOW GO!
while(j<postsize+1){
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
    CID <- as.integer(slopeframe[cg,"cageNumbs"]) #make variable for cage ID
    if(CID<11){ #cristatellus
      slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciescristatellus"]
      #print(paste0(CID,": cristatellus picked"))
    } else { #if the species is a green anole (range 11-20)
      if(CID>10 & CID<21){ #sagrei
        slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciessagrei"]
        #print(paste0(CID,": sagrei picked"))
      } else {
        if(CID>20 & CID<31){ #carolinensis
          slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"]
          #print(paste0(CID,": carolinensis picked"))
        } else {
          if(CID>30 & CID<43){ #chlorocyanus
            slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Specieschlorocyanus"]
            #print(paste0(CID,": chlorocyanus picked"))
          } else {
            if(CID>42 & CID<53){ #distichus
              slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciesdistichus"]
              #print(paste0(CID,": distichus picked"))
            } else {
              if(CID>56 & CID<68){ #cybotes
                slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciescybotes"]
                #print(paste0(CID,": cybotes picked"))
              } else {
                if(CID>67 & CID<78){ #equestris
                  slopeframe[cg,"islps"] <- slopeframe[cg,"devs"] + tmpMod$Sol[i,"nTreatment"] + tmpMod$Sol[i,"nTreatment:Speciesequestris"]
                  #print(paste0(CID,": equestris picked"))
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
  iframe$scislps <- scale(iframe$islps)
  iframe$sislps <- scale(iframe$islps, center=FALSE)
  
  ## Run Model 2 with phylogeny 
  set.seed(330)
  bwumodel2 <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                        random=~Species,
                        ginverse=list(Species=Ainv), 
                        data=iframe,
                        prior=prior2,
                        family = rep("gaussian", 1),
                        nitt=1500000, burnin=500000, thin=500, verbose=FALSE) #drop to 5000 burn-in
  summary(bwumodel2)
  autocorr.diag(bwumodel2$VCV)[2,]
  
  gremwu <- gremlin(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
          random=~Species + AbSpec,
          ginverse=list(Species=Ainv), 
          data=iframe,
          v=0)
  
  #approximate standard errors for each variance component - sometimes differences in scale if it looks weired
  #inverse info matrix shows sampling variance and covariance among variance params
  #then convert covariation to correlation matrix, gets rid of sampling variance in each variance component, converts to correlation
  
  #they are very highly negatively correlated, meaning when model tries variance estimate that's high for first term,
  # way to increase likelihood is to decrease the second term's variance
  
  #posteriors are heavily influenced by prior in MCMCglmm, so sampling correlation isn't too high
  # between variance components bc sampling diff distributions, but w/in variance commponents getting
  # high autocorrelation because just sampling prior density near zero
  
  #other problem is estimating variance from 7 datapoints and correlated 7 points in tree
  
  
  
  #tree variance - variance among spp due to shared evo history fits Brownian Motion
  # if estimated with precision such that residuals represent within-species variance plus any among-spp variance
  # not due to correlation in BM model and can't be explained by climate, then q is
  # what are among-spp variance not captured by hist clim, current clim, invtime
  
  
  ## Run Model 2 without phylogeny
  set.seed(330)
  bwumodel2poly <- MCMCglmm(sislps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                            random=~AbSpec,
                            #ginverse=list(Species=Ainv), 
                            data=iframe,
                            prior=prior2,
                            family = rep("gaussian", 1),
                            nitt=1500000, burnin=500000, thin=500, verbose=FALSE)
  autocorr.diag(bwumodel2poly$VCV)[2,]

  ## save posteriors
  tmpbig <- as.data.frame(cbind(bwumodel2$Sol, bwumodel2$VCV))
  tmpbigpoly <- as.data.frame(cbind(bwumodel2poly$Sol,bwumodel2poly$VCV))
  
  #calculate h2, s2, and coefficients of variation
  tmpbig$h2phylo <- s2fun(tmpbig$Species, tmpbig$units)
  tmpbigpoly$s2 <- s2fun(tmpbigpoly$AbSpec, tmpbigpoly$units) #this was a mistake, wasn't pulling from poly mcmc
  
  tmpbig$CVphylo <- CVfun(tmpbig$Species, tmpbig$`(Intercept)`)
  #tmpbig$CVspp <- CVfun(tmpbig$AbSpec, tmpbig$`(Intercept)`)
  tmpbig$CVresid <- CVfun(tmpbig$units, tmpbig$`(Intercept)`)
  tmpbigpoly$CVspp <- CVfun(tmpbigpoly$AbSpec, tmpbigpoly$`(Intercept)`)
  tmpbigpoly$CVresid <- CVfun(tmpbigpoly$units, tmpbigpoly$`(Intercept)`)

  
  #add summary stats to growing net posterior distribution
  wumod2.mcmc[rlow:rhigh,] <- as.data.frame(tmpbig)
  wumod2poly.mcmc[rlow:rhigh,] <- as.data.frame(tmpbigpoly)

  #grab relevant output
  tmpsum <- summary(bwumodel2)
  tmpsumpoly <- summary(bwumodel2poly)
  tmpgrem <- summary(gremwu)
  tmpheifx <- heidel.diag(bwumodel2$Sol)
  tmpheiran <- heidel.diag(bwumodel2$VCV)
  tmpheifxpoly <- heidel.diag(bwumodel2poly$Sol)
  tmpheiranpoly <- heidel.diag(bwumodel2poly$VCV)
  tmpacrand <- autocorr.diag(bwumodel2$VCV)[2,]
  tmpacrandpoly <- autocorr.diag(bwumodel2poly$VCV)[2,]
  tmpacfix <- autocorr.diag(bwumodel2$Sol)[2,]
  tmpacfixpoly <- autocorr.diag(bwumodel2poly$Sol)[2,]
  
  
  ##--pull effects and convergence diagnostics and store 
  #random effect sampling correlation for full model
  sampcor[j,1] <- cor(bwumodel2$VCV[,1], bwumodel2$VCV[,2])
  gremcor[j,1] <- cov2cor(solve(gremwu$grMod$AI))[2] #correlation btw tree and spp
  
  ##--gremlin initial summary statistics
  gremfout[j,"Intercept.Estimate"] <- tmpgrem$coefficients[1,1]
  gremfout[j,"Intercept.StdError"] <- tmpgrem$coefficients[1,2]
  gremfout[j,"Intercept.z.value"] <- tmpgrem$coefficients[1,3]
  gremfout[j,"scPreClim.Estimate"] <- tmpgrem$coefficients[2,1]
  gremfout[j,"scPreClim.StdError"] <- tmpgrem$coefficients[2,2]
  gremfout[j,"scPreClim.z.value"] <- tmpgrem$coefficients[2,3]
  gremfout[j,"scInvClim.Estimate"] <- tmpgrem$coefficients[3,1]
  gremfout[j,"scInvClim.StdError"] <- tmpgrem$coefficients[3,2]
  gremfout[j,"scInvClim.z.value"] <- tmpgrem$coefficients[3,3]
  gremfout[j,"InvTime.Estimate"] <- tmpgrem$coefficients[4,1]
  gremfout[j,"InvTime.StdError"] <- tmpgrem$coefficients[4,2]
  gremfout[j,"InvTime.z.value"] <- tmpgrem$coefficients[4,3]
  
  gremrout[j,"Tree.Estimate"] <- tmpgrem$varcompSummary["G.Species","Estimate"]
  gremrout[j,"SpeciesCat.Estimate"] <- tmpgrem$varcompSummary["G.AbSpec","Estimate"]
  gremrout[j,"Residual.Estimate"] <- tmpgrem$varcompSummary["ResVar1","Estimate"]
  gremrout[j,"Tree.StdError"] <- tmpgrem$varcompSummary["G.Species","Std. Error"]
  gremrout[j,"SpeciesCat.StdError"] <- tmpgrem$varcompSummary["G.AbSpec","Std. Error"]
  gremrout[j,"Residual.StdError"] <- tmpgrem$varcompSummary["ResVar1","Std. Error"]
  gremrout[j,"h2phylo"] <- h2fun(tmpgrem$varcompSummary["G.Species","Estimate"],tmpgrem$varcompSummary["G.AbSpec","Estimate"], tmpgrem$varcompSummary["ResVar1","Estimate"])
  gremrout[j,"s2"] <- s2fun(tmpgrem$varcompSummary["G.AbSpec","Estimate"],tmpgrem$varcompSummary["ResVar1","Estimate"])
  
  #DIC scores
  dicscore[j,1] <- i
  dicscore[j,2] <- tmpsum$DIC
  dicscore[j,3] <- tmpsumpoly$DIC
  dicscore[j,4] <- tmpsum$DIC - tmpsumpoly$DIC
  
  #autocorrelation
  autofix[j,1] <- i
  autofix[j,c(2:5)] <- tmpacfix[c(1:4)]
  autorand[j,1] <- i
  autorand[j,c(2:3)] <- tmpacrand[c(1:2)] #fix
  
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
  heidelrout[j,4] <- tmpheiran[2,1] #residuals pass/fail, stationarity
  heidelrout[j,5] <- tmpheiran[2,3] #residuals p-value, stationarity
  heidelrout[j,6] <- tmpheiran[1,4] #species pass/fail, halfwidth
  heidelrout[j,7] <- (tmpheiran[1,6] / tmpheiran[1,5]) #species ratio
  heidelrout[j,8] <- tmpheiran[2,4] #residuals pass/fail, halfwidth
  heidelrout[j,9] <- (tmpheiran[2,6] / tmpheiran[2,5]) #residuals ratio
  
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

write.csv(wumod2.mcmc, file="results/run1/fullmod_mcmc_out.csv")
write.csv(wumod2poly.mcmc, file="results/run1/reducemod_mcmc_out.csv")
write.csv(gremfout, file="results/run1/gremlin_fixed_out.csv")
write.csv(gremrout, file="results/run1/gremlin_rand_out.csv")

#perform mean calculation while still data frames
for(n in rownames(m2sum)){
  m2sum[n,"post.mean"] <- mean(wumod2.mcmc[,n])
}

for(n in rownames(polysum)){
  polysum[n,"post.mean"] <- mean(wumod2poly.mcmc[,n])
}

for(n in rownames(gremsum)){
  gremsum[n,"post.mean"] <- mean(cbind(gremfout, gremrout)[,n])
}

#convert to mcmc objects to calculate HPD intervals and posterior modes
wumod2.mcmc <- as.mcmc(wumod2.mcmc)
wumod2poly.mcmc <- as.mcmc(wumod2poly.mcmc)
wugrem.mcmc <- as.mcmc(cbind(gremfout, gremrout))

#calculate HPD and modes
for(n in rownames(m2sum)){
  m2sum[n,"lower.HPD.95"] <- HPDinterval(wumod2.mcmc[,n])[1]
  m2sum[n,"upper.HPD.95"] <- HPDinterval(wumod2.mcmc[,n])[2]
  m2sum[n,"post.mode"] <- posterior.mode(wumod2.mcmc[,n])
}

for(n in rownames(gremsum)){
  polysum[n,"lower.HPD.95"] <- HPDinterval(wumod2poly.mcmc[,n])[1]
  polysum[n,"upper.HPD.95"] <- HPDinterval(wumod2poly.mcmc[,n])[2]
  polysum[n,"post.mode"] <- posterior.mode(wumod2poly.mcmc[,n])
}

for(n in rownames(gremsum)){
  gremsum[n,"lower.HPD.95"] <- HPDinterval(wugrem.mcmc[,n])[1]
  gremsum[n,"upper.HPD.95"] <- HPDinterval(wugrem.mcmc[,n])[2]
  gremsum[n,"post.mode"] <- posterior.mode(wugrem.mcmc[,n])
}

#calculate pMCMC
for (n in rownames(m2sum)){
  m2sum[n,"pMCMC"] <- 2 * pmax(0.5/dim(wumod2.mcmc)[1], pmin(colSums(wumod2.mcmc[,n, drop = FALSE] > 0)/dim(wumod2.mcmc)[1],
 1 - colSums(wumod2.mcmc[, n, drop = FALSE] > 0)/dim(wumod2.mcmc)[1]))
}

for (n in rownames(polysum)){
  polysum[n,"pMCMC"] <- 2 * pmax(0.5/dim(wumod2poly.mcmc)[1], pmin(colSums(wumod2poly.mcmc[,n, drop = FALSE] > 0)/dim(wumod2poly.mcmc)[1],
 1 - colSums(wumod2poly.mcmc[, n, drop = FALSE] > 0)/dim(wumod2poly.mcmc)[1]))
}

for (n in rownames(gremsum)){
  gremsum[n,"pMCMC"] <- 2 * pmax(0.5/dim(wugrem.mcmc)[1], pmin(colSums(wugrem.mcmc[,n, drop = FALSE] > 0)/dim(wugrem.mcmc)[1],
                                                                   1 - colSums(wugrem.mcmc[, n, drop = FALSE] > 0)/dim(wugrem.mcmc)[1]))
}


##--export results
#Model 1 Output
wumodel1out <- summary(bwumodel1)
sink("results/run1/wu_summ_model1.txt")
print(wumodel1out)
sink()
wumode1 <- posterior.mode(bwumodel1$VCV)
write.csv(wumode1, file="results/run1/wu_mode_1.csv")
wuautocor1 <- autocorr.diag(bwumodel1$VCV)
write.csv(wuautocor1, file="results/run1/wu_autocorr_VCV1.csv")
wuheidelrout1 <- heidel.diag(bwumodel1$VCV)
write.csv(wuheidelrout1, file="results/run1/wu_heidelrout_1.csv")
wuheidelfout1 <- heidel.diag(bwumodel1$Sol)
write.csv(wuheidelfout1, file="results/run1/wu_heidelfout_1.csv")
write.csv(rbind(HPDinterval(bwumodel1$Sol), HPDinterval(bwumodel1$VCV)), file="results/run1/wu_hpds_model1")

#Model 2 Output
write.csv(m2sum, file="results/run1/fullmod_summary.csv")
write.csv(polysum, file="results/run1/reducemod_summary.csv")
write.csv(gremsum, file="results/run1/gremlinmod_summary.csv")

write.csv(dicscore, file="results/run1/wu_dics.csv")

write.csv(heidelfout, file="results/run1/wu_heidelfout_longrun.csv")
write.csv(heidelrout, file="results/run1/wu_heidelrout_longrun.csv")
write.csv(autofix, file="results/run1/wu_autocorr_fix_longrun.csv")
write.csv(autorand, file="results/run1/wu_autocorr_rand_longrun.csv")


write.csv(polyheidelfout, file="results/run1/poly_wu_heidelfout_longrun.csv")
write.csv(polyheidelrout, file="results/run1/poly_wu_heidelrout_longrun.csv")
write.csv(polyautofix, file="results/run1/poly_wu_autocorr_fix_longrun.csv")
write.csv(polyautorand, file="results/run1/poly_wu_autocorr_rand_longrun.csv")

write.csv(sampcor, file="results/run1/rf_sampling_correlations_longrun.csv")
write.csv(gremcor, file="results/run1/randeff_sampling_correlations_gremlin.csv")
