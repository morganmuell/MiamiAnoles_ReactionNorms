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
set.seed(254)
bincmodel1 <- MCMCglmm(logInc ~ 1 + EggMass + nTreatment + Species
                      + Species*nTreatment, 
                      random = ~us(1 + nTreatment):Cage, pr = TRUE,
                      data = incp, nitt = 275000, thin = 50, burnin = 25000,
                      verbose = FALSE, prior = prior1)
print("Model 1 Complete")

##--Housekeeping prior to looop
postsize <- 100

## Make df for ESS values
fullessvals <- matrix(NA, nrow=postsize, ncol=7)
colnames(fullessvals) <- c("Phylogeny","Species","Residuals","Intercept","InvClim","PreClim","InvTime")

redessvals <- matrix(NA, nrow=postsize, ncol=6)
colnames(redessvals) <- c("Species","Residuals","Intercept","InvClim","PreClim","InvTime")

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
  
  ## Run Model 2 with phylogeny 
  set.seed(254)
  bincmodel2 <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                        random=~Species + AbSpec,
                        ginverse=list(Species=Ainv), 
                        data=iframe,
                        prior=pr2A,
                        family = rep("gaussian", 1),
                        nitt=550000, burnin=50000, thin=250, verbose=FALSE)
  
  
  
  ## Run Model 2 without phylogeny
  set.seed(254)
  bincmodel2poly <- MCMCglmm(islps ~ scPreClim + scInvClim + at.level(dummyInv, "Invasive"):scInvTime,
                            random=~AbSpec,
                            #ginverse=list(Species=Ainv), 
                            data=iframe,
                            prior=pr2B,
                            family = rep("gaussian", 1),
                            nitt=550000, burnin=50000, thin=250, verbose=FALSE)


  #grab ESS objects
  fullsum <- summary(bincmodel2)
  polysum <- summary(bincmodel2poly)

  fullessvals[j,"Phylogeny"] <- fullsum$Gcovariances["Species","eff.samp"]
  fullessvals[j, "Species"] <- fullsum$Gcovariances["AbSpec", "eff.samp"]
  fullessvals[j, "Residuals"] <- fullsum$Rcovariances["units", "eff.samp"]
  fullessvals[j, "Intercept"] <- fullsum$solutions["(Intercept)", "eff.samp"]
  fullessvals[j, "InvClim"] <- fullsum$solutions["scInvClim", "eff.samp"]
  fullessvals[j, "PreClim"] <- fullsum$solutions["scPreClim", "eff.samp"]
  fullessvals[j, "InvTime"] <- fullsum$solutions[4, "eff.samp"]

  redessvals[j, "Species"] <- polysum$Gcovariances["AbSpec", "eff.samp"]
  redessvals[j, "Residuals"] <- polysum$Rcovariances["units", "eff.samp"]
  redessvals[j, "Intercept"] <- polysum$solutions["(Intercept)", "eff.samp"]
  redessvals[j, "InvClim"] <- polysum$solutions["scInvClim", "eff.samp"]
  redessvals[j, "PreClim"] <- polysum$solutions["scPreClim", "eff.samp"]
  redessvals[j, "InvTime"] <- polysum$solutions[4, "eff.samp"]

  
  #last but not least
  i<- i+50
  j<- j+1
  rlow <- rlow + 2000
  rhigh <- rhigh + 2000
  
}



##--export results
write.csv(fullessvals, file="results/run2/fullmod_ess.csv")
write.csv(redessvals, file="results/run2/polymod_ess.csv")

