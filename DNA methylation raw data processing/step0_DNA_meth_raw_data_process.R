##### description: Processing of raw DNA methylation data from Illumina 850K platform.

##### 01Start: Load raw DNA methylation data #####

#load raw DNA methylation data of M90
library(ChAMP)
myload01=champ.load(directory='90-day spaceflight_M90/DNA methylation/rawdata',
 methValue = 'B',filterXY = FALSE,
 filterDetP = TRUE,detPcut = 0.01,
 filterBeads=TRUE,beadCutoff=0.05,
 filterNoCG=FALSE,filterSNPs=TRUE,
 filterMultiHit=TRUE,arraytype='EPIC') 
saveRDS(myload01, file="90-day spaceflight_M90/DNA methylation/beta/beta.Rds")

#load raw DNA methylation data of M180-1
myload02=champ.load(directory='180-day spaceflight_M180-1/DNA methylation/rawdata',
 methValue = 'B',filterXY = FALSE,
 filterDetP = TRUE,detPcut = 0.01,
 filterBeads=TRUE,beadCutoff=0.05,
 filterNoCG=FALSE,filterSNPs=TRUE,
 filterMultiHit=TRUE,arraytype='EPIC') 
saveRDS(myload02, file="180-day spaceflight_M180-1/DNA methylation/beta/beta.Rds")
##### 01End: load raw DNA methylation data #####


##### 02Start: Normalization of DNA methylation data #####

#Data normalization of M90
myload01 <- readRDS("90-day spaceflight_M90/DNA methylation/beta/beta.Rds")
myNorm01 <- champ.norm(beta=myload01$beta,
                     rgSet=myload01$rgSet,
                     mset=myload01$mset,
                     resultsDir="./CHAMP_Normalization/",
                     method="BMIQ",
                     plotBMIQ=TRUE,
                     arraytype="850K",
                     cores=3)
saveRDS(myNorm01, file="90-day spaceflight_M90/DNA methylation/beta/Norm_beta.Rds")

#Data normalization of M180-1
myload02 <- readRDS("180-day spaceflight_M180-1/DNA methylation/beta/beta.Rds")
myNorm02 <- champ.norm(beta=myload02$beta,
                       rgSet=myload02$rgSet,
                       mset=myload02$mset,
                       resultsDir="./CHAMP_Normalization/",
                       method="BMIQ",
                       plotBMIQ=TRUE,
                       arraytype="850K",
                       cores=3)
saveRDS(myNorm02, file="180-day spaceflight_M180-1/DNA methylation/beta/Norm_beta.Rds")
##### 02End: Normalization of DNA methylation data #####
