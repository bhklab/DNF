###############################################################################################################
## Function intersects the benchmark data and the layers of the network and then generates drug pairs for the 
## final (intersection) subset of drugs
##
## input: 
##     benchDat
##     benchDat
##     strcAff
##     sensAff
##     pertAff
##     integration  
##
## output: 
##     list of drug pairs from benchmark, structure, perturbation, sensitivtiy and the integration of all
##
## 
###############################################################################################################

generateDrugPairs <- function( benchDat, strcAff, sensAff, pertAff, integration) {
    
## intersection of the bechmark set and the network layers
intx <- intersect(colnames(benchDat), colnames(strcAff))
dataStr <- strcAff[intx,intx]
dataStr[upper.tri(dataStr, diag=TRUE)] <- NA
dataSens <- sensAff[intx,intx]
dataSens[upper.tri(dataSens, diag=TRUE)] <- NA
dataPert <- pertAff[intx,intx]
dataPert[upper.tri(dataPert, diag=TRUE)] <- NA
dataIntegrAll <- integration[intx,intx]
dataIntegrAll[upper.tri(dataIntegrAll, diag=TRUE)] <- NA
    

##create pairs of drugs from benchmark 1 if same drug set from GMT and 0 otherwise 
benchDat[upper.tri(benchDat, diag=TRUE)] <- NA
benchPairs <- melt(benchDat)
benchPairs <- na.omit(benchPairs)
colnames(benchPairs)[3] <- "bench"

## create pairs of drugs from structure
strcPairs <- melt(dataStr)
strcPairs <- na.omit(strcPairs)
colnames(strcPairs)[3] <- "obs.str"

## create pairs of drugs from sensitivty
sensPairs <- melt(dataSens)
sensPairs <- na.omit(sensPairs)
colnames(sensPairs)[3] <- "obs.sens"

## create pairs of drugs from perturbation
pertPairs <- melt(dataPert)
pertPairs <- na.omit(pertPairs)
colnames(pertPairs)[3] <- "obs.pert"

## create pairs of drugs from combination
integrPairs <- melt(dataIntegrAll)
integrPairs <- na.omit(integrPairs)
colnames(integrPairs)[3] <- "obs.combiall"


res <- list(benchPairs=benchPairs, strcPairs=strcPairs, sensPairs=sensPairs, pertPairs=pertPairs, integrPairs=integrPairs)


}