

## The code assumes that the working directory is "drugSNF" folder
rm(list=ls())

source('./RCode/preprocessInput.R')
source('./RCode/sensitivityData.R')
source('./RCode/perturbationData.R')
source('./RCode/structureData.R')
source('./RCode/constStructureLayer.R')
source('./RCode/constSensitivityLayer.R')
source('./RCode/constPerturbationLayer.R')
source('./RCode/integrateStrctSensPert.R')
source('./RCode/drugTargetBench.R')
source('./RCode/ATCBench.R')
source('./RCode/generateDrugPairs.R')
source('./RCode/compConcordIndx.R')
source('./RCode/generateRocPlot.R')
source('./RCode/generatePRPlot.R')
source('./RCode/predPerf.R')
source('./RCode/communityGen.R')


library(PharmacoGx) 
library(apcluster)
library(rcdk)
library(fingerprint)
library(annotate)
library(org.Hs.eg.db)
library(SNFtool)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
if ( ! file.exists("Output")) {
  dir.create("Output")
}

cDrugs <- preprocessInput(dname="nci60", "lincs")
dim(cDrugs$lincsboth)  ##  238 X 28
dim(cDrugs$nciboth)  ##  238 X 66


## loading and cleaning data for the layers 
sensData <- sensitivityData("nci60", cDrugs$nciboth)
dim(sensData)  ##[1]  60 238
pertData <- perturbationData("lincs", cDrugs$lincsboth, "nci60")
dim(pertData)  ##[1]  978 238
strcData <- structureData("lincs", cDrugs$lincsboth)  ## a vector
length(strcData) ## 238

#Reduce All Matrices to lowest common set of drugs across all 3
# Get 237 drugs now in the reduced sets
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                     sort(colnames(pertData))))

##Reduce the data to the set of common drugs across the three layers
strcData<- strcData[commonDrugs]
sensData <- sensData[,commonDrugs]
pertData <- pertData[,commonDrugs]
## filter out cDrugs$lincsboth "dataframe" as well accordingly
cDrugs$lincsboth <- cDrugs$lincsboth[cDrugs$lincsboth$pert_iname %in% commonDrugs, ]


#Sanity Checks
if (ncol(sensData) != ncol(pertData)) stop(sprintf("error!"))
if (ncol(sensData) !=length(strcData)) stop(sprintf("error!"))
if (all(colnames(pertData) != colnames(sensData))) stop(sprintf("error!"))
if (all(colnames(pertData) != names(strcData))) stop(sprintf("error!"))


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)
integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/nci60-Integrated.RData")

# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)


## Benchmarking and validation
## 1- DRUG-TARGET 
dataBench1.5 <- drugTargetBench("uniprot", commonDrugs)
dim(dataBench1.5) ##[1] 86
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)

pairs1 <- generateDrugPairs(dataBench1.5, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL)
## compare cindices of integration layer vs. all single layers (e.g., structure, perturbation, sensitivity )
res1 <- compConcordIndx(pairs1)
cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res1$cindxLst$integrCindex$c.index, "\n structure: ", res1$cindxLst$structureLayerCindex$c.index,
    "\n perturbation: ",  res1$cindxLst$perturbationLayerCindex$c.index, "\n sensitivity: ", res1$cindxLst$sensitivityLayerCindex$c.index)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res1$pVals$intgrStrcPVal,"\n perturbation: ", res1$pVals$intgrPertPVal,
    "\n sensitivity: ", res1$pVals$intgrSensPVal, "\n Iorio: ", res1$pVals$intgrIorioPVal, "\n Iskar: ", res1$pVals$intgrIskarPVal)
## ROC and PR plots
generateRocPlot(pairs1, d1Name="nci60", d2Name="lincs", benchNam="drug-target(UNIPROT)-sep16")
generatePRPlot(pairs1, d1Name="nci60", d2Name="lincs", benchNam="drug-target(UNIPROT)-sep16")


## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs$lincsboth)
dim(dataBench3) ##[1] 72 72
##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/SuperPredsimil-NCI60.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
##load "iskar" results here ... (see iskar.R)
load("Data/averageIskarFinal.RData")

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, SuperPredsimil)

## compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2)
cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res2$cindxLst$integrCindex$c.index, "\n structure: ", res2$cindxLst$structureLayerCindex$c.index,
    "\n perturbation: ",  res2$cindxLst$perturbationLayerCindex$c.index, "\n sensitivity: ", res2$cindxLst$sensitivityLayerCindex$c.index)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res2$pVals$intgrStrcPVal,"\n perturbation: ", res2$pVals$intgrPertPVal,
    "\n sensitivity: ", res2$pVals$intgrSensPVal, "\n Iorio: ", res2$pVals$intgrIorioPVal, "\n Iskar: ", res2$pVals$intgrIskarPVal)
## ROC and PR plots
generateRocPlot(pairs2, d1Name="nci60", d2Name="lincs", benchNam="ATC(CHEMBL)-Zscore-sep16")
generatePRPlot(pairs2, d1Name="nci60", d2Name="lincs", benchNam="ATC(CHEMBL)-Zscore-sep16")


## generate communities
load("Output/gmt_targ_chembl.RData")
communityGen(integrtStrctSensPert, "nci60", GMT_TARG)







