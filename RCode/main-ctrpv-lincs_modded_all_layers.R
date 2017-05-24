

## The code assumes that the working directory is "drugSNF" folder
rm(list=ls())

source('./RCode/preprocessInputModded.R')
source('./RCode/sensitivityDataModded.R')
source('./RCode/perturbationData.R')
source('./RCode/perturbationDataModded.R')
source('./RCode/structureData.R')
source('./RCode/constStructureLayer.R')
source('./RCode/constSensitivityLayer.R')
source('./RCode/constSensitivityLayerCombined.R')
source('./RCode/constPerturbationLayer.R')
source('./RCode/integrateStrctSensPert.R')
source('./RCode/integrateStrctSensPertModded.R')
source('./RCode/drugTargetBench.R')
source('./RCode/drugTargetBenchModded.R')
source('./RCode/generateDrugPairs.R')
source('./RCode/generateDrugPairsModded.R')
source('./RCode/compConcordIndx.R')
source('./RCode/compConcordIndxModded.R')
source('./RCode/generateRocPlot.R')
source('./RCode/generateRocPlotModded.R')
source('./RCode/generatePRPlot.R')
source('./RCode/predPerf.R')
source('./RCode/ATCbench.R')
source('./RCode/communityGen.R')
source('./RCode/communityGenModded.R')
source('./RCode/cindexComp2.R')

source("./RCode/luminexData.R")
source("./RCode/constLuminexLayer.R")

source("./RCode/imagingData.R")
source("./RCode/constImagingLayer.R")

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
library(PRROC)

library(apcluster)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
#outputDir <- paste(getwd(), "/outputDir", sep="")
if ( ! file.exists("Output")) {
    dir.create("Output")
}

# Find common drugs between combined sensitivity and LINCS dataset
cDrugs <- preprocessInputModded(dname="combined", "lincs", "Data/combined_sens_adjusted_diag_datasets_with.RData")
dim(cDrugs)  ##239 X 28 for just CTRPv2, 309 x 28 for combined data

# Process Sensitivity, Perturbation, and Structure layers for set of common drugs
sensData <- sensitivityDataModded("combined", cDrugs, "Data/combined_sens_adjusted_diag_datasets_with.RData")  ## 645 X 239
dim(sensData) # 309 x 309 for combined data
pertData <- perturbationDataModded("Data/L1000_compound_signatures.RData", cDrugs)  ## 978 X 239
dim(pertData) # 978 x 237 for 
strcData <- structureData("lincs", cDrugs)  ## a vector  --> 239 elemnts
length(strcData)

luminexData <- LuminexData(cDrugs, badchars)
dim(luminexData)
imagingData <- ImagingData(cDrugs, badchars)
dim(imagingData)

## Get the common drugs (239) among the 3 datasets/layers
#commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
#                                     sort(colnames(pertData))))
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                     sort(colnames(pertData))))
length(commonDrugs) ## 239 ..

strcData<- strcData[commonDrugs] # 239 drugs
sensData <- sensData[commonDrugs, commonDrugs] # 645 x 239 drugs
pertData<- pertData[, commonDrugs] #978 genes x  239 
luminexData <- luminexData[, commonDrugs]
imagingData <- imagingData[, commonDrugs]


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayerCombined(sensData)
pertAffMat <- constPerturbationLayer(pertData)
luminexAffMat <- constLuminexLayer(luminexData)
imagingAffMat <- constImagingLayer(imagingData)
integrtStrctSensPert <- integrateStrctSensPertModded(sensAffMat, strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/ctrpv2-Integrated.RData")


## Benchmarking and validation
## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141
dataBench <- drugTargetBenchModded("ctrpv",  colnames(sensData), "gmt_targ_combined_with_gdsc_adjusted_diag.RData") # 141 x 141 drug-drug adjacency matrix --> 141
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Data/drugERankSimil-CTRPV2-kendall.Rdata")

#pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2)
pairs <- generateDrugPairsModded(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)


## compare cindices of combiantion layer vs. a single layer (e.g., structure)
#res <- compConcordIndx(pairs)
res <- compConcordIndxModded(pairs)

#cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res$cindxLst$integrCindex, "\n structure: ", res$cindxLst$structureLayerCindex,
#    "\n perturbation: ",  res$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res$cindxLst$iorioCindex, 
#    "\n Iskar: ", res$cindxLst$iskarCindex)
#cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
#    "\n sensitivity: ", res$pVals$intgrSensPVal, "\n Iorio: ", res$pVals$intgrIorioPVal, "\n Iskar: ", res$pVals$intgrIskarPVal)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res$cindxLst$integrCindex, "\n structure: ", res$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res$cindxLst$sensitivityLayerCindex)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
    "\n sensitivity: ", res$pVals$intgrSensPVal)

## ROC and PR plots
generateRocPlotModded(pairs, d1Name="ctrpv2", d2Name="lincs", benchName="drug-target_combined_adjusted_diagonal_with")
generatePRPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target_combined_adjusted_diagonal_with")

## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs)
dim(dataBench3) ##[1]  51 51
##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/superPredSimil-ctrp-kendall.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, superPredSimil2, NULL)

## compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res2$cindxLst$integrCindex, "\n structure: ", res2$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res2$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res2$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res2$cindxLst$iorioCindex, 
    "\n Iskar: ", res2$cindxLst$iskarCindex, "\n superPred: ", res2$cindxLst$superPredCindex)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res2$pVals$intgrStrcPVal,"\n perturbation: ", res2$pVals$intgrPertPVal,
    "\n sensitivity: ", res2$pVals$intgrSensPVal, "\n Iorio: ", res2$pVals$intgrIorioPVal, "\n Iskar: ", res2$pVals$intgrIskarPVal, "\n superPred: ", res2$pVals$intgrSuperPVal)

## ROC and PR plots
generateRocPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)_combined_without_gdsc")
generatePRPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)_combined_without_gdsc")


## generate communities
load("Output/gmt_targ_combined_with_gdsc.RData") # Files is generated in drugTargetBench.R
communityGenModded(integrtStrctSensPert, "ctrpv2", GMT_TARG, "communities_combined_with_gdsc.csv")


