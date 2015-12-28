

rm(list=ls())

source('~/RCode/DRUG_SNF/preprocessInput.R')
source('~/RCode/DRUG_SNF/sensitivityData.R')
source('~/RCode/DRUG_SNF/perturbationData.R')
source('~/RCode/DRUG_SNF/structureData.R')
source('~/RCode/DRUG_SNF/constStructureLayer.R')
source('~/RCode/DRUG_SNF/constSensitivityLayer.R')
source('~/RCode/DRUG_SNF/constPerturbationLayer.R')
source('~/RCode/DRUG_SNF/integrateStrctSensPert.R')
source('~/RCode/DRUG_SNF/drugTargetBench.R')
source('~/RCode/DRUG_SNF/ATCBench.R')
source('~/RCode/DRUG_SNF/generateDrugPairs.R')
source('~/RCode/DRUG_SNF/compConcordIndx.R')
source('~/RCode/DRUG_SNF/generateRocPlot.R')
source('~/RCode/DRUG_SNF/predPerf.R')


library(PharmacoGx)  ## NOTE: requires version 1.0.6
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
dim(cDrugs)  ## 353 X 5
intersec <- cDrugs


## loading and cleaning data for the layers 
sensData <- sensitivityData("nci60", cDrugs)
dim(sensData)  ##[1]  48 348
pertData <- perturbationData("lincs", cDrugs)
dim(pertData)  ##[1] 978 350
strcData <- structureData("lincs", cDrugs)  ## a vector
length(strcData) ## 353


#Reduce All Matrices to lowest common set of drugs across all 3
# Get 346 drugs now in the reduced sets
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                     sort(colnames(pertData))))
commonDrugs <- commonDrugs[commonDrugs!="EPOTHILONEA"] #Remove double match with EPOTHILONE

##Reduce the data to the set of common drugs across the three layers
strcData<- strcData[commonDrugs]
sensData <- sensData[,commonDrugs]
pertData<- pertData[,commonDrugs]
## filter out cDrugs "dataframe" as well accordingly
cDrugs <- cDrugs[cDrugs$pert_iname %in% commonDrugs, ]


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

# Sanity Check - should all have the same dimensions: 345x345
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)


## USING 3 DIFFERENT BENCHMARK SETS TO EVALUATE ###
## 1- CHMEMBL -> DRUG-TARGET
## loading and cleaning benchmark dataset
dataBench1 <- drugTargetBench("chembl", commonDrugs)
dim(dataBench1) ##[1] 89 89
pairs1 <- generateDrugPairs(dataBench1, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combination layer vs. a single layer (e.g., structure)
res1 <- compConcordIndx(pairs1, "structure")
paste("c.index, combination of layers (integrative method): ", res1$c1$c.index)
paste("c.index, structure layer only: ", res1$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs1, d1Name="nci60", d2Name="lincs", benchNam="drug-target(CHEMBL)")


## 2- CHMEMBL -> ATC
## loading and cleaning benchmark dataset
dataBench2 <- ATCBench("chembl", cDrugs)
dim(dataBench2) ##[1] 94 94
pairs2 <- generateDrugPairs(dataBench2, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2, "structure")
paste("c.index, combination of layers (integrative method): ", res2$c1$c.index)
paste("c.index, structure layer only: ", res2$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs2, d1Name="nci60", d2Name="lincs", benchNam="ATC(CHEMBL)")


## 3- STITICH -> DRUG-TARGET
## loading and cleaning benchmark dataset
dataBench3 <- drugTargetBench("stitch", cDrugs)
dim(dataBench3)  ##[1] 108 108
pairs3 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combination layer vs. a single layer (e.g., structure)
res3 <- compConcordIndx(pairs3, "structure")
paste("c.index, combination of layers (integrative method): ", res3$c1$c.index)
paste("c.index, structure layer only: ", res3$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs3, d1Name="nci60", d2Name="lincs", benchNam="drug-target(STITCH)")

