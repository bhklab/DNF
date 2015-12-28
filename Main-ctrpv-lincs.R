


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
source('~/RCode/DRUG_SNF/generateDrugPairs.R')
source('~/RCode/DRUG_SNF/compConcordIndx.R')
source('~/RCode/DRUG_SNF/generateRocPlot.R')
source('~/RCode/DRUG_SNF/predPerf.R')


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
#outputDir <- paste(getwd(), "/outputDir", sep="")
if ( ! file.exists("Output")) {
    dir.create("Output")
}

commonDrugs <- preprocessInput(dname="ctrpv2", "lincs")
length(commonDrugs)  ##239

sensData <- sensitivityData("ctrpv2", commonDrugs)
pertData <- perturbationData("lincs", commonDrugs)
commonDrugs <- intersect(colnames(sensData), colnames(pertData))
length(commonDrugs) ##237
strcData <- structureData("lincs", commonDrugs)  ## a vector

sensData <- sensData[,commonDrugs,drop=F]
sensData <- sensData[,order(colnames(sensData)),drop=F]


## Get the common drugs (237) among the 3 datasets/layers
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                      sort(colnames(pertData))))
##Reduce the data to the set of common drugs across the three layers
strcData<- strcData[commonDrugs]
sensData <- sensData[,commonDrugs]
pertData<- pertData[,commonDrugs]


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)
integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)

## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv", commonDrugs) # 139
## intersecting the layers with the benchmark 
pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)

## validation: 1) compare cindices of combiantion layer vs. a single layer (e.g., structure)
res <- compConcordIndx(pairs, "structure")
paste("c.index, combination of layers (integrative method): ", res$c1$c.index)
paste("c.index, structure layer only: ", res$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target")



