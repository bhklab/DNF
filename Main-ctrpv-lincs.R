#Newest Version of PharmacoGx pkg:
library(devtools)
#devtools::install_github(repo="bhaibeka/PharmacoGx-private", ref="public", 
#                         auth_token="3cfdf7dfa99718ac738896cb34fe3ec91d264c10", build_vignettes=FALSE)

setwd("/Users/DeenaGendoo/Desktop/drugSNF/")

rm(list=ls())

source('./RCode/DRUG_SNF/preprocessInput.R')
source('./RCode/DRUG_SNF/sensitivityData.R')
source('./RCode/DRUG_SNF/perturbationData.R')
source('./RCode/DRUG_SNF/structureData.R')
source('./RCode/DRUG_SNF/constStructureLayer.R')
source('./RCode/DRUG_SNF/constSensitivityLayer.R')
source('./RCode/DRUG_SNF/constPerturbationLayer.R')
source('./RCode/DRUG_SNF/integrateStrctSensPert.R')
source('./RCode/DRUG_SNF/drugTargetBench.R')
source('./RCode/DRUG_SNF/generateDrugPairs.R')
source('./RCode/DRUG_SNF/compConcordIndx.R')
source('./RCode/DRUG_SNF/generateRocPlot.R')
source('./RCode/DRUG_SNF/predPerf.R')
source('./RCode/DRUG_SNF/ATCBench.R')


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

# Find common drugs between CTRPV2 and LINCS dataset
commonDrugs <- preprocessInput(dname="ctrpv2", "lincs")
length(commonDrugs)  ##239

# Process Sensitivity, Perturbation, and Structure layers for set of common drugs
sensData <- sensitivityData("ctrpv2", commonDrugs)
pertData <- perturbationData("lincs", commonDrugs)
strcData <- structureData("lincs", commonDrugs)  ## a vector

# Re-extract common drugs (some drugs may have lost information) - can probably remove
#commonDrugs <- intersect(colnames(sensData), colnames(pertData))
#sensData <- sensData[,commonDrugs,drop=F]
#sensData <- sensData[,order(colnames(sensData)),drop=F]
##Reduce the data to the set of common drugs across the three layers
#sensData <- sensData[,order(colnames(sensData)),drop=F] #Do we need to keep this? Do we need to order the rest?

## Get the common drugs (237) among the 3 datasets/layers
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                      sort(colnames(pertData))))
length(commonDrugs) ##237

strcData<- strcData[commonDrugs] # 237 drugs
sensData <- sensData[,commonDrugs] # 645 x 237 drugs
pertData<- pertData[,commonDrugs] #978 genes x 237


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)
integrtStrctSensPert <- integrateStrctSensPert(sensAffMat, strcAffMat, pertAffMat)

## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv", commonDrugs) # 139 x 139 drug-drug adjacency matrix

## intersecting the SNF layers (ie, SNF adjacency matrix of each layer and the integration) with the benchmark 
## Returns: list of 5 containing scores of drug-drug pairs for each of the layers and the integration and the benchmark
pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)

## validation: 1) compare cindices of combiantion layer vs. a single layer (e.g., structure)
res <- compConcordIndx(pairs, "structure")
paste("c.index, combination of layers (integrative method): ", res$c1$c.index)
paste("c.index, structure layer only: ", res$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target")





## 2- CHMEMBL -> ATC
## loading and cleaning benchmark dataset
cDrugs<-as.data.frame(commonDrugs)
names(cDrugs)<-"pert_iname"
dataBench2 <- ATCBench("chembl", cDrugs)

dim(dataBench2) ##[1] 43 43
pairs2 <- generateDrugPairs(dataBench2, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2, "structure")
paste("c.index, combination of layers (integrative method): ", res2$c1$c.index)
paste("c.index, structure layer only: ", res2$c2$c.index)
## validation: 2) ROC plots
generateRocPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)")


