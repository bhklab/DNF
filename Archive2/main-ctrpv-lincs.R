

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
source('./RCode/generateDrugPairs.R')
source('./RCode/compConcordIndx.R')
source('./RCode/generateRocPlot.R')
source('./RCode/predPerf.R')
source('./RCode/ATCBench.R')
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

library(apcluster)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
#outputDir <- paste(getwd(), "/outputDir", sep="")
if ( ! file.exists("Output")) {
    dir.create("Output")
}

# Find common drugs between CTRPV2 and LINCS dataset
cDrugs <- preprocessInput(dname="ctrpv2", "lincs")
dim(cDrugs)  ##239 X 28

# Process Sensitivity, Perturbation, and Structure layers for set of common drugs
sensData <- sensitivityData("ctrpv2", cDrugs)  ## 645 X 239
dim(sensData)
pertData <- perturbationData("lincs", cDrugs)  ## 978 X 237
dim(pertData)
strcData <- structureData("lincs", cDrugs)  ## a vector  --> 239 elemnts
length(strcData)

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
save(integrtStrctSensPert, file="Data/ctrpv2-Integrated.RData")

## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv", commonDrugs) # 139 x 139 drug-drug adjacency matrix
## intersecting the SNF layers (ie, SNF adjacency matrix of each layer and the integration) with the benchmark 
## Returns: list of 5 containing scores of drug-drug pairs for each of the layers and the integration and the benchmark
pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combiantion layer vs. a single layer (e.g., structure)
res <- compConcordIndx(pairs)
cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res$cindxLst$integrCindex$c.index, "\n structure: ", res$cindxLst$structureLayerCindex$c.index,
    "\n perturbation: ",  res$cindxLst$perturbationLayerCindex$c.index, "\n sensitivity: ", res$cindxLst$sensitivityLayerCindex$c.index)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
    "\n sensitivity: ", res$pVals$intgrSensPVal)
## validation: 2) ROC plots
generateRocPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target")


cDrugs<-as.data.frame(commonDrugs)
names(cDrugs)<-"pert_iname"
## 2- CHEMBL -> ATC
## loading and cleaning benchmark dataset
dataBench2 <- ATCBench("chembl", cDrugs)
dim(dataBench2) ##[1] 43 43
pairs2 <- generateDrugPairs(dataBench2, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert)
## validation: 1) compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2)
cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res2$cindxLst$integrCindex$c.index, "\n structure: ", res2$cindxLst$structureLayerCindex$c.index,
    "\n perturbation: ",  res2$cindxLst$perturbationLayerCindex$c.index, "\n sensitivity: ", res2$cindxLst$sensitivityLayerCindex$c.index)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res2$pVals$intgrStrcPVal,"\n perturbation: ", res2$pVals$intgrPertPVal,
    "\n sensitivity: ", res2$pVals$intgrSensPVal)
## validation: 2) ROC plots
generateRocPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)")


## generate communities
load("Output/gmt_targ_ctrpv.RData")
communityGen(integrtStrctSensPert, "ctrpv2", GMT_TARG)


