

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

pert.file.name <- "Data/L1000_compound_signatures.RData"
sensitivity.file.names <- c("Data/combined_sens_adjusted_diag_datasets_with.RData", 
                            "Data/combined_sens_adjusted_diag_datasets_without.RData",
                            "Data/combined_sens_datasets_with.RData",
                            "Data/combined_sens_datasets_without.RData")

imaging.combinations <- list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE), c(FALSE, FALSE))
remove.list <- c("roc.file.name", "cDrugs", "sensData", "pertData", "strcData",
                 "commonDrugs", "luminexData", "imagingData", "strcAffMat",
                 "sensAffMat", "pertAffMat", "integrtStrctSensPert",
                 "luminexAffMat", "imagingAffMat", "dataBench", "pairs", "res")

benchmark.names <- c("ctrpv2", "combined")

for (i in 1:length(sensitivity.file.names)) {
    for (j in 1:length(benchmark.names)) {
        for (k in 1:length(imaging.combinations)) {
            Main(sensitivity.file.name = sensitivity.file.names[i],
                 pert.file.name = pert.file.name,
                 benchmark.name = benchmark.names[j],
                 use.luminex = imaging.combinations[[k]][1],
                 use.imaging = imaging.combinations[[k]][2]
            )            
        }
    }
}

Main <- function(sensitivity.file.name, pert.file.name, benchmark.name, use.luminex=FALSE, use.imaging=FALSE) {
    roc.file.name <- CreatROCFileName(sensitivity.file.name = sensitivity.file.name,
                                      use.luminex = use.luminex, use.imaging = use.imaging,
                                      benchmark.name = benchmark.name)
    
    cDrugs <- preprocessInputModded(dname="combined", "lincs", sensitivity.file.name)
    dim(cDrugs)  ##239 X 28 for just CTRPv2, 309 x 28 for combined data
    
    sensData <- sensitivityDataModded("combined", cDrugs, sensitivity.file.name)  ## 645 X 239
    dim(sensData) # 309 x 309 for combined data
    pertData <- perturbationDataModded(pert.file.name, cDrugs)  ## 978 X 239
    dim(pertData) # 978 x 237 for 
    strcData <- structureData("lincs", cDrugs)  ## a vector  --> 239 elemnts
    length(strcData)
    
    luminexData <- NULL
    imagingData <- NULL
    
    if (use.luminex == FALSE && use.imaging == FALSE) {
        commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                             sort(colnames(pertData))))        
    } else if (use.luminex == TRUE && use.imaging == FALSE) {
        luminexData <- LuminexData(cDrugs, badchars)
        dim(luminexData)        
        
        commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                             sort(colnames(pertData)), sort(colnames(luminexData))))
    } else if (use.imaging == TRUE && use.luminex == FALSE) {
        imagingData <- ImagingData(cDrugs, badchars)
        dim(imagingData)        
        
        commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                         sort(colnames(pertData)), sort(colnames(imagingData))))
    } else {
        luminexData <- LuminexData(cDrugs, badchars)
        dim(luminexData)        

        imagingData <- ImagingData(cDrugs, badchars)
        dim(imagingData)        
    
        commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                         sort(colnames(pertData)), sort(colnames(imagingData)),
                                         sort(colnames(luminexData))))
    }
    
    length(commonDrugs) ## 239 ..
    
    strcData<- strcData[commonDrugs] # 239 drugs
    sensData <- sensData[commonDrugs, commonDrugs] # 645 x 239 drugs
    pertData<- pertData[, commonDrugs] #978 genes x  239 
    
    strcAffMat <- constStructureLayer(strcData)
    sensAffMat <- constSensitivityLayerCombined(sensData)
    pertAffMat <- constPerturbationLayer(pertData)
    luminexAffMat <- NULL
    imagingAffMat <- NULL
    
    if (use.luminex == FALSE && use.imaging == FALSE) {
        integrtStrctSensPert <- integrateStrctSensPertModded(sensAffMat, strcAffMat, 
                                                             pertAffMat)   
    } else if (use.luminex == TRUE && use.imaging == FALSE) {
        luminexData <- luminexData[, commonDrugs]
        luminexAffMat <- constLuminexLayer(luminexData)
        integrtStrctSensPert <- integrateStrctSensPertModded(sensAffMat, strcAffMat, 
                                                             pertAffMat, luminexAff = luminexAffMat)
    } else if (use.imaging == TRUE && use.luminex == FALSE) {
        imagingData <- imagingData[, commonDrugs]
        imagingAffMat <- constImagingLayer(imagingData)
        integrtStrctSensPert <- integrateStrctSensPertModded(sensAffMat, strcAffMat, 
                                                             pertAffMat, imagingAff = imagingAffMat)
    } else {
        luminexData <- luminexData[, commonDrugs]
        luminexAffMat <- constLuminexLayer(luminexData)
        
        imagingData <- imagingData[, commonDrugs]
        imagingAffMat <- constImagingLayer(imagingData)
        
        integrtStrctSensPert <- integrateStrctSensPertModded(sensAffMat, strcAffMat, 
                                                             pertAffMat, imagingAff = imagingAffMat,
                                                             luminexAff = luminexAffMat)
    }
    
    load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
    load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
    load("Data/drugERankSimil-CTRPV2-kendall.Rdata")
    
    if (benchmark.name == "ctrpv2") {
        dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141
    } else if (benchmark.name == "combined") {
        dataBench <- drugTargetBenchModded("ctrpv",  colnames(sensData), 
                                           "temp.RData") # 141 x 141 drug-drug adjacency matrix --> 141
    }
    
    if (use.luminex == FALSE && use.imaging == FALSE && benchmark.name == "ctrpv2") {
        pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2)
        res <- compConcordIndx(pairs)
        
        cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res$cindxLst$integrCindex, "\n structure: ", res$cindxLst$structureLayerCindex,
            "\n perturbation: ",  res$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res$cindxLst$iorioCindex, 
            "\n Iskar: ", res$cindxLst$iskarCindex)
        cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
            "\n sensitivity: ", res$pVals$intgrSensPVal, "\n Iorio: ", res$pVals$intgrIorioPVal, "\n Iskar: ", res$pVals$intgrIskarPVal)
        
        generateRocPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchName="drug-target_combined_adjusted_diagonal_with",
                        roc.file.name, nrow(dataBench))
        
    } else {
        pairs <- generateDrugPairsModded(dataBench, strcAffMat, sensAffMat, 
                                         pertAffMat, integrtStrctSensPert, 
                                         luminexAff=luminexAffMat, imagingAff=imagingAffMat)    
        res <- compConcordIndxModded(pairs)
        
        cat("c.indexes values from each layer vs. the benchmark: \n integration: ", 
            res$cindxLst$integrCindex, "\n structure: ", res$cindxLst$structureLayerCindex,
            "\n perturbation: ",  res$cindxLst$perturbationLayerCindex, "\n sensitivity: ", 
            res$cindxLst$sensitivityLayerCindex, "\n luminex: ",
            res$cindxLst$luminexLayerCindex, "\n imaging: ", 
            res$cindxLst$imagingLayerCindex)
        cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", 
            res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
            "\n sensitivity: ", res$pVals$intgrSensPVal,
            "\n luminex: ", res$pVals$intgrLuminexPVal,
            "\n imaging: ", res$pVals$intgrImagingPVal)
        
        generateRocPlotModded(pairs, d1Name="ctrpv2", d2Name="lincs", benchName="drug-target_combined_adjusted_diagonal_with",
                              roc.file.name, nrow(dataBench))
    }
}

CreatROCFileName <- function(sensitivity.file.name, use.luminex, use.imaging, benchmark.name) {
    base.dir <- "Output/auc_plots"
    
    if(!file.exists(base.dir)) {
        dir.create(base.dir)
    }
    
    sensitivity.file.name <- strsplit(sensitivity.file.name, "/")[[1]][2]
    sensitivity.file.name <- strsplit(sensitivity.file.name, "[.]")[[1]][1]
    print(sensitivity.file.name)
    
    file.name <- paste(base.dir, sensitivity.file.name, sep="/")
    
    if (use.luminex == TRUE) {
        file.name <- paste(file.name, "luminex", sep="_")
    }
    
    if (use.imaging == TRUE) {
        file.name <- paste(file.name, "imaging", sep="_")
    }
    
    file.name <- paste(file.name, benchmark.name, sep="_")
    file.name <- paste(file.name, "pdf", sep=".")
}
