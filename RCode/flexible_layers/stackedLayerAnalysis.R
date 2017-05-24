rm(list=ls())

source("RCode/flexible_layers/structureDataFlexible.R")
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/luminexDataFlexible.R")
source("RCode/flexible_layers/imagingDataFlexible.R")
source("RCode/flexible_layers/integrateLayersFlexible.R")
source("RCode/flexible_layers/generateDrugPairsFlexible.R")
source("RCode/flexible_layers/compConcordIndxFlexible.R")
source("RCode/flexible_layers/printCIndices.R")
source("RCode/flexible_layers/printPVals.R")
source("RCode/flexible_layers/generateROCPlotFlexible.R")

source("RCode/cindexComp2.R")
source("RCode/constImagingLayer.R")
source("RCode/constPerturbationLayer.R")
source("RCode/constSensitivityLayerCombined.R")
source("RCode/constStructureLayer.R")
source("RCode/constLuminexLayer.R")
source("RCode/drugTargetBenchModded.R")
source("RCode/drugTargetBench.R")
source("RCode/predPerf.R")

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
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
sensitivity.file.name <- "Data/combined_sens_adjusted_diag_datasets_with.RData"
    
res <- Main(use.sensitivity = TRUE, use.perturbation=FALSE, use.structure = TRUE, 
     use.imaging = TRUE, use.luminex = FALSE, sensitivity.file.name = sensitivity.file.name,
     pert.file.name = pert.file.name, lincs.meta = lincs.meta)

Main <- function(use.sensitivity, use.perturbation, use.structure, use.imaging, use.luminex, 
                 sensitivity.file.name="", pert.file.name="", lincs.meta=NULL, benchmark.name="ctrpv2") {
    roc.file.name <- CreateROCFileName(sensitivity.file.name=sensitivity.file.name, 
                                       benchmark.name=benchmark.name, use.sensitivity=use.sensitivity,
                                       use.perturbation=use.perturbation, use.structure=use.structure,
                                       use.luminex=use.luminex, use.imaging=use.imaging)
    sensData <- NULL
    pertData <- NULL
    strcData <- NULL
    luminexData <- NULL
    imagingData <- NULL
    
    if (use.sensitivity) {
        sensData <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
        dim(sensData) # 309 x 309 for combined data    
        rownames(sensData) <- toupper(rownames(sensData))
        rownames(sensData) <- gsub(badchars, "", rownames(sensData))
        
        colnames(sensData) <- toupper(colnames(sensData))
        colnames(sensData) <- gsub(badchars, "", colnames(sensData))
    }

    if (use.perturbation) {
        pertData <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 239
        print(dim(pertData)) # 978 x 237 for         
        
        colnames(pertData) <- toupper(colnames(pertData))
        colnames(pertData) <- gsub(badchars, "", colnames(pertData))
        pertNames <- colnames(pertData)
    } else {
        pertNames <- lincs.meta$pert_iname
    }

    if (use.luminex) {
        luminexData <- LuminexDataFlexible(badchars)
    }
    
    if (use.imaging) {
        imagingData <- ImagingDataFlexible(badchars)
    }
    
    layers <- list(sensNames = sort(colnames(sensData)), pertNames=pertNames,
                   luminexNames = sort(colnames(luminexData)), imagingNames = sort(colnames(imagingData)))

    commonDrugs <- Reduce(intersect, Filter(Negate(is.null),layers))
    
    # Subset the datasets for the different layers based on common drugs
    sensData <- sensData[commonDrugs, commonDrugs] # 645 x 239 drugs
    pertData <- pertData[, commonDrugs] #978 genes x  239 
    luminexData <- luminexData[, commonDrugs]
    imagingData <- imagingData[, commonDrugs]
    
    lincs.meta.subset <- lincs.meta[match(commonDrugs, lincs.meta$pert_iname),]
    lincs.meta.subset <- lincs.meta.subset[!is.na(lincs.meta.subset$X),]
    
    # Create placeholder variables for the affinity matrices
    sensAffMat <- NULL
    strcAffMat <- NULL
    pertAffMat <- NULL
    imagingAffMat <- NULL
    luminexAffMat <- NULL

    # Since structure data is ubiquitous, we use this layer at the very end, once we've
    # determined the intersection of drugs that is relevant.    
    if (use.structure) {
        strcData <- StructureDataFlexible(lincs.meta.subset)  ## a vector  --> 239 elemnts
        length(strcData)     
        strcAffMat <- constStructureLayer(strcData)
    }
    
    if (use.sensitivity) {
        sensAffMat <- constSensitivityLayerCombined(sensData)   
    }
    
    if (use.perturbation) {
        pertAffMat <- constPerturbationLayer(pertData)
    }
    
    if (use.luminex) {
        luminexAffMat <- constLuminexLayer(luminexData)
    }
    
    if (use.imaging) {
        imagingAffMat <- constImagingLayer(imagingData)
    }
    
    integrated <- IntegrateLayersFlexible(sensAff=sensAffMat, strcAff=strcAffMat, pertAff=pertAffMat, 
                                          luminexAff=luminexAffMat, imagingAff=imagingAffMat)
    
    if (benchmark.name == "ctrpv2") {
        dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141
    } else if (benchmark.name == "combined") {
        dataBench <- drugTargetBenchModded("ctrpv",  commonDrugs, 
                                           "temp.RData") # 141 x 141 drug-drug adjacency matrix --> 141
    }
    
    pairs <- GenerateDrugPairsFlexible(dataBench, strcAff=strcAffMat, sensAff=sensAffMat,
                                       pertAff=pertAffMat, integration=integrated,
                                       luminexAff=luminexAffMat, imagingAff=imagingAffMat)
    
    res <- CompConcordIndxFlexible(pairs)
    
    PrintCIndices(res$c.index.list)
    PrintPVals(res$p.vals.list)
    
    GenerateROCPlotFlexible(pairs, roc.file.name, nrow(dataBench))
}

CreateROCFileName <- function(sensitivity.file.name="", benchmark.name="", use.sensitivity=FALSE,
                             use.perturbation=FALSE, use.structure=FALSE, use.luminex=FALSE,
                             use.imaging=FALSE) {
    base.dir <- "Output/auc_p_flex"
    
    if(!file.exists(base.dir)) {
        dir.create(base.dir)
    }
    
    file.name <- paste(base.dir, "/", sep="")
    file.name <- paste(file.name, "sens", sep="")
    
    if (use.sensitivity) {
        sensitivity.file.name <- strsplit(sensitivity.file.name, "/")[[1]][2]
        sensitivity.file.name <- strsplit(sensitivity.file.name, "[.]")[[1]][1]        
        file.name <- paste(file.name, sensitivity.file.name, sep="_")
    }
    
    if (use.perturbation) {
        file.name <- paste(file.name, "pert", sep="_")
    }
    
    if (use.structure) {
        file.name <- paste(file.name, "strc", sep="_")
    }
    
    if (use.luminex) {
        file.name <- paste(file.name, "luminex", sep="_")
    }
    
    if (use.imaging) {
        file.name <- paste(file.name, "imaging", sep="_")
    }

    file.name <- paste(file.name, benchmark.name, sep="_")
    file.name <- paste(file.name, "pdf", sep=".")
}
