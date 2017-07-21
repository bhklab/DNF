# This script is similar to the main-ctrpv-lincs.R script in the original version of DNF except that
# it is flexible in terms of what layers are being used. Arbitrary combinations of layers, as well
# as arbitrary combinations of drug target datasets can be used in the analysis.

rm(list=ls())

source("RCode/flexible_layers/structureDataFlexible.R")
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/luminexDataFlexible.R")
source("RCode/flexible_layers/imagingDataFlexible.R")

source("RCode/flexible_layers/constSensitivityLayerFlexible.R")
source("RCode/flexible_layers/constStructureLayerFlexible.R")
source("RCode/flexible_layers/constLuminexLayerFlexible.R")
source("RCode/flexible_layers/constImagingLayerFlexible.R")
source("RCode/flexible_layers/constPerturbationLayerFlexible.R")

source("RCode/flexible_layers/integrateLayersFlexible.R")
source("RCode/flexible_layers/generateDrugPairsFlexible.R")
source("RCode/flexible_layers/drugTargetBenchFlexible.R")
source("RCode/flexible_layers/compConcordIndxFlexible.R")
source("RCode/flexible_layers/printCIndices.R")
source("RCode/flexible_layers/printPVals.R")
source("RCode/flexible_layers/generateROCPlotFlexible.R")
source("RCode/flexible_layers/atcBenchFlexible.R")
source("RCode/flexible_layers/compDrugTargetBenchmarkFlexible.R")
source("RCode/flexible_layers/computeATCBenchmarkFlexible.R")
source("RCode/flexible_layers/stackedLayerAnalysisHelpers.R")
source("RCode/flexible_layers/communityGenFlexible.R")

source("RCode/cindexComp2.R")
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

# Load lincs metadata and clean the pert_inames in it
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
sensitivity.file.name <- "Data/combined_sensitivity//combined_sens_iname_replaced.RData"
    
res <- Main(use.sensitivity = TRUE, use.perturbation=TRUE, use.structure = TRUE, 
     use.imaging = FALSE, use.luminex = FALSE, sensitivity.file.name = sensitivity.file.name,
     pert.file.name = pert.file.name, lincs.meta = lincs.meta, use.ctrpv2=TRUE,
     use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE,
     create.communities=FALSE)

Main <- function(use.sensitivity, use.perturbation, use.structure, use.imaging, use.luminex, 
                 sensitivity.file.name="", pert.file.name="", lincs.meta=NULL,
                 atc.benchmark.name="chembl-new", compute.atc=FALSE, use.ctrpv2=TRUE, use.clue=FALSE,
                 use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE, create.communities=FALSE) {
    target.roc.file.name <- CreateTargetROCFileName(sensitivity.file.name=sensitivity.file.name, 
                                       use.sensitivity=use.sensitivity,
                                       use.perturbation=use.perturbation, use.structure=use.structure,
                                       use.luminex=use.luminex, use.imaging=use.imaging,
                                       use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl, use.dbank=use.dbank,
                                       use.dtc=use.dtc)
    
    atc.roc.file.name <- CreateATCROCFileName(sensitivity.file.name=sensitivity.file.name, 
                                           atc.benchmark.name=atc.benchmark.name, use.sensitivity=use.sensitivity,
                                           use.perturbation=use.perturbation, use.structure=use.structure,
                                           use.luminex=use.luminex, use.imaging=use.imaging)
    
    gmt.file.name <- CreateGMTFileName(use.sensitivity=use.sensitivity,
                                       use.perturbation=use.perturbation, use.structure=use.structure,
                                       use.luminex=use.luminex, use.imaging=use.imaging,
                                       use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl, use.dbank=use.dbank,
                                       use.dtc=use.dtc)
    
    sens.data <- NULL
    pert.data <- NULL
    strc.data <- NULL
    luminex.data <- NULL
    imaging.data <- NULL
    
    if (use.sensitivity) {
        # Load in the sensitivity data. Note that unlike before, this is now
        # already a correlation matrix.
        sens.data <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
        dim(sens.data) # 309 x 309 for combined data    
        rownames(sens.data) <- toupper(rownames(sens.data))
        rownames(sens.data) <- gsub(badchars, "", rownames(sens.data))
        
        colnames(sens.data) <- toupper(colnames(sens.data))
        colnames(sens.data) <- gsub(badchars, "", colnames(sens.data))
    }

    if (use.perturbation) {
        # If using the perturbation layer, use the names from the signature
        # file as the names to be intersected with the sensitivity layer
        pert.data <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 239
        print(dim(pert.data)) # 978 x 237 for         
        
        colnames(pert.data) <- toupper(colnames(pert.data))
        colnames(pert.data) <- gsub(badchars, "", colnames(pert.data))
        pert.names <- colnames(pert.data)
        
        saveRDS(pert.data, "Data/uploading_features/perturbation/pert_features.RData")
    } else if (use.structure) {
        # If using the structure layer, use the pert_iname column from the LINCS
        # metadata file as the names to be intersected with the sensitivity layer.
        # The reason for this is that the LINCS metadata file also has a column
        # with SMILES in it which will later on be used to create the fingerprints.
        pert.names <- lincs.meta$pert_iname
    } else {
        # If not using the perturbation or structure layer, then there is no intersection
        # to be performed between sensitivity layer and perturbation or structure layer.
        pert.names <- NULL
    }

    if (use.luminex) {
        luminex.data <- LuminexDataFlexible(badchars)
    }
    
    if (use.imaging) {
        imaging.data <- ImagingDataFlexible(badchars)
    }
    
    # Find the common drugs between the selected layers
    layers <- list(sens.names = sort(colnames(sens.data)), pert.names=pert.names,
                   luminex.names = sort(colnames(luminex.data)), imaging.names = sort(colnames(imaging.data)))
    common.drugs <- Reduce(intersect, Filter(Negate(is.null),layers))
    print(length(common.drugs))
    
    # Subset the datasets for the different layers based on common drugs
    sens.data <- sens.data[common.drugs, common.drugs] # 645 x 239 drugs
    pert.data <- pert.data[, common.drugs] #978 genes x  239 
    luminex.data <- luminex.data[, common.drugs]
    imaging.data <- imaging.data[, common.drugs]
    
    lincs.meta.subset <- lincs.meta[match(common.drugs, lincs.meta$pert_iname),]
    lincs.meta.subset <- lincs.meta.subset[!is.na(lincs.meta.subset$X),]
    
    # Create placeholder variables for the affinity matrices
    sens.aff.mat <- NULL
    strc.aff.mat <- NULL
    pert.aff.mat <- NULL
    imaging.aff.mat <- NULL
    luminex.aff.mat <- NULL

    # Since structure data is ubiquitous, we use this layer at the very end, once we've
    # determined the intersection of drugs that is relevant.    
    if (use.structure) {
        strc.data <- StructureDataFlexible(lincs.meta.subset)  ## a vector  --> 239 elemnts
        length(strc.data)     
        
        saveRDS(strc.data, "Data/uploading_features/structure/structure_features.RData")
        
        strc.aff.mat <- ConstStructureLayerFlexible(strc.data)
    }
    
    if (use.sensitivity) {
        sens.aff.mat <- ConstSensitivityLayerFlexible(sens.data)   
    }
    
    if (use.perturbation) {
        pert.aff.mat <- ConstPerturbationLayerFlexible(pert.data)
    }
    
    if (use.luminex) {
        luminex.aff.mat <- ConstLuminexLayerFlexible(luminex.data)
    }
    
    if (use.imaging) {
        imaging.aff.mat <- ConstImagingLayerFlexible(imaging.data)
    }

    # Combine the selected layers via SNF
    integrated <- IntegrateLayersFlexible(sens.aff=sens.aff.mat, strc.aff=strc.aff.mat, pert.aff=pert.aff.mat, 
                                          luminex.aff=luminex.aff.mat, imaging.aff=imaging.aff.mat)
    
    saveRDS(integrated, "integrated.RData")
    print("Integration done")
    
    # Compute P-Values and create an AUC plot for the drug target benchmark
    CompDrugTargetBenchmarkFlexible(common.drugs=common.drugs, gmt.file.name=gmt.file.name,
                            strc.aff.mat=strc.aff.mat, sens.aff.mat=sens.aff.mat, pert.aff.mat=pert.aff.mat,
                            integration=integrated, luminex.aff.mat=luminex.aff.mat,
                            imaging.aff.mat=imaging.aff.mat, target.roc.file.name=target.roc.file.name,
                            use.ctrpv2=use.ctrpv2, use.clue=use.clue, use.chembl=use.chembl, 
                            use.dbank=use.dbank, use.dtc=use.dtc)
    if (compute.atc) {
        # Compute P-Values and create an AUC plot for the ATC benchmark
        ComputeATCBenchmarkFlexible(atc.benchmark.name=atc.benchmark.name, common.drugs=common.drugs,
                                    strc.aff.mat=strc.aff.mat, sens.aff.mat=sens.aff.mat, pert.aff.mat=pert.aff.mat,
                                    integration=integrated, luminex.aff.mat=luminex.aff.mat, 
                                    imaging.aff.mat=imaging.aff.mat, atc.roc.file.name=atc.roc.file.name)
    }
    
    if (create.communities) {
        # Use Affinity Propagation clustering to determine the clusters
        # formed by the integrated layer.
        load(gmt.file.name)
        CommunityGenFlexible(integrated, GMT_TARG)
    }
}
