# This script is similar to the main-ctrpv-lincs.R script in the original version of DNF except that
# it is flexible in terms of what layers are being used. Arbitrary combinations of layers, as well
# as arbitrary combinations of drug target datasets can be used in the analysis.
EvaluateModelROC <- function(use.sensitivity, use.perturbation, use.structure, use.imaging, use.luminex, 
                 sensitivity.file.name="", pert.file.name="", lincs.meta=NULL,
                 atc.benchmark.name="chembl-new", compute.atc=FALSE, use.ctrpv2=TRUE, use.clue=FALSE,
                 use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE, create.communities=FALSE, 
                 base.dir="", pert.new=FALSE, common.drugs=NULL, compute.drug.target.bench=TRUE,
                 use.subsetted.pert=FALSE) {
    target.roc.file.name <- CreateTargetROCFileName(base.dir=base.dir, sensitivity.file.name=sensitivity.file.name, 
                                       use.sensitivity=use.sensitivity,
                                       use.perturbation=use.perturbation, use.structure=use.structure,
                                       use.luminex=use.luminex, use.imaging=use.imaging,
                                       use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl, use.dbank=use.dbank,
                                       use.dtc=use.dtc)
    
    atc.roc.file.name <- CreateATCROCFileName(base.dir=base.dir, sensitivity.file.name=sensitivity.file.name, 
                                           atc.benchmark.name=atc.benchmark.name, use.sensitivity=use.sensitivity,
                                           use.perturbation=use.perturbation, use.structure=use.structure,
                                           use.luminex=use.luminex, use.imaging=use.imaging)
    
    gmt.file.name <- CreateGMTFileName(base.dir=base.dir, use.sensitivity=use.sensitivity,
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
        if (pert.new) {
            pert.data <- PerturbationDataFlexibleCustomSig(pert.file.name)
        } else {
            pert.data <- PerturbationDataFlexible(pert.file.name, lincs.meta, use.subsetted.pert)  ## 978 X 239
            print(dim(pert.data)) # 978 x 237 for
            
            colnames(pert.data) <- toupper(colnames(pert.data))
            colnames(pert.data) <- gsub(badchars, "", colnames(pert.data))
        }
        
        pert.names <- colnames(pert.data)

        # saveRDS(pert.data, "Data/uploading_features/perturbation/pert_features.RData")
        
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

    if (is.null(common.drugs)) {
        common.drugs <- Reduce(intersect, Filter(Negate(is.null),layers))
    } 
    
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
        
        # saveRDS(strc.data, "Data/uploading_features/structure/structure_features.RData")
        
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
    
    # saveRDS(integrated, "integrated.RData")
    print("Integration done")
    
    # Compute P-Values and create an AUC plot for the drug target benchmark
    if (compute.drug.target.bench) {
        CompDrugTargetBenchmarkFlexible(common.drugs=common.drugs, gmt.file.name=gmt.file.name,
                                        strc.aff.mat=strc.aff.mat, sens.aff.mat=sens.aff.mat, pert.aff.mat=pert.aff.mat,
                                        integration=integrated, luminex.aff.mat=luminex.aff.mat,
                                        imaging.aff.mat=imaging.aff.mat, target.roc.file.name=target.roc.file.name,
                                        use.ctrpv2=use.ctrpv2, use.clue=use.clue, use.chembl=use.chembl, 
                                        use.dbank=use.dbank, use.dtc=use.dtc)
        
    }
    
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
    
    res <- list(integrated=integrated)
    return(res)
}
