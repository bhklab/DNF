# This file contains the function that is responsible for extracting relevant drug targets for the 
# specified drugs in order to evaluate the performance of the KNN method.

DrugTargetsKNN <- function(common.drugs, gmt_file_name="temp", use.ctrpv2=FALSE,
                                    use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE) {
    drug.targets <- GetDrugTargetsFromDatasets(common.drugs=common.drugs, use.ctrpv2=use.ctrpv2, use.clue=use.clue,
                                               use.chembl=use.chembl, use.dbank=use.dbank,
                                               use.dtc=use.dtc)
    # Create a dataframe containing drugs and their corresponding targets. The dataframe is comprised 
    # of two columns: one column indicates drug names, and the other column indicates target names.
    #
    # Args:
    #   common.drugs: A character vector containing the drugs for which to obtain targets.
    #   gmt_file_name: Legacy.
    #   use.ctrpv2: A logical value indicating whether or not to include targets from the CTRPv2 dataset.
    #   use.clue. A logical value indicating whether or not to include targets from the clue.io repurposing drugs 
    #             dataset.
    #   use.chembl: A logical value indicating whether or not to include targets from the CHEMBL dataset.
    #   use.dbank: A logical value indicating whether or not to include targets from the drug bank (uniprot) dataset.
    #   use.dtc: A logical value indicating whether or not to include targets from the drug target commons dataset.
    list.of.targs <- list()
    for(targName in unique(drug.targets$TARGET_NAME)){
        targGmt <- with(drug.targets, MOLECULE_NAME[TARGET_NAME==targName]) 
        list.of.targs[[targName]] <- targGmt
    }
    
    common.targs <- sapply(list.of.targs, function(x) length(x) >= 2)
    GMT_TARG<- list.of.targs[common.targs]
    
    #Number of drugs with TARGETS: length(unique(chembl_Targets.common.all$MOLECULE_NAME))
    drug.targets <- unique(drug.targets)
    drug.targets <- drug.targets[!is.na(drug.targets[,2]),]
    drug.targets <- drug.targets[drug.targets$TARGET_NAME %in% names(GMT_TARG), ]

    return(drug.targets)
}


