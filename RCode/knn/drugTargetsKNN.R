###############################################################################################################
## Function reads in the benchmark sets and filter them out according to the intersection of the input datasets
## 
## input: 
##     cdrugs: a vector of common drugs between the input datasets
## output: 
##        a drug x drug adjacency/similarity matrix
##
## 
###############################################################################################################


DrugTargetsKNN <- function(cdrugs, gmt_file_name="communities_flexible", use.ctrpv2=FALSE,
                                    use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE) {
    drug.targets <- GetDrugTargetsFromDatasets(cdrugs=cdrugs, use.ctrpv2=use.ctrpv2, use.clue=use.clue,
                                               use.chembl=use.chembl, use.dbank=use.dbank,
                                               use.dtc=use.dtc)
    
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


