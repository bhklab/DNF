GetTargetDiscrepanciesLenient <- function(d1, d2, common.drugs, discrepancy.counts, pair.name) {
    # Given a pair of drug target datasets, calculates the number of drugs between the 
    # two datasets that have targets which disagree. This is based on the following:
    # for an arbitrary drug x, if dataset 1 says that the targets for drug x are
    # gene 1, gene 2, gene 3, and dataset 2 says that the targets for drug x
    # are gene 4, gene 5, then there is no overlap between the targets in the two 
    # datasets for drug x, so drug x is counted as a "disagreeing" drug. Overlap
    # between the drug targets is based on just the first 3 characters of the 
    # targets such that CDK4 and CDK5 would be considered to be overlapping.
    #
    # Args:
    #   d1: A dataframe for the first dataset. This has MOLECULE_NAME column
    #       and a TARGET_NAME column.
    #   d2: A dataframe for the second dataset. This has a MOLECULE_NAME column
    #       and a TARGET_NAME column.
    #   common.drugs: A charcter vector of drugs that are in common between d1 and d2.
    #   discrepancy.counts: A list of all dataset pairs. Each element in the list
    #                       is a character vector containing names of disagreeing drugs.
    #   pair.name: A string indicating the name of the current dataset pair.
    #
    # Returns:
    #   discrepancy.counts, a list of drug target pairs, modified to include the disagreeing
    #   drugs for the current dataset pair specified by pair.name
    
    for (k in 1:length(common.drugs)) {
        drug.name <- common.drugs[k]
        d1.subset <- d1[d1$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        d2.subset <- d2[d2$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        
        # d1.diff.d2 is the targets in dataset 1 that are not in dataset 2
        # d2.diff.d1 is the targets in dataset 2 that are not in dataset 1
        d1.diff.d2 <- setdiff(d1.subset[,1], d2.subset[,1])
        d2.diff.d1 <- setdiff(d2.subset[,1], d1.subset[,1])
        
        # Print the targets that are in the intersection between the two datasets
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
            print("Targets in common:")
            print(intersect(d1.subset[,1], d2.subset[,1]))
            print(paste("Drug: ", drug.name, sep=""))
            print("Target Discrepancies: ")                    
        }
        
        # Print the targets in dataset 1 that are not in dataset 2
        if (length(d1.diff.d2) > 0) {
            print("Targets in d1 not in d2: ")
            print(d1.diff.d2)                    
        }
        
        # Print the targets in dataset 2 that are not in dataset 1
        if (length(d2.diff.d1) > 0) {
            print("Targets in d2 not in d1: ")
            print(d2.diff.d1)   
        }
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
        }
        
        # Trim the target names to just the first 3 characters
        d1.subset.trimmed <- substring(d1.subset[,1], 1, 3)
        d2.subset.trimmed <- substring(d2.subset[,1], 1, 3)
        
        # Determine the intersection between targets
        target.intersection <- intersect(d1.subset.trimmed, d2.subset.trimmed)
        
        # If there are no targets in common betweent the two datasets for the given
        # drug, then add the current drug to the list of disagreeing drugs
        if (length(target.intersection) == 0) {
            discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            #discrepancy.counts[[pair.name]][[drug.name]] <- 
            #    list(dataset.1.targets=d1.subset[,1], dataset.2.targets=d2.subset[,1])
        }
    }
    
    discrepancy.counts
}

GetTargetDiscrepanciesPathway <- function(d1, d2, common.drugs, 
                                          discrepancy.counts, pair.name, xx) {
    # Given a pair of drug target datasets, calculates the number of drugs between the 
    # two datasets that have targets which disagree. This is based on the following:
    # for an arbitrary drug x, if dataset 1 says that the targets for drug x are
    # gene 1, gene 2, gene 3, and dataset 2 says that the targets for drug x
    # are gene 4, gene 5, first determine the pathways that each of these targets
    # belong to. If there is no common pathway between the two datasets, then drug x
    # is counted as a disagreeing drug.
    #
    # Args:
    #   d1: A dataframe for the first dataset. This has MOLECULE_NAME column
    #       and a TARGET_NAME column.
    #   d2: A dataframe for the second dataset. This has a MOLECULE_NAME column
    #       and a TARGET_NAME column.
    #   common.drugs: A charcter vector of drugs that are in common between d1 and d2.
    #   discrepancy.counts: A list of all dataset pairs. Each element in the list
    #                       is a character vector containing names of disagreeing drugs.
    #   pair.name: A string indicating the name of the current dataset pair.
    #
    # Returns:
    #   discrepancy.counts, a list of drug target pairs, modified to include the disagreeing
    #   drugs for the current dataset pair specified by pair.name
    for (k in 1:length(common.drugs)) {
        drug.name <- common.drugs[k]
        d1.subset <- d1[d1$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        d2.subset <- d2[d2$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        
        # d1.diff.d2 is the targets in dataset 1 that are not in dataset 2
        # d2.diff.d1 is the targets in dataset 2 that are not in dataset 1
        d1.diff.d2 <- setdiff(d1.subset[,1], d2.subset[,1])
        d2.diff.d1 <- setdiff(d2.subset[,1], d1.subset[,1])
        
        # Print the targets that are in the intersection between the two datasets
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
            print("Targets in common:")
            print(intersect(d1.subset[,1], d2.subset[,1]))
            print(paste("Drug: ", drug.name, sep=""))
            print("Target Discrepancies: ")                    
        }
        
        # Print the targets in dataset 1 that are not in dataset 2
        if (length(d1.diff.d2) > 0) {
            print("Targets in d1 not in d2: ")
            print(d1.diff.d2)                    
        }
        
        # Print the targets in dataset 2 that are not in dataset 1
        if (length(d2.diff.d1) > 0) {
            print("Targets in d2 not in d1: ")
            print(d2.diff.d1)   
        }
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
        }
        
        # Get the targets that are not NA
        d1.genes <- d1.subset[, 1]
        d1.genes <- d1.genes[!is.na(d1.genes)]
        d2.genes <- d2.subset[, 1]
        d2.genes <- d2.genes[!is.na(d2.genes)]
        
        # Subset targets to only those that have corresponding ENTREZ IDs
        d1.genes <- d1.genes[d1.genes %in% keys(org.Hs.egSYMBOL2EG)]
        d2.genes <- d2.genes[d2.genes %in% keys(org.Hs.egSYMBOL2EG)]
        
        # If one of the datasets doesn't have targets with an ENTREZ ID mapping,
        # then add the drug to the list of disagreeing drugs
        if (length(d1.genes) == 0 | length(d2.genes) == 0) {
            discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            next()
        }
        
        # Get the ENTREZ IDs for the relevant targets
        d1.entrez <- unlist(mget(d1.genes, org.Hs.egSYMBOL2EG))
        d2.entrez <- unlist(mget(d2.genes, org.Hs.egSYMBOL2EG))
        
        if (!(d1.entrez %in% keys(reactome.db))) {
            next()
        }
        
        if (!(d2.entrez %in% keys(reactome.db))) {
            next()
        }
        
        if (any(is.na(d1.entrez)) | any(is.na(d2.entrez))) {
            next()
        }
        
        # Get the pathways for the targets in dataset 1
        d1.pathways <- AnnotationDbi::select(reactome.db,
                              keys=d1.entrez,
                              keytype="ENTREZID",
                              columns=c("ENTREZID","REACTOMEID","PATHNAME")
        )
        
        # Get the pathways for the targets in dataset 2
        d2.pathways <- AnnotationDbi::select(reactome.db,
                             keys=d2.entrez,
                             keytype="ENTREZID",
                             columns=c("ENTREZID","REACTOMEID","PATHNAME")
        )
        
        # Determine the intersection between the pathway identifiers for dataset 1
        # and dataset 2
        pathway.intersection <- intersect(d1.pathways$REACTOMEID, d2.pathways$REACTOMEID)
        
        # If no common pathway between the targets in dataset 1 and the targets in dataset 2
        # then add the current drug to the list of disagreeing drugs
        if (length(pathway.intersection) == 0) {
            #discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            discrepancy.counts[[pair.name]][[drug.name]] <- list(dataset.1.targets=)
        }
    }
    
    discrepancy.counts
}

