getTargetDiscrepanciesLenient <- function(d1, d2, common.drugs, discrepancy.counts, pair.name) 
    {
    for (k in 1:length(common.drugs)) {
        drug.name <- common.drugs[k]
        d1.subset <- d1[d1$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        d2.subset <- d2[d2$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        
        d1.diff.d2 <- setdiff(d1.subset[,1], d2.subset[,1])
        d2.diff.d1 <- setdiff(d2.subset[,1], d1.subset[,1])
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
            print("Targets in common:")
            print(intersect(d1.subset[,1], d2.subset[,1]))
            print(paste("Drug: ", drug.name, sep=""))
            print("Target Discrepancies: ")                    
        }
        
        if (length(d1.diff.d2) > 0) {
            print("Targets in d1 not in d2: ")
            print(d1.diff.d2)                    
        }
        
        if (length(d2.diff.d1) > 0) {
            print("Targets in d2 not in d1: ")
            print(d2.diff.d1)   
        }
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
        }
        
        d1.subset.trimmed <- substring(d1.subset[,1], 1, 3)
        d2.subset.trimmed <- substring(d2.subset[,1], 1, 3)
        
        target.intersection <- intersect(d1.subset.trimmed, d2.subset.trimmed)
        
        if (length(target.intersection) == 0) {
            discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            #discrepancy.counts[[pair.name]][[drug.name]] <- 
            #    list(dataset.1.targets=d1.subset[,1], dataset.2.targets=d2.subset[,1])
        }
    }
    
    discrepancy.counts
}

getTargetDiscrepanciesPathway <- function(d1, d2, common.drugs, 
                                          discrepancy.counts, pair.name, xx) {
    for (k in 1:length(common.drugs)) {
        drug.name <- common.drugs[k]
        d1.subset <- d1[d1$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        d2.subset <- d2[d2$MOLECULE_NAME == drug.name, "TARGET_NAME", drop=FALSE]
        
        d1.diff.d2 <- setdiff(d1.subset[,1], d2.subset[,1])
        d2.diff.d1 <- setdiff(d2.subset[,1], d1.subset[,1])
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
            print("Targets in common:")
            print(intersect(d1.subset[,1], d2.subset[,1]))
            print(paste("Drug: ", drug.name, sep=""))
            print("Target Discrepancies: ")                    
        }
        
        if (length(d1.diff.d2) > 0) {
            print("Targets in d1 not in d2: ")
            print(d1.diff.d2)                    
        }
        
        if (length(d2.diff.d1) > 0) {
            print("Targets in d2 not in d1: ")
            print(d2.diff.d1)   
        }
        
        if (length(d1.diff.d2) > 0 | length(d2.diff.d1) > 0) {
            print("________________________________________________________________")
        }
        
        d1.genes <- d1.subset[, 1]
        d1.genes <- d1.genes[!is.na(d1.genes)]
        d2.genes <- d2.subset[, 1]
        d2.genes <- d2.genes[!is.na(d2.genes)]
        
        d1.genes <- d1.genes[d1.genes %in% keys(org.Hs.egSYMBOL2EG)]
        d2.genes <- d2.genes[d2.genes %in% keys(org.Hs.egSYMBOL2EG)]
        
        if (length(d1.genes) == 0 | length(d2.genes) == 0) {
            discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            next()
        }
        
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
        
        d1.pathways <- AnnotationDbi::select(reactome.db,
                              keys=d1.entrez,
                              keytype="ENTREZID",
                              columns=c("ENTREZID","REACTOMEID","PATHNAME")
        )
        
        d2.pathways <- AnnotationDbi::select(reactome.db,
                             keys=d2.entrez,
                             keytype="ENTREZID",
                             columns=c("ENTREZID","REACTOMEID","PATHNAME")
        )
        
        pathway.intersection <- intersect(d1.pathways$REACTOMEID, d2.pathways$REACTOMEID)
        
        
        if (length(pathway.intersection) == 0) {
            #discrepancy.counts[[pair.name]] <- c(discrepancy.counts[[pair.name]], drug.name)
            discrepancy.counts[[pair.name]][[drug.name]] <- list(dataset.1.targets=)
        }
    }
    
    discrepancy.counts
}

