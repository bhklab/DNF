###############################################################################################################
## Function computes and compares concordance indices and p-values for the integrative method vs. a single layered
## method, e.g., structure
##
## input: 
##     allPairs: list of all pairs obtained for the benchmark, structure, sensitivity, perturbation, and integrative method 
##     
## output: 
##     two lists containing the "c-indices"	of comparing all layers with the benchmark and "p-vals" of comaprison between each single 
##     layer vs. the integration layer
##
## 
###############################################################################################################


CompConcordIndxFlexible <- function(allPairs, compute.p.val=TRUE)
{
    c.index.obj.list <- list()
    
    # Create a list of c index objects
    for (i in 1:length(allPairs)) {
        pair.name <- strsplit(names(allPairs)[i], "Pairs")[[1]][1]
        
        if (pair.name != "bench") {
            c.index.name <- paste(pair.name, "Cindex", sep="")
            
            temp <- Hmisc::rcorr.cens(x=allPairs[[i]][, paste("obs", pair.name, sep=".")], S=allPairs$benchPairs$bench)
            c.index.obj.list[[c.index.name]] <- temp            
        }
    }
    
    c.index.list <- list()
    
    # Use the c.index.obj.list to obtain the actual C index values from the objects
    for (i in 1:length(c.index.obj.list)) {
        c.index.list[[names(c.index.obj.list)[i]]] <- c.index.obj.list[[i]]["C Index"]
    }
    
    p.vals.list <- list()
    
    if (compute.p.val) {
        for (i in 1:length(c.index.obj.list)) {
            c.index.name <- names(c.index.obj.list)[i]
            
            if (c.index.name != "integrCindex") {
                prefix <- strsplit(c.index.name, "Cindex")[[1]][1]
                
                pair.name <- paste(prefix, "Pairs", sep="")
                col.name <- paste("obs", prefix, sep=".")
                
                temp <- cindexComp2(c.index.obj.list$integrCindex, c.index.obj.list[[c.index.name]],
                                    allPairs$integrPairs$obs.integr, allPairs[[pair.name]][col.name])
                temp <- temp$p.value
                
                p.val.name <- paste(prefix, "PVal", sep="")
                p.vals.list[[p.val.name]] <- temp
            }
        }        
    }
    
    ## return both lists of results
    r <- list(c.index.list=c.index.list, p.vals.list=p.vals.list) 
    
    return(r)
}


