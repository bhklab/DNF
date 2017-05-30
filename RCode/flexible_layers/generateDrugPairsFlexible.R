###############################################################################################################
## Function intersects the benchmark data and the layers of the network and then generates drug pairs for the 
## final (intersection) subset of drugs
##
## input: 
##     bench.data
##     bench.data
##     strcAff
##     sensAff
##     pertAff
##     integration  
##
## output: 
##     list of drug pairs from benchmark, structure, perturbation, sensitivtiy and the integration of all
##
## 
###############################################################################################################

GenerateDrugPairsFlexible <- function(bench.data, strc.aff=NULL, sens.aff=NULL, 
                                    pert.aff=NULL, integration=NULL, 
                                    luminex.aff=NULL, imaging.aff=NULL) {
    
    layers <- list(strcAff=strc.aff, sensAff=sens.aff, pertAff=pert.aff, integr=integration,
                   luminexAff=luminex.aff, imagingAff=imaging.aff)
    
    layers <- layers[!sapply(layers, is.null)]
    intx <- intersect(colnames(bench.data), colnames(layers[[1]]))
    
    data.list <- list()
    
    for (i in 1:length(layers)) {
        prefix <- strsplit(names(layers)[[i]], split="[[:upper:]]")[[1]][1]
        prefix <- paste0(toupper(substr(prefix, 1, 1)), substr(prefix, 2, nchar(prefix)))
        
        dataset.name <- paste("data", prefix, sep="")
        aff.mat <- layers[[i]]
        
        temp <- aff.mat[intx, intx] 
        temp[upper.tri(temp, diag=TRUE)] <- NA
        data.list[[dataset.name]] <- temp 
    }
    
    pairs.list <- list()
    
    bench.data[upper.tri(bench.data, diag=TRUE)] <- NA
    bench.pairs <- melt(bench.data)
    bench.pairs <- na.omit(bench.pairs)
    colnames(bench.pairs)[3] <- "bench"
    
    for (i in 1:length(data.list)) {
        prefix <- strsplit(names(data.list)[i], split="data")[[1]][2]
        prefix <- tolower(prefix)
        
        pairs.name <- paste(prefix, "Pairs", sep="")
        temp <- melt(data.list[[i]])
        temp <- na.omit(temp)
        colnames(temp)[3] <- paste("obs", prefix, sep=".")
        
        pairs.list[[pairs.name]] <-  temp
    }
    
    pairs.list[["benchPairs"]] <- bench.pairs
    
    return(pairs.list)
}