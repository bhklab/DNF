CreateLeaveOneOutNetworks <- function(integrated, correlation.matrices, all.drugs, data.bench) {
    c1 <- makeCluster(8, outfile="")
    registerDoParallel(c1)
    
    functions <- ListFunctions()
    
    loo.networks <- foreach(drug=unique(data.bench$MOLECULE_NAME), .export=functions) %dopar% {
        common.drugs.one.out <- setdiff(all.drugs, drug)
        
        temp.correlation.matrices <- lapply(correlation.matrices, function(x) {
            x[intersect(rownames(x), common.drugs.one.out), intersect(colnames(x), common.drugs.one.out)]
        })    
        
        net.1 <- IntegrateCorrelationMatrices(correlation.matrices, common.drugs.one.out)
        
        net.2 <- integrated
        net.2[rownames(net.1), colnames(net.1)] <- net.1
        
        net.2
    }
    
    stopCluster(c1)
    
    names(loo.networks) <- unique(data.bench$MOLECULE_NAME)
    
    return(loo.networks)
}