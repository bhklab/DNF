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


compConcordIndxCombined <- function(allPairs)
{
    sensitivityLayerCindex <- Hmisc::rcorr.cens(x=allPairs$sensPairs$obs.sens, S=allPairs$benchPairs$bench)
    iorioCindex <- Hmisc::rcorr.cens(x=allPairs$iorioPairs$obs.iorio, S=allPairs$benchPairs$bench)  
    iskarCindex <- Hmisc::rcorr.cens(x=allPairs$iskarPairs$obs.iskar, S=allPairs$benchPairs$bench)  
    
    if (length(allPairs)==8) {
        
        superPredCindex <- Hmisc::rcorr.cens(x=allPairs$superPairs$obs.superPred, S=allPairs$benchPairs$bench)  
        cindxLst <- list(sensitivityLayerCindex=sensitivityLayerCindex['C Index'], iorioCindex=iorioCindex['C Index'], iskarCindex=iskarCindex['C Index'], superPredCindex=superPredCindex['C Index'])
    } else {
        cindxLst <- list(sensitivityLayerCindex=sensitivityLayerCindex['C Index'], iorioCindex=iorioCindex['C Index'], iskarCindex=iskarCindex['C Index'])
    }
    
    
    ## return both lists of results
    r <- list(cindxLst=cindxLst) 
    
    
    return(r)
}


