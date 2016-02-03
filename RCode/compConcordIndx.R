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


compConcordIndx <- function(allPairs)
{
  
     ## compute concordance-index for each single layer and integration vs. the benchmark
     integrCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$integrPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                              surv.event=rep(1, nrow(allPairs$integrPairs)), method="noether")
   
     structureLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$strcPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                          surv.event=rep(1, nrow(allPairs$strcPairs)), method="noether")
   
     perturbationLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$pertPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                      surv.event=rep(1, nrow(allPairs$pertPairs)), method="noether")
   
     sensitivityLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$sensPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                      surv.event=rep(1, nrow(allPairs$sensPairs)), method="noether")
   
     cindxLst <- list(integrCindex=integrCindex, structureLayerCindex=structureLayerCindex, perturbationLayerCindex=perturbationLayerCindex, 
                       sensitivityLayerCindex=sensitivityLayerCindex) 
     
     ## perform c-index comparison and return the p-vals
     intgrStrcPVal <- cindex.comp(integrCindex, structureLayerCindex)
     intgrPertPVal <- cindex.comp(integrCindex, perturbationLayerCindex)
     intgrSensPVal <- cindex.comp(integrCindex, sensitivityLayerCindex)
     pVals <- list(intgrStrcPVal=intgrStrcPVal$p.value, intgrPertPVal=intgrPertPVal$p.value, intgrSensPVal=intgrSensPVal$p.value)
     
     ## return both lists of results
     r <- list(cindxLst=cindxLst, pVals=pVals) 

     
     
   return(r)
}


