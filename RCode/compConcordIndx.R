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
     

     iorioCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$iorio[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                           surv.event=rep(1, nrow(allPairs$iorio)), method="noether")
     
     iskarCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$iskar[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                surv.event=rep(1, nrow(allPairs$iskar)), method="noether")
     if (length(allPairs)==8) {
         superPredCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$superPred[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                          surv.event=rep(1, nrow(allPairs$superPred)), method="noether")
         cindxLst <- list(integrCindex=integrCindex, structureLayerCindex=structureLayerCindex, perturbationLayerCindex=perturbationLayerCindex, 
                          sensitivityLayerCindex=sensitivityLayerCindex , iorioCindex=iorioCindex, iskarCindex=iskarCindex, superCindex=superPredCindex)
     } else { 
       cindxLst <- list(integrCindex=integrCindex, structureLayerCindex=structureLayerCindex, perturbationLayerCindex=perturbationLayerCindex, 
                       sensitivityLayerCindex=sensitivityLayerCindex , iorioCindex=iorioCindex, iskarCindex=iskarCindex)
     }
     ## perform c-index comparison and return the p-vals
     intgrStrcPVal <- cindex.comp(integrCindex, structureLayerCindex)
     intgrPertPVal <- cindex.comp(integrCindex, perturbationLayerCindex)
     intgrSensPVal <- cindex.comp(integrCindex, sensitivityLayerCindex)
     intgrIorioPVal <- cindex.comp(integrCindex, iorioCindex)
     intgrIskarPVal <- cindex.comp(integrCindex, iskarCindex)
     if (length(allPairs)==8) {
         intgrSuperPVal <- cindex.comp(integrCindex, superPredCindex)
         pVals <- list(intgrStrcPVal=intgrStrcPVal$p.value, intgrPertPVal=intgrPertPVal$p.value, 
                       intgrSensPVal=intgrSensPVal$p.value, intgrIorioPVal=intgrIorioPVal$p.value, intgrIskarPVal=intgrIskarPVal$p.value, 
                       intgrSuperPVal=intgrSuperPVal$p.value)
     } else { 
         pVals <- list(intgrStrcPVal=intgrStrcPVal$p.value, intgrPertPVal=intgrPertPVal$p.value, 
                   intgrSensPVal=intgrSensPVal$p.value, intgrIorioPVal=intgrIorioPVal$p.value, intgrIskarPVal=intgrIskarPVal$p.value) 
     }
     
     ## return both lists of results
     r <- list(cindxLst=cindxLst, pVals=pVals) 

   return(r)
}


