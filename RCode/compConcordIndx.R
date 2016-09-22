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
  
  
     integrCindex <- Hmisc::rcorr.cens(x=allPairs$integrPairs$obs.combiall, S=allPairs$benchPairs$bench)
     structureLayerCindex <- Hmisc::rcorr.cens(x=allPairs$strcPairs$obs.str, S=allPairs$benchPairs$bench) 
     perturbationLayerCindex <- Hmisc::rcorr.cens(x=allPairs$pertPairs$obs.pert, S=allPairs$benchPairs$bench)
     sensitivityLayerCindex <- Hmisc::rcorr.cens(x=allPairs$sensPairs$obs.sens, S=allPairs$benchPairs$bench)
     iorioCindex <- Hmisc::rcorr.cens(x=allPairs$iorioPairs$obs.iorio, S=allPairs$benchPairs$bench)  
     iskarCindex <- Hmisc::rcorr.cens(x=allPairs$iskarPairs$obs.iskar, S=allPairs$benchPairs$bench)  
       
     if (length(allPairs)==8) {
         
         superPredCindex <- Hmisc::rcorr.cens(x=allPairs$superPairs$obs.superPred, S=allPairs$benchPairs$bench)  
         cindxLst <- list(integrCindex=integrCindex['C Index'], structureLayerCindex=structureLayerCindex['C Index'], perturbationLayerCindex=perturbationLayerCindex['C Index'], 
                          sensitivityLayerCindex=sensitivityLayerCindex['C Index'], iorioCindex=iorioCindex['C Index'], iskarCindex=iskarCindex['C Index'], superPredCindex=superPredCindex['C Index'])
     } else {
         cindxLst <- list(integrCindex=integrCindex['C Index'], structureLayerCindex=structureLayerCindex['C Index'], perturbationLayerCindex=perturbationLayerCindex['C Index'],
                       sensitivityLayerCindex=sensitivityLayerCindex['C Index'], iorioCindex=iorioCindex['C Index'], iskarCindex=iskarCindex['C Index'])
     }
     
     ## perform c-index comparison and return the p-vals
     intgrStrcPVal <- cindexComp2(integrCindex, structureLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$strcPairs$obs.str)
     intgrPertPVal <- cindexComp2(integrCindex, perturbationLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$pertPairs$obs.pert)
     intgrSensPVal <- cindexComp2(integrCindex, sensitivityLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$sensPairs$obs.sens)
     intgrIorioPVal <- cindexComp2(integrCindex, iorioCindex, allPairs$integrPairs$obs.combiall, allPairs$iorioPairs$obs.iorio)
     intgrIskarPVal <- cindexComp2(integrCindex, iskarCindex, allPairs$integrPairs$obs.combiall, allPairs$iskarPairs$obs.iskar)
     if (length(allPairs)==8) {
           intgrSuperPVal <- cindexComp2(integrCindex, superPredCindex, allPairs$integrPairs$obs.combiall, allPairs$superPairs$obs.superPred)
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


