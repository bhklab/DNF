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


compConcordIndxModded <- function(allPairs)
{
    integrCindex <- Hmisc::rcorr.cens(x=allPairs$integrPairs$obs.combiall, S=allPairs$benchPairs$bench)
    structureLayerCindex <- Hmisc::rcorr.cens(x=allPairs$strcPairs$obs.str, S=allPairs$benchPairs$bench) 
    perturbationLayerCindex <- Hmisc::rcorr.cens(x=allPairs$pertPairs$obs.pert, S=allPairs$benchPairs$bench)
    sensitivityLayerCindex <- Hmisc::rcorr.cens(x=allPairs$sensPairs$obs.sens, S=allPairs$benchPairs$bench)
    
    
    luminexLayerCindex <- NULL
    imagingLayerCindex <- NULL
    
    intgrLuminexPVal <- NULL
    intgrImagingPVal <- NULL
    
    if (!is.null(allPairs$luminexPairs$obs.luminex) && !is.null(allPairs$imagingPairs$obs.imaging)) {
        luminexLayerCindex <- Hmisc::rcorr.cens(x=allPairs$luminexPairs$obs.luminex, S=allPairs$benchPairs$bench)
        imagingLayerCindex <- Hmisc::rcorr.cens(x=allPairs$imagingPairs$obs.imaging, S=allPairs$benchPairs$bench)
        
        intgrLuminexPVal <- cindexComp2(integrCindex, luminexLayerCindex, allPairs$integrPairs$obs.combiall,
                                        allPairs$luminexPairs$obs.luminex)
        intgrImagingPVal <- cindexComp2(integrCindex, imagingLayerCindex, allPairs$integrPairs$obs.combiall,
                                        allPairs$imagingPairs$obs.imaging)
    } else if (!is.null(allPairs$luminexPairs$obs.luminex) && is.null(allPairs$imagingPairs$obs.imaging)) {
        luminexLayerCindex <- Hmisc::rcorr.cens(x=allPairs$luminexPairs$obs.luminex, S=allPairs$benchPairs$bench)
        
        intgrLuminexPVal <- cindexComp2(integrCindex, luminexLayerCindex, allPairs$integrPairs$obs.combiall,
                                        allPairs$luminexPairs$obs.luminex)
    } else if (!is.null(allPairs$imagingPairs$obs.imaging) && is.null(allPairs$luminexPairs$obs.luminex)) {
        imagingLayerCindex <- Hmisc::rcorr.cens(x=allPairs$imagingPairs$obs.imaging, S=allPairs$benchPairs$bench)
        
        intgrImagingPVal <- cindexComp2(integrCindex, imagingLayerCindex, allPairs$integrPairs$obs.combiall,
                                        allPairs$imagingPairs$obs.imaging)
    }
    
    cindxLst <- list(integrCindex=integrCindex['C Index'], structureLayerCindex=structureLayerCindex['C Index'], 
                     perturbationLayerCindex=perturbationLayerCindex['C Index'],
                     sensitivityLayerCindex=sensitivityLayerCindex['C Index'],
                     luminexLayerCindex=luminexLayerCindex['C Index'],
                     imagingLayerCindex=imagingLayerCindex['C Index'])
    
    ## perform c-index comparison and return the p-vals
    intgrStrcPVal <- cindexComp2(integrCindex, structureLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$strcPairs$obs.str)
    intgrPertPVal <- cindexComp2(integrCindex, perturbationLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$pertPairs$obs.pert)
    intgrSensPVal <- cindexComp2(integrCindex, sensitivityLayerCindex, allPairs$integrPairs$obs.combiall, allPairs$sensPairs$obs.sens)
    
    
    pVals <- list(intgrStrcPVal=intgrStrcPVal$p.value, intgrPertPVal=intgrPertPVal$p.value,
                  intgrSensPVal=intgrSensPVal$p.value, intgrLuminexPVal=intgrLuminexPVal$p.value,
                  intgrImagingPVal=intgrImagingPVal$p.value)

    
    
    ## return both lists of results
    r <- list(cindxLst=cindxLst, pVals=pVals) 
    
    
    return(r)
}


