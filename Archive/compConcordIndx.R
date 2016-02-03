###############################################################################################################
## Function computes and compares concordance indices and p-values for the integrative method vs. a single layered
## method, e.g., structure
##
## input: 
##     allPairs: list of all pairs obtained for the benchmark, structure, sensitivity, perturbation, and integrative method 
##     singleLayerNam: name ("character") of a layer 
## output: 
##     list containing the c-indices	for two method and comaprison information
##
## 
###############################################################################################################


compConcordIndx <- function(allPairs, singleLayerNam)
{
  
    integrCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$integrPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                              surv.event=rep(1, nrow(allPairs$integrPairs)), method="noether")

   if (singleLayerNam == "structure") {
        singleLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$strcPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                          surv.event=rep(1, nrow(allPairs$strcPairs)), method="noether")
   }else if (singleLayerNam == "perturbation") {
     singleLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$pertPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                      surv.event=rep(1, nrow(allPairs$pertPairs)), method="noether")
   } else if (singleLayerNam == "sensitivity") {
     singleLayerCindex <- survcomp::concordance.index(x=1-as.numeric(allPairs$sensPairs[ , 3]), surv.time=as.numeric(allPairs$benchPairs[ , 3]), 
                                                      surv.event=rep(1, nrow(allPairs$sensPairs)), method="noether")
   } else { stop("Error!")}
    
    
   r <- list(c1=integrCindex, c2=singleLayerCindex, comp=cindex.comp(cindex1=integrCindex, cindex2=singleLayerCindex)) 

   return(r)
}