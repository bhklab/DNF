

predPerf <- function(layerPair, benchPair, plotType="ROC") {
  
  
  if (plotType == "ROC") {
    pred <- prediction(layerPair, benchPair)
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    auc <- unlist(slot(auc, "y.values"))
    return(list(auc=auc, pred=pred, perf=perf))
    
  }
  
  ## NOTE  
  if (plotType == "PR") {
       PR <- PRROC::pr.curve(scores.class0=layerPair, weights.class0=benchPair, curve=TRUE, 
                             max.compute=TRUE, min.compute=TRUE, rand.compute=TRUE)
       return(PR)
  }
  

}
