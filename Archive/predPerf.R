

predPerf <- function(layerPair, benchPair) {
  
  pred <- prediction(layerPair, benchPair)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  auc <- unlist(slot(auc, "y.values"))
  
  return(list(auc=auc, pred=pred, perf=perf))
  
}
