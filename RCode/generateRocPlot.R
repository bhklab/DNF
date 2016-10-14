###############################################################################################################
## Function generates ROC plots for every single layer separately and for the integration of all layers
## 
## input: 
##     allPairs: list of all pairs obtained for the benchmark, structure, sensitivity, perturbation, and integrative method 
##     d1Name: name of the dataset which is used for sensitivity information, e.g., "nci60" 
##     d2Name: currently must be set to "lincs"
##     benchName : name of the benchmark set, e.g., "stitch", "chembl"
## output: 
##     
##
## 
###############################################################################################################

generateRocPlot <- function(allPairs, d1Name, d2Name="lincs", benchName) {
  
  
  predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench)
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench)
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench)
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench)
  iorio <- predPerf(allPairs$iorio$iorio, allPairs$benchPairs$bench)
  iskar <- predPerf(allPairs$iskar$iskar, allPairs$benchPairs$bench)
  
  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
      predSuper <- predPerf(allPairs$superPairs$obs.superPred, allPairs$benchPairs$bench)
  }
  
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
      predDrugE <- predPerf(allPairs$drugePairs$obs.drugerank, allPairs$benchPairs$bench)
  }
  
  # changing params for the ROC plot - width, etc
  filename = paste(getwd(), "/Output/", "ROC_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
  pdf(filename, width=5, height = 5)
  par(mar=c(5,5,2,2), xaxs = "i", yaxs = "i", cex.axis=1.3, cex.lab=1.4)
  # plotting the ROC curve
  plot(predIntegr$perf, col="black", lwd=2)
  plot(predStrc$perf, col="#d7191c", lwd=2,add = TRUE)
  plot(predSens$perf, col = "#41ab5d", lwd=2,add = TRUE)
  plot(predPert$perf, col = "#2b83ba", lwd=2,add = TRUE)
  plot(iorio$perf, col = "pink", lwd=2,add = TRUE)
  plot(iskar$perf, col = "purple", lwd=2,add = TRUE)
  
  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
      plot(predSuper$perf, col = "orange", lwd=2,add = TRUE)
  }
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
      plot(predDrugE$perf, col = "orange", lwd=2,add = TRUE)
  }
  
  aucLegIntegr <- paste(c("Integration = "), round(predIntegr$auc,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc,3), sep="")
  aucLegIorio <- paste(c("IorioPGX = "), round(iorio$auc,3), sep="")
  aucLegIskar <- paste(c("Iskar = "), round(iskar$auc,3), sep="")
  
  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
      aucLegSuper<- paste(c("SuperPred = "), round(predSuper$auc,3), sep="")
      legend(0.5,0.4,c(aucLegIntegr, aucLegIorio, aucLegIskar, aucLegSuper), border="white", cex=0.75,
      box.col = "white",fill=c("black", "pink", "purple", "orange"))
  }
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
      aucLegDrugE <- paste(c("DrugERank = "), round(predDrugE$auc,3), sep="")
      legend(0.5,0.4,c(aucLegIntegr, aucLegIorio, aucLegIskar, aucLegDrugE), border="white", cex=0.75,
      box.col = "white",fill=c("black", "pink", "purple", "orange"))
  }
  
  abline(0,1, col = "gray")
  dev.off()
  
}
