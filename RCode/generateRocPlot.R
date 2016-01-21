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
  
  # changing params for the ROC plot - width, etc
  filename = paste(getwd(), "/Output/", "ROC_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
  pdf(filename, width=5, height = 5)
  par(mar=c(5,5,2,2), xaxs = "i", yaxs = "i", cex.axis=1.3, cex.lab=1.4)
  # plotting the ROC curve
  plot(predIntegr$perf, col="black", lty=3, lwd=3)
  plot(predStrc$perf, col="#d7191c", lwd=2,add = TRUE)
  plot(predSens$perf, col = "#41ab5d", lwd=2,add = TRUE)
  plot(predPert$perf, col = "#2b83ba", lwd=2,add = TRUE)
  
  aucLegIntegr <- paste(c("Integration = "), round(predIntegr$auc,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc,3), sep="")
  
  legend(0.6,0.3,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer), border="white", cex=0.75, 
         box.col = "white",fill=c("black","#d7191c","#41ab5d","#2b83ba"))
  dev.off()
  
  
  
  
}
