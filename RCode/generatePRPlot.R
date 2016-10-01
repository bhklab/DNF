


generatePRPlot <- function(allPairs, d1Name, d2Name="lincs", benchName) {
  
  
  predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench, "PR")
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench, "PR")
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench, "PR")
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench, "PR")
  iorio <- predPerf(allPairs$iorio$iorio, allPairs$benchPairs$bench, "PR")
  iskar <- predPerf(allPairs$iskar$iskar, allPairs$benchPairs$bench, "PR")
  
  if (length(allPairs)==8) {
   super <- predPerf(allPairs$superPred$superPred, allPairs$benchPairs$bench, "PR")
  }
  
  # changing params for the ROC plot - width, etc
  filename = paste(getwd(), "/Output/", "NEW_PR_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
  pdf(filename, width=5, height = 5)
  par(mar=c(5,5,2,2), xaxs = "i", yaxs = "i", cex.axis=1.3, cex.lab=1.4)
  # plotting the ROC curve
  plot(predIntegr, col="black", lwd=2 , main="")
  plot(predStrc, col="#d7191c", lwd=2,add = TRUE)
  plot(predSens, col = "#41ab5d", lwd=2,add = TRUE)
  plot(predPert, col = "#2b83ba", lwd=2,add = TRUE)
  plot(iorio, col = "pink", lwd=2,add = TRUE)
  plot(iskar, col = "purple", lwd=2,add = TRUE)
  plot(predIntegr$rand, col = "gray", lwd=2,add = TRUE)
  
  if (length(allPairs)==8) {
    plot(super, col = "cyan", lwd=2,add = TRUE)
  }

  aucLegIntegr <- paste(c("Integration = "), round(predIntegr$auc.integral,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc.integral,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc.integral,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc.integral,3), sep="")
  aucLegIorio <- paste(c("IorioPGX = "), round(iorio$auc.integral,3), sep="")
  aucLegIskar <- paste(c("Iskar = "), round(iskar$auc.integral,3), sep="")
  
  rand <- paste(c("rand = "), round(predIntegr$rand$auc.integral,3), sep="")
  
  if (length(allPairs)==8) {
  aucLegSuper<- paste(c("SuperPred = "), round(super$auc.integral,3), sep="")
  legend(0.5,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper, rand), border="white", cex=0.75,
         box.col = "white",fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "cyan", "gray"))
  } else {
  legend(0.5,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, rand), border="white", cex=0.75,
         box.col = "white",fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "gray"))
  }
  
  
  dev.off()
  
  
}
