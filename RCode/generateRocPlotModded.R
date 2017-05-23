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

generateRocPlotModded <- function(allPairs, d1Name, d2Name="lincs", benchName, file.name, num.drugs) {
    
    
    predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench)
    predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench)
    predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench)
    predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench)
    predLuminex <- NULL
    predImaging <- NULL
    
    if (!is.null(allPairs$luminexPairs)) {
        predLuminex <- predPerf(allPairs$luminexPairs$obs.luminex, allPairs$benchPairs$bench)    
    }
    
    if (!is.null(allPairs$imagingPairs)) {
        predImaging <- predPerf(allPairs$imagingPairs$obs.imaging, allPairs$benchPairs$bench)    
    }
    
    # changing params for the ROC plot - width, etc
    #filename = paste(getwd(), "/Output/", "ROC_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
    #pdf(filename, width=5, height = 5)
    pdf(file.name, width=5, height = 5)
    
    par(mar=c(5,5,2,2), xaxs = "i", yaxs = "i", cex.axis=1.3, cex.lab=1.4)
    
    aucLegIntegr <- paste(c("Integration = "), round(predIntegr$auc,3), sep="")
    aucLegStr <- paste(c("Structure = "), round(predStrc$auc,3),sep="")
    aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc,3) , sep="")
    aucLegPer <- paste(c("Perturbation = "), round(predPert$auc,3), sep="")
    aucLegLuminex <- NULL
    aucLegImaging <- NULL
    
    # plotting the ROC curve
    plot(predIntegr$perf, col="black", lwd=2)
    plot(predStrc$perf, col="#d7191c", lwd=2,add = TRUE)
    plot(predSens$perf, col = "#41ab5d", lwd=2,add = TRUE)
    plot(predPert$perf, col = "#2b83ba", lwd=2,add = TRUE)
    
    if (!is.null(predLuminex)) {
        plot(predLuminex$perf, col="#ce0aad", lwd=2, add=TRUE)
        aucLegLuminex <- paste(c("Luminex = "), round(predLuminex$auc,3), sep="")
    }
    
    if (!is.null(predImaging)) {
        plot(predImaging$perf, col="#ffc300", lwd=2, add=TRUE)
        aucLegImaging <- paste(c("Imaging = "), round(predImaging$auc,3), sep="")
    }
    
    rand <- paste(c("rand = "), 0.5, sep="")
    
    legend(0.5,0.4,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, rand, aucLegLuminex, aucLegImaging),
           border="white", cex=0.75,
           box.col = "white",fill=c("black","#d7191c","#41ab5d","#2b83ba", "gray", "#ce0aad", "#ffc300"))
    
    title(paste("AUC Benchmark Number of Drugs: ", num.drugs, sep=" "))

    abline(0,1, col = "gray")
    dev.off()
    
}
