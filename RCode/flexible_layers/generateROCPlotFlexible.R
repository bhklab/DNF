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

GenerateROCPlotFlexible <- function(allPairs, file.name, num.drugs) {
    name.mapping <- list(sens="Sensitivity", pert="Perturbation", strc="Structure",
                         integr="Integration", luminex="Luminex", imaging="Imaging")
    color.mapping <- list(sens="#41ab5d", pert="#2b83ba", strc="#d7191c",
                          integr="black", luminex="#ce0aad", imaging="#ffc300")
    pred.list <- list()
    
    for (i in 1:length(allPairs)) {
        prefix <- strsplit(names(allPairs)[i], "Pairs")[[1]][1]
        if (prefix != "bench") {
            prefix.capital <- paste0(toupper(substr(prefix, 1, 1)), substr(prefix, 2, nchar(prefix)))
            pred.name <- paste("pred", prefix.capital, sep="")
            col.name <- paste("obs", prefix, sep=".")
            
            temp <- predPerf(allPairs[[i]][col.name], allPairs$benchPairs$bench)
            pred.list[[pred.name]] <- temp   
        }
    }

    # windows(width = 5, height = 5)
    # changing params for the ROC plot - width, etc
    # pdf(file.name, width=5, height = 5)
    par(mar=c(5,5,2,2), xaxs = "i", yaxs = "i", cex.axis=1.3, cex.lab=1.4)
    
    aucs.list <- list()
    legend.vals <- c()
    fill.vals <- c()
    add.to.plot <- FALSE
    
    for (i in 1:length(pred.list)) {
        prefix <- strsplit(names(pred.list)[i], "pred")[[1]][2]
        prefix <- tolower(prefix)
        
        temp <- paste(c(paste(name.mapping[[prefix]], " = ", sep="")), round(pred.list[[i]]$auc, 3), sep="")
        legend.vals <- c(legend.vals, temp)
        fill.vals <- c(fill.vals, color.mapping[[prefix]])
        
        if (i > 1) {
            add.to.plot <- TRUE
        }
        
        plot(pred.list[[i]]$perf, col=color.mapping[[prefix]], lwd=2, add=add.to.plot)
    }

    rand <- paste(c("rand = "), 0.5, sep="")
    
    legend(0.5,0.4,legend.vals,
           border="white", cex=0.75,
           box.col = "white",fill=fill.vals)
    
    title(paste("AUC Benchmark Number of Drugs: ", num.drugs, sep=" "))
    
    abline(0,1, col = "gray")
    # dev.off()
    dev.copy2pdf(file=file.name, width=5, height=5)
    # dev.off()
}
