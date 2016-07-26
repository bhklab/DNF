######## get ATC enrichments for the single and combination layers

ttx <- matrix(NA)
for (i in 1:length(HypergeomCalculations)){
  tty <-  data.frame(HypergeomCalculations[[i]]$ATC[,2]) #Pvalues of ATCs
  ttx <- cbind(ttx,tty)
}
ttx <- ttx[,-1]
row.names(ttx) <- gsub("_.*","",HypergeomCalculations[[1]]$ATC[,1])
names(ttx) <- names(HypergeomCalculations)
ttx.atc <- ttx

##### correct p values to account for empty cases
ttx.atc <- apply(ttx.atc, 2, function(x) p.adjust(x, method="fdr"))
ttx.atc <- ttx.atc[,c("Sensitivity","Structure","Perturbation", "SensitivityStructurePerturbation"),drop=F]

ttx.atc <- -log10(ttx.atc)


########### plot ATC enrichment

library(gplots)
myCol <- c("#ffffff","#bdbdbd","#fbb4b9","#f768a1","#c51b8a","#7a0177") 
# Defining breaks for the color scale
myBreaks <- c(0, 1.3, 2.5, 4, 6, 8, 1)
pdf("ATC_LAYER_ENRICHMENT_CTRP.pdf", width = 30, height = 45)

hm <- heatmap.2(as.matrix(ttx.atc), scale="none", Rowv=T, Colv=T,col = myCol, sepwidth=c(0.01,0.01),
                sepcolor="black", dendrogram = "none", density.info = "none",
                colsep=1:ncol(ttx.atc),
                rowsep=1:nrow(ttx.atc),
                breaks = myBreaks, margins=c(16,16), cexRow=1.5, cexCol=1.2, key=FALSE, keysize=1.5,trace="none")
legend("topright", fill = myCol, cex=1.3,legend = c("0", "<= 1.3", "1.31 to 3.9", "4 to 5.9", "6 to 7.9", ">=8")) ### cutoffs for -log10 pval

dev.off()



######## Get Target enrichments for the single and combination layers

ttx <- matrix(NA)
for (i in 1:length(HypergeomCalculations)){
  tty <-  data.frame(HypergeomCalculations[[i]]$Targets[,2]) #Pvalues of ATCs
  ttx <- cbind(ttx,tty)
}
ttx <- ttx[,-1]
row.names(ttx) <- gsub("_.*","",HypergeomCalculations[[1]]$Targets[,1])
names(ttx) <- names(HypergeomCalculations)
ttx.targ <- ttx

##### correct p values to account for empty cases
ttx.targ <- apply(ttx.targ, 2, function(x) p.adjust(x, method="fdr"))
ttx.targ <- ttx.targ[,c("Sensitivity","Structure","Perturbation", "SensitivityStructurePerturbation"),drop=F]

ttx.targ <- -log10(ttx.targ)


########### plot Target enrichment

library(gplots)
myCol <- c("#ffffff","#bdbdbd","#fbb4b9","#f768a1","#c51b8a","#7a0177") 
# Defining breaks for the color scale
myBreaks <- c(0, 1.3, 2.5, 4, 6, 8, 1)
pdf("Target_LAYER_ENRICHMENT_CTRP.pdf", width = 30, height = 45)

hm <- heatmap.2(as.matrix(ttx.targ), scale="none", Rowv=T, Colv=T,col = myCol, sepwidth=c(0.01,0.01),
                sepcolor="black", dendrogram = "none", density.info = "none",
                colsep=1:ncol(ttx.targ),
                rowsep=1:nrow(ttx.targ),
                breaks = myBreaks, margins=c(16,16), cexRow=1.5, cexCol=1.2, key=FALSE, keysize=1.5,trace="none")
legend("topright", fill = myCol, cex=1.3,legend = c("0", "<= 1.3", "1.31 to 3.9", "4 to 5.9", "6 to 7.9", ">=8")) ### cutoffs for -log10 pval

dev.off()