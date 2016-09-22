# Deena M.A. Gendoo
# Determining the Correlation between different types of data (single-layer drug taxonomy)
# Do this by computing the spearman correlation between all pairs of similarity matrices
# September 22, 2016
################################################################################################
################################################################################################

bp<-matrix(ncol = 1,nrow = 3)
bp[1]<-bp[2]<-bp[3]<-0.15

#CTRPv2-based DNF
################################
#load Structure Affinity Matrix
load("strcSimilarityMat-ctrpv2.RData")
#load Sensitivity Affinity Matrix
load("sensSimilarityMat-ctrpv2.RData")
#load Perturbation Affinity Matrix
load("pertSimilarityMat-ctrpv2.RData")

#compute correlations between pairs of data types
# Structure-Sensitivity
Struct.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Structure - Perturbation
Struct.Vs.Pert<-cor.test(strcAffMat,pertAffMat,method = "spearman",alternative = "two.sided")
# Perturbation - Sensitivity
Pert.Vs.Sens<-cor.test(pertAffMat,sensAffMat,method = "spearman",alternative = "two.sided")

CTRPv2_Cors<-c(Struct.Vs.Pert$estimate,Struct.Vs.Sens$estimate,Pert.Vs.Sens$estimate)
names(CTRPv2_Cors)<-c("Struct vs Pert","Struct vs Sens","Pert vs Sens")

CTRPv2_PVals<-c(prettyNum(Struct.Vs.Pert$p.value,digits=3),prettyNum(Struct.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Sens$p.value,digits=3))

barplot(CTRPv2_Cors,col = "#d9ef8b",las=2,ylim=c(0,1))
text(bp,labels = CTRPv2_PVals)


#NCI60-based DNF
################################
#load Structure Affinity Matrix
load("strcSimilarityMat-nci60.RData")
#load Sensitivity Affinity Matrix
load("sensSimilarityMat-nci60.RData")
#load Perturbation Affinity Matrix
load("pertSimilarityMat-nci60.RData")

#compute correlations between pairs of data types
# Structure-Sensitivity
Struct.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Structure - Perturbation
Struct.Vs.Pert<-cor.test(strcAffMat,pertAffMat,method = "spearman",alternative = "two.sided")
# Perturbation - Sensitivity
Pert.Vs.Sens<-cor.test(pertAffMat,sensAffMat,method = "spearman",alternative = "two.sided")

NCI60_Cors<-c(Struct.Vs.Pert$estimate,Struct.Vs.Sens$estimate,Pert.Vs.Sens$estimate)
names(NCI60_Cors)<-c("Struct vs Pert","Struct vs Sens","Pert vs Sens")

NCI60_PVals<-c(prettyNum(Struct.Vs.Pert$p.value,digits=3),prettyNum(Struct.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Sens$p.value,digits=3))

barplot(NCI60_Cors,col = "#fee08b",las=2,ylim=c(0,1))
text(bp,labels = NCI60_PVals)


# Plot all Correlations For Both DNF taxonomies
pdf("CorrelationBetweenSingleLayers.pdf",width = 10)
par(mfrow=c(1,2),mai=c(1.82,0.82,0.82,0.82))
barplot(CTRPv2_Cors,col = "#d9ef8b",las=2,main="CTRPv2",ylab = "Spearman Correlation",ylim=c(0,1))
text(bp,labels = CTRPv2_PVals)
barplot(NCI60_Cors,col = "#fee08b",las=2,main="NCI60",ylab = "Spearman Correlation",ylim=c(0,1))
text(bp,labels = NCI60_PVals,las=2)
dev.off()


