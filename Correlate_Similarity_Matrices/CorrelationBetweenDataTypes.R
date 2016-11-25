# Deena M.A. Gendoo
# Determining the Correlation between different types of data (single-layer drug taxonomy)
# Do this by computing the spearman correlation between all pairs of similarity matrices
# September 22, 2016
################################################################################################
################################################################################################

bp<-matrix(ncol = 1,nrow = 6)
bp[1]<-bp[2]<-bp[3]<-bp[4]<-bp[5]<-bp[6]<-0.15

#CTRPv2-based DNF
################################
#load Structure Affinity Matrix
load("strcSimilarityMat-ctrpv2.RData")
#load Sensitivity Affinity Matrix
load("sensSimilarityMat-ctrpv2.RData")
#load Perturbation Affinity Matrix
load("pertSimilarityMat-ctrpv2.RData")
#load DNF itself!
load("integrationSimilarityMat-ctrpv2.RData")

#compute correlations between pairs of data types
#Single Layers against themselves
# Structure-Sensitivity
Struct.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Structure - Perturbation
Struct.Vs.Pert<-cor.test(strcAffMat,pertAffMat,method = "spearman",alternative = "two.sided")
# Perturbation - Sensitivity
Pert.Vs.Sens<-cor.test(pertAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Single layers against DNF
# DNF- Pert
DNF.Vs.Pert<-cor.test(integrtStrctSensPert,pertAffMat,method = "spearman",alternative = "two.sided")
# DNF - Sens
DNF.Vs.Sens<-cor.test(integrtStrctSensPert,sensAffMat,method = "spearman",alternative = "two.sided")
# DNF - Struct
DNF.Vs.Struct<-cor.test(integrtStrctSensPert,strcAffMat,method = "spearman",alternative = "two.sided")

CTRPv2_Cors<-c(Struct.Vs.Pert$estimate,Struct.Vs.Sens$estimate,Pert.Vs.Sens$estimate,NA,
               DNF.Vs.Pert$estimate,DNF.Vs.Sens$estimate,DNF.Vs.Struct$estimate)
names(CTRPv2_Cors)<-c("Struct vs Pert","Struct vs Sens","Pert vs Sens","", 
                      "DNF vs Pert","DNF vs Sens", "DNF vs Struct")

CTRPv2_PVals<-c(prettyNum(Struct.Vs.Pert$p.value,digits=3),prettyNum(Struct.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Sens$p.value,digits=3),
                prettyNum(DNF.Vs.Pert$p.value,digits=3),prettyNum(DNF.Vs.Sens$p.value,digits=3),prettyNum(DNF.Vs.Struct$p.value,digits=3))

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
#load DNF itself!
load("integrationSimilarityMat-nci60.RData")

#compute correlations between pairs of data types
#Single Layers against themselves
# Structure-Sensitivity
Struct.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Structure - Perturbation
Struct.Vs.Pert<-cor.test(strcAffMat,pertAffMat,method = "spearman",alternative = "two.sided")
# Perturbation - Sensitivity
Pert.Vs.Sens<-cor.test(pertAffMat,sensAffMat,method = "spearman",alternative = "two.sided")
# Single layers against DNF
# DNF- Pert
DNF.Vs.Pert<-cor.test(integrtStrctSensPert,pertAffMat,method = "spearman",alternative = "two.sided")
# DNF - Sens
DNF.Vs.Sens<-cor.test(integrtStrctSensPert,sensAffMat,method = "spearman",alternative = "two.sided")
# DNF - Struct
DNF.Vs.Struct<-cor.test(integrtStrctSensPert,strcAffMat,method = "spearman",alternative = "two.sided")



NCI60_Cors<-c(Struct.Vs.Pert$estimate,Struct.Vs.Sens$estimate,Pert.Vs.Sens$estimate,NA,
              DNF.Vs.Pert$estimate,DNF.Vs.Sens$estimate,DNF.Vs.Struct$estimate)
names(NCI60_Cors)<-c("Struct vs Pert","Struct vs Sens","Pert vs Sens","", 
                     "DNF vs Pert","DNF vs Sens", "DNF vs Struct")

NCI60_PVals<-c(prettyNum(Struct.Vs.Pert$p.value,digits=3),prettyNum(Struct.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Sens$p.value,digits=3),
               prettyNum(DNF.Vs.Pert$p.value,digits=3),prettyNum(DNF.Vs.Sens$p.value,digits=3),prettyNum(DNF.Vs.Struct$p.value,digits=3))

barplot(NCI60_Cors,col = "#fee08b",las=2,ylim=c(0,1))
text(bp,labels = NCI60_PVals)


# Plot all Correlations For Both DNF taxonomies
pdf("CorrelationBetweenSingleLayers.pdf",width = 10)
par(mfrow=c(1,2),mai=c(1.82,0.82,0.82,0.82))
barplot(CTRPv2_Cors,col = "#d9ef8b",las=2,main="CTRPv2",ylab = "Spearman Correlation",ylim=c(0,1))
#text(bp,labels = CTRPv2_PVals)
barplot(NCI60_Cors,col = "#fee08b",las=2,main="NCI60",ylab = "Spearman Correlation",ylim=c(0,1))
#text(bp,labels = NCI60_PVals,las=2)
dev.off()


