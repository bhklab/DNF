# Deena M.A. Gendoo
# Determining Community Enrichment for DNF communities, against Drug Targets and ATC classes. 


# integrtStrctSensPert - DNFmatrix
# GMT_TARG - GMT

library(pheatmap)
library(apcluster)

### Load DNF datasets
A=load("old_ctrpv2-Integrated.RData")
CTRP_DNF = get(A)

B=load("old_nci60-Integrated.RData")
NCI60_DNF = get(B)

### Load Benchmarks for Drug Targets
# CTRP - Internal benchmark (unchanged from previous iteration)
x=load("gmt_targ_ctrpv.RData")
CTRP_GMT_TARG = get(x)

# NCI60 - Uniprot instead of Chembl
# Previously, y=load("./Output/gmt_targ_chembl.RData")
y=load("gmt_targ_uniprot.RData")
NCI60_GMT_TARG= get(y)

### Load Benchmarks for ATC
# ATC - New Updated chembl as opposed to older chembl!
# NB: The GMT files (RData) for the ATC were externally renamed to reflect CTRP or NCI DNF

# Previously,z=load("./Output/gmt_atc_chembl_CTRP.RData")
z=load("gmt_atc_chembl-new_CTRP.RData")
CTRP_GMT_ATC = get(z)

# Previously,k=load("./Output/gmt_atc_chembl_NCI60.RData")
k=load("gmt_atc_chembl-new_NCI60.RData")
NCI60_GMT_ATC = get(k)

### THE FUNCTION!
communityEnrichment <- function(DNFmatrix, DrugDNFName,BenchmarkName, GMT) {

  HypergeomCalculations<-list()
  
  #Generate a list of drug communities from AP Clustering
  apcomb <- apcluster(DNFmatrix, q=0.9)
  DrugClusters <- list()
  for(i in 1:length(apcomb)){
    xx <- names(apcomb[[i]])
    DrugClusters[[i]] <- xx
  }
  #Add community numbers to the AP clusters
  names(DrugClusters)<-paste("Community_",rep(1:length(DrugClusters)),sep="")
  #Number of drugs per community
  indx <- sapply(DrugClusters, length)

  #Rename Gene lists of the Drug Targets GMT file
  if(BenchmarkName=="ATC"){names(GMT) <- paste("ATC",unlist(names(GMT)),sep="_")}
  if(BenchmarkName=="Target"){names(GMT) <- paste("Target",unlist(names(GMT)),sep="_")}
  
  #Check if the genesets from the Drug Clusters intersect with the drugs in the GMT
  # If not, filter out the non-intersecting drug communities
  ClusterList <-apcomb@clusters
  names(ClusterList)<-paste("Community_",rep(1:length(apcomb@clusters)),sep="")
  
  Bkrd_Drugs<-length((unique(rownames(DNFmatrix))))
  Bkrd_Benchmark<- length(unique(c(as.character(unlist(GMT))))) #All the drugs in the GMT benchmark, N=141
  Bkrd_Benchmark_Drugs<- unique(c(as.character(unlist(GMT)))) #All the drugs in the GMT benchmark, N=141  
  
  # Reduce the community to keep only the drugs, in each community, which intersect with the Bkrd Benchmark
    ClusterList <- lapply(ClusterList,function(x){if (length(intersect(Bkrd_Benchmark_Drugs,names(x[]))) >= 1) {intersect(Bkrd_Benchmark_Drugs,names(x[]))} else {NULL}})
    ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 

    EnrichmentMatrix <- matrix(NA, nrow=length(DrugClusters), ncol=length(GMT), 
                             dimnames=list(names(DrugClusters), names(GMT)))
  
  ###### hypergeometric
  for (i in 1:length(DrugClusters)) { #For every drug cluster/community
    for (j in 1:length(GMT)) { # For every geneset of drugs associated with a particular target
      
      DrugClus <- unique(unlist(DrugClusters[[i]])) #drugs per population
      DrugTarg <-  unique(as.character(unlist(GMT[[j]]))) #drugs per geneset
      
      LenClus <- length(DrugClus)
      LenTarg <- length(DrugTarg)
      
      ClusAndTarg <- length(intersect(DrugClus, DrugTarg)) #intersected drugs between the two sets (+ves)
      ClusAndNotTarg<-length(DrugClus)-ClusAndTarg #exist in drug cluster, but not found in target geneset
      TargAndNotClus<-length(DrugTarg)-ClusAndTarg #exist in target geneset, but not found in drug cluster
      NotTargAndNotClus<-Bkrd_Benchmark-(ClusAndTarg+ClusAndNotTarg+TargAndNotClus)

      #Conduct the test only if there is at least 1 drug overlap between a community and drug target set
      if(ClusAndTarg>=1){
        minxx3 <- fisher.test(matrix(c(ClusAndTarg,ClusAndNotTarg,TargAndNotClus,NotTargAndNotClus),
                                     nrow=2,ncol=2),alternative="greater")$p.value
        EnrichmentMatrix[i, j] <- minxx3 }
    }
  }
  
  # Replace the 'NA' hits (ie, non-overlaps) with p-value of 1
  EnrichmentMatrix[is.na(EnrichmentMatrix)]<-1
  # Do a correction for multiple testing
  EnrichmentMatrix_FDRcorrected<-apply(EnrichmentMatrix,1,p.adjust,method="fdr")
  
  # Find out how many communities have at least one calculated p-value
  # (ie, how many communities have calculations on thems)
  EnrichmentMatrix_FDRcorrected<-EnrichmentMatrix_FDRcorrected[!apply(EnrichmentMatrix_FDRcorrected,1,function(x){all(x==1)}),]
  EnrichmentMatrix_FDRcorrected<-EnrichmentMatrix_FDRcorrected[,!apply(EnrichmentMatrix_FDRcorrected,2,function(x){all(x==1)})]
  dim(EnrichmentMatrix_FDRcorrected) #rows=targets, columns = communities
  
  # FYI: -log10(0.01)=2
  pdf(paste("Community_Enrichment_",DrugDNFName,"_",BenchmarkName,".pdf",sep = ""),width = 12,height = 10,onefile=FALSE)
  pheatmap(-log10(t(EnrichmentMatrix_FDRcorrected)),cluster_rows = F,cluster_cols = F,
           color = c('#ffffff',rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))) #higher number = more significant
  
  dev.off()
} 

### GENERATE PLOTS

#CTRPv2 with Internal Drug Target Benchmark
communityEnrichment(DNFmatrix = CTRP_DNF,GMT = CTRP_GMT_TARG,DrugDNFName = "CTRP",BenchmarkName = "Target")
#NCI60 with Drug Target Benchmark (Chembl)
communityEnrichment(DNFmatrix = NCI60_DNF,GMT = NCI60_GMT_TARG,DrugDNFName = "NCI60",BenchmarkName = "Target")

# CTRP with ATC Benchmark (Chembl)
communityEnrichment(DNFmatrix = CTRP_DNF,GMT = CTRP_GMT_ATC,DrugDNFName = "CTRP",BenchmarkName = "ATC")
# NCI60 with ATC Benchmark (Chembl)
communityEnrichment(DNFmatrix = NCI60_DNF,GMT = NCI60_GMT_ATC,DrugDNFName = "NCI60",BenchmarkName = "ATC")


#pheatmap(EnrichmentMatrix,cluster_rows = F,cluster_cols = F,
#         color = rev(c("#ffffff","#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177","#49006a")))

#pheatmap(-log10(t(EnrichmentMatrix_FDRcorrected)),cluster_rows = F,cluster_cols = F,
#         color = colorRampPalette(c("white","white","blue", "skyblue", "skyblue", "yellow", "yellow", "yellow", "yellow", "orange", "orange", "orange", "red", "darkred"))(length(seq(-0.4,1,by=0.01))-1))
