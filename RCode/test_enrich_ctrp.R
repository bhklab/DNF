###########################################################

# Clustering to get populations/clusters for SENSITIVITY
ap.ctrp.sens <- apcluster(cor.ctrp,q=0.9) #### 
length(ap.ctrp.sens) #38 populations/clusters

# Clustering to get populations/clusters for STRUCTURE
ap.smiles.common <- apcluster(fp.sim, q=0.9) 
length(ap.smiles.common) #58 populations/clusters

#Clustering to get populations/clusters for PERTURBATION
ap.lincs <- apcluster(cormat.lincs,q=0.9) 
length(ap.lincs) #31 populations/clusters


###########################################################
###########################################################
## HYPERGEOMETRIC TESTS ON BENCHMARKS - SINGLE LAYERS
###########################################################
###########################################################
# list of Benchmarked-GMT
ll_list<-list(GMT_ATC,GMT_TARG)
names(ll_list)<-c("ATC","Targets")
#Sensitivity    #Structure    #Perturbation
ApClustersList<-list(ap.ctrp.sens,ap.smiles.common,ap.lincs)
names(ApClustersList)<-c("Sensitivity","Structure","Perturbation")
HypergeomCalculations<-list()

#Conduct analysis for every layer, for every benchmark
for(Layer in 1:length(ApClustersList))
{
  HypergeomCalculations[[Layer]]<-list()
  
  ApCluster<-ApClustersList[[Layer]]
  message("Layer Analyzed: ",names(ApClustersList)[Layer])
  
  for(Benchmarks in 1:length(ll_list))
  {
    ll1<-ll_list[[Benchmarks]]
    message("Benchmark Analyzed: ",names(ll_list)[Benchmarks])
    
    Bkrd_Layer<-names(unlist(ApCluster@clusters)) #All the drugs in the test Layer
    Bkrd_Benchmark<- c(as.character(unlist(ll1))) #All the drugs in the GMT benchmark
    #totalx<-unique(c(Bkrd_Benchmark,Bkrd_Layer)) #Union of both
    
    ## test hypergeometric against Benchmark gmt for the sensitivity data
    ## test enrichment for ATC and Drug Targets in the common 191 drugs from NCI60 based on sensitivity data
    names(ApCluster@clusters)<-names(ApCluster@exemplars)
    ClusterList <-ApCluster@clusters
    
    #Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs
    ## ClusterList<-lapply(ClusterList,function(x){names(x[])<-intersect(names(x[]),Bkrd_Benchmark)})
    
    
    
    ###### if intersect with background > 1 keep the whole population
    ###### if intersect with background > 1 keep the whole population, new sept 20
    ClusterList <- lapply(ClusterList,function(x){if (length(intersect(Bkrd_Benchmark,names(x[]))) >= 2) {names(x[])} else {NULL}})
    ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 
    ClusterListRefined <- ClusterList[ClusterByCondition]
    
    ##### new added on sept 19 3:32 am
    totalx <- union(unique(as.character(unlist(ll1))), unique(unlist(ClusterListRefined)))
    
    EnrichmentMatrix <- matrix(NA, nrow=length(ll1), ncol=length(ClusterListRefined), 
                               dimnames=list(paste(names(ll1), 1:length(ll1),sep="_"), paste(names(ClusterListRefined), 1:length(ClusterListRefined), sep="_")))
    
    ###### hypergeometric
    # Compare every Benchmark geneset from ll1 to every population from ll2
    for (i in 1:length(ll1)) { #For all genesets in benchmark (52 for ATC, 26 for Targets)
      for (j in 1:length(ClusterListRefined)) {
        unx <-  unique(as.character(unlist(ll1[[i]]))) #drugs per ATC geneset
        unx2 <- unique(unlist(ClusterListRefined[[j]])) #drugs per population
        l1 <- length(unx)
        l2 <- length(unx2)
        intxx <- length(intersect(unx, unx2)) #intersect drugs
        if(intxx>1){
          #minxx3 <- 1-phyper(intxx-1, l1, length(totalx)-l1, l2)
          minxx3 <- fisher.test(matrix(c(intxx,l1-intxx,l2-intxx,(length(totalx)-l1)-(l2-intxx)),nrow=2,ncol=2),alternative="greater")$p.value
          EnrichmentMatrix[i, j] <- minxx3 }
      }
    }
    
    EnrichmentMatrix <- apply(EnrichmentMatrix,2, function(x) ifelse(!is.na(x),x,1)) #All NA's get assigned to 1 temporarily
    # data.frame(V2=apply(EnrichmentMatrix,1,function(x){ min(x[!is.na(x)])}))
    
    ## get the minimal p value from each row/column combination
    Result.EnrichmentMatrix <- t(sapply(seq(nrow(EnrichmentMatrix)), function(i) {
      j <- which.min(EnrichmentMatrix[i,])
      c(paste(rownames(EnrichmentMatrix)[i], colnames(EnrichmentMatrix)[j], sep='_'), EnrichmentMatrix[i,j])}))
    
    colnames(Result.EnrichmentMatrix) <- c("Enriched_Geneset", "PValue")
    Result.EnrichmentMatrix <- data.frame(Result.EnrichmentMatrix)
    
    ## transform factor to numeric (had some issue with that, you can improve the code)
    Result.EnrichmentMatrix$PValue <- as.numeric(levels(Result.EnrichmentMatrix$PValue))[Result.EnrichmentMatrix$PValue]
    # Result.EnrichmentMatrix[Result.EnrichmentMatrix == 1]<-NA #Check for true results. Convert 1's to NA
    
    HypergeomCalculations[[Layer]][[Benchmarks]]<-Result.EnrichmentMatrix
    
  } #end of loop across benchmarks
} #end of loop across layers

###########################################################
###########################################################
## HYPERGEOMETRIC TESTS ON BENCHMARKS - MULTIPLE LAYERS WITH SNF
###########################################################
###########################################################
ll_list<-list(GMT_ATC,GMT_TARG)
names(ll_list)<-c("ATC","Targets")

P1<-P2<-c("Sensitivity","Structure","Perturbation")
CombinationTwoLayers<-combn(unique(c(P2, P1)),m=2) 
colnames(CombinationTwoLayers)<-c("SensitivityStructure","SensitivityPerturbation","StructurePerturbation")


Structure = fp.sim #### structure
Sensitivity = cor.ctrp #### sensitivity
Perturbation = cormat.lincs #### perturbation

for(combo in 1:4)
{
  HypergeomCalculations[[combo+3]]<-list()
  for(Counter in 1:length(ll_list))
  {
    message("Benchmark Analyzed: ",names(ll_list)[Counter])
    ll1<-ll_list[[Counter]]
    message("COMBO Pair: ",combo)
    
    if(combo < 4)
    {
      First<-get(CombinationTwoLayers[1,combo])
      Second<-get(CombinationTwoLayers[2,combo]) 
      
      #Get SNF calculation for two layers
      W<- SNF(list(First,Second), 20, 20)
      colnames(W) <- colnames(cor.ctrp)
      rownames(W) <- rownames(cor.ctrp)
    }
    
    if(combo == 4) #last layer
    {
      W<- SNF(list(Sensitivity,Structure,Perturbation), 20, 20)
      colnames(W) <- colnames(cor.ctrp)
      rownames(W) <- rownames(cor.ctrp)
    }
    
    #Do Clustering to get Population/Clusters
    ap.combo1 <- apcluster(W)
    
    Bkrd_Layer<-names(unlist(ap.combo1@clusters)) #All the drugs in the test layer after SNF
    Bkrd_Benchmark<-c(as.character(unlist(ll1))) #All the drugs in the GMT benchmark
    #totalx <- unique(c(Bkrd_Benchmark,Bkrd_Layer)) #Intersection of both
    
    names(ap.combo1@clusters)<-names(ap.combo1@exemplars)
    ClusterList <- ap.combo1@clusters
    
    #Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs
    #ClusterList<-lapply(ClusterList,function(x){names(x[])<-intersect(names(x[]),Bkrd_Benchmark)})
    
    
    
    
    ###### if intersect with background > 1 keep the whole population, new sept 20
    ClusterList <- lapply(ClusterList,function(x){if (length(intersect(Bkrd_Benchmark,names(x[]))) >= 2) {names(x[])} else {NULL}})
    ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 
    ClusterListRefined <- ClusterList[ClusterByCondition]
    
    ##### new added on sept 19 3:32 am
    totalx <- union(unique(as.character(unlist(ll1))), unique(unlist(ClusterListRefined)))
    
    EnrichmentMatrix <- matrix(NA, nrow=length(ll1), ncol=length(ClusterListRefined), 
                               dimnames=list(paste(names(ll1), 1:length(ll1),sep="_"), 
                                             paste(names(ClusterListRefined), 1:length(ClusterListRefined), sep="_")))
    
    ###### hypergeometric
    for (i in 1:length(ll1)) {
      for (j in 1:length(ClusterListRefined)) {
        unx <-  unique(as.character(unlist(ll1[[i]])))
        unx2 <- unique(unlist(ClusterListRefined[[j]]))
        inters <- intersect(unx, unx2)
        l1 <- length(unx)
        l2 <- length(unx2)
        intxx <- length(inters)
        if(intxx>1){
          #minxx3 <- 1-phyper(intxx-1, l1, length(totalx)-l1, l2)
          minxx3 <- fisher.test(matrix(c(intxx,l1-intxx,l2-intxx,(length(totalx)-l1)-(l2-intxx)),nrow=2,ncol=2),alternative="greater")$p.value
          EnrichmentMatrix[i, j] <- minxx3  }
      }
    }
    
    EnrichmentMatrix <- apply(EnrichmentMatrix,2, function(x) ifelse(!is.na(x),x,1)) #All NA's get assigned to 1 temporarily
    
    ###### get the minimal p value from each row/column combination (you can improve the code!)
    Result.EnrichmentMatrix <- t(sapply(seq(nrow(EnrichmentMatrix)), function(i) {
      j <- which.min(EnrichmentMatrix[i,])
      c(paste(rownames(EnrichmentMatrix)[i], colnames(EnrichmentMatrix)[j], sep='_'), EnrichmentMatrix[i,j])
    }))
    
    colnames(Result.EnrichmentMatrix) <- c("Enriched_Geneset", "PValue")
    Result.EnrichmentMatrix <- data.frame(Result.EnrichmentMatrix)
    
    ### transform factor to numeric 
    Result.EnrichmentMatrix[,2] <- as.numeric(levels(Result.EnrichmentMatrix[,2]))[Result.EnrichmentMatrix[,2]]
    
    HypergeomCalculations[[combo+3]][[Counter]]<-Result.EnrichmentMatrix
    
  }
}

names(HypergeomCalculations)<-c("Sensitivity","Structure","Perturbation","SensitivityStructure","SensitivityPerturbation","StructurePerturbation","SensitivityStructurePerturbation")

for(i in 1:length(HypergeomCalculations))
{
  names(HypergeomCalculations[[i]])<-c("ATC","Targets")
}



