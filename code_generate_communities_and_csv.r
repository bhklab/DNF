library(apcluster)

################ This is an example for NCI60/L1000 combination please do the same for CTRPv2/L1000
############# all communities from nci60/l1000 combi
set.seed(12345)
apcomb <- apcluster(combi, q=0.9)
ll <- list()
for(i in 1:length(apcomb)){
  xx <- names(apcomb[[i]])
  ll[[i]] <- xx
}

indx <- sapply(ll, length)
#indx <- lengths(lst) 
res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
row.names(llx) <- names(apcomb@exemplars)
write.csv(llx," communities_combiNCI60.csv", row.names=TRUE)
dim(res)

##########
############ Generate a CSV file from GMT_TARGET (chembl if set is NCI60)

ll <- list()
for(i in 1:length(GMT_TARGchembl)){
  xx <- GMT_TARGchembl[[i]]
  ll[[i]] <- xx
}

indx <- sapply(ll, length)
#indx <- lengths(lst) 
res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
row.names(llx) <- unlist(names(GMT_TARGchembl))
write.csv(llx," GMT_targ_chembl.csv", row.names=TRUE)
dim(res)


############### Keep communities with at least 2 drugs showing a known mechanism of action from GMT

GMT_TARG2 <- c(as.character(unlist(GMT_TARGchembl)))

#Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs

Clust <- apcomb@clusters
ClusterList <- lapply(Clust,function(x){if (length(intersect(GMT_TARG2 ,names(x[]))) >= 2) {names(x[])} else {NULL}})
ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 
ClusterListRefined <- ClusterList[ClusterByCondition]

indx <- sapply(ClusterListRefined, length)
res <- as.data.frame(do.call(rbind,lapply(ClusterListRefined, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
write.csv(llx,"ClusterListRefined_combiNCI60.csv", row.names=TRUE)