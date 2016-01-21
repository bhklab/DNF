#library(apcluster)


communityGen <- function(combi, dname, GMT_TARG) {

## all communities from integrative analysis results for nci60/l1000 or ctrpv2/l1000 

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
filename = paste(getwd(), "/Output/", "communities_combi_", dname, ".csv", sep="")
write.csv(llx, filename, row.names=TRUE)
dim(res)

ll <- list()
for(i in 1:length(GMT_TARG)){
   xx <- GMT_TARG[[i]]
   ll[[i]] <- xx
}

indx <- sapply(ll, length)
#indx <- lengths(lst) 
res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
row.names(llx) <- unlist(names(GMT_TARG))
#write.csv(llx," GMT_targ_chembl.csv", row.names=TRUE)
dim(res)


## Keep communities with at least 2 drugs showing a known mechanism of action from GMT
GMT_TARG2 <- c(as.character(unlist(GMT_TARG)))

#Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs
Clust <- apcomb@clusters
ClusterList <- lapply(Clust,function(x){if (length(intersect(GMT_TARG2 ,names(x[]))) >= 2) {names(x[])} else {NULL}})
ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 
ClusterListRefined <- ClusterList[ClusterByCondition]

indx <- sapply(ClusterListRefined, length)
res <- as.data.frame(do.call(rbind,lapply(ClusterListRefined, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
filename = paste(getwd(), "/Output/", "clusterListRefined_combi_", dname, ".csv", sep="")
write.csv(llx, filename, row.names=TRUE)
}

