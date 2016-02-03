#################
### Code to construct a network from community clustering and keeping exemplar drugs as nodes
#### please do the same for ctrp
#################
library(igraph)
library(netbiov)

############### transform the similarity network into a distance matrix where the diagonals are 0s and similarity values are transformed into shortest distances
combiall <- combi
diag(combiall) <- 0
cormd <- 0.1/combiall
diag(cormd) <- 0
colnames(cormd) <- colnames(combiall)
rownames(cormd) <- rownames(combiall)

##### if you have forget to run the apclust object
set.seed(12345)
apcomb <- apcluster(combi.ctrp, q=0.9)

###### names of the exemplar drugs
exemp <- names(apcomb@exemplars)

########## build the network
graph <- graph.adjacency(cormd, weighted=TRUE, mode="undirected")
dfd <- get.data.frame(graph, what="edges")
nexemp <- dfd$from %in% exemp   &  dfd$to %in%  exemp 
dfd <- dfd[nexemp,]

dataframe2adjacency<-function(dfd, names.var=NULL){
  if(is.null(names.var)){
    names.var <- sort(union(dfd[,1],dfd[,2]))
  }
  adj <- matrix(0,nrow=length(names.var),ncol=length(names.var),dimnames=list(names.var,names.var))
  for(i in 1:nrow(dfd)){
    # print(dfd[i,3])
    adj[as.character(dfd[i,1]),as.character(dfd[i,2])]<-dfd[i,3]
  }
  return(adj)
}
mynet <- dataframe2adjacency(dfd)
g <- graph.adjacency(mynet,weighted=T,mode="undirected")

mycols <- rep("darkgrey",nrow(mynet))

names(mycols) <- rownames(mynet) 
v.colors <- mycols

###### order the names of exemplars according to network g and assign size of nodes based on number of drugs in each community
names(apcomb@clusters) <- exemp
apres2 <- apcomb@clusters[order(names(apcomb@clusters))]
ll <- unlist(lapply(apres2, function(x) length(apres2[x])))

set.seed(12345)
pdf("plot_network_nci60_targchembl_exemp.pdf",height=10,width=10, useDingbats =F)
mst.plot(g, mst.edge.col="black", vertex.color=mycols, colors="gray95",tkplot=FALSE, bg="white", v.size=ll,
         layout.function=layout.kamada.kawai,e.arrow=0, e.size=dfd$weight <= 30,v.lab =names(apres2),
         v.lab.col="black", lab.dist=5, v.lab.cex=0.45)

dev.off()
