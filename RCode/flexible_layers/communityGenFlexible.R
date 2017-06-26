CommunityGenFlexible <- function(integrated, GMT_TARG) {
    # Creates two csv files containing communities. The communities are
    # determied via Affinity Propagation clustering on the integration layer.
    # The second file created only differs from the first file in that
    # it makes use of the GMT_TARG file to remove communities which don't
    # have at least 2 drugs in them for which there is drug target info available.
    #
    # Args:
    #   integrated: A similarity matrix resulting from combining layers via SNF.
    #   GMT_TARG: A list where the name of each list element is a drug target
    #             and the value of each list element is a character vector
    #             of all drugs that have that particular target.
    
    # Apply AP clustering on the integrated layer
    apcomb <- apcluster(integrated, q=0.9)
    # Extract the drugs in from each cluster and put them in a list
    list.of.communities <- list()
    for(i in 1:length(apcomb)){
        xx <- names(apcomb[[i]])
        list.of.communities[[i]] <- xx
    }
    
    # Create the unrefined communities file which contains communities with 
    # drugs that don't necessarily have drug target info available. 
    # Compute the size of each community
    communities.sizes <- sapply(list.of.communities, length)
    # Create a data frame of communities from the list of communities created previously
    community.data.frame <- as.data.frame(do.call(rbind,lapply(list.of.communities, `length<-`,max(communities.sizes))))
    # Add a column indicating community number, and another column indicating the size of each
    # community to community.data.frame
    final.data.frame <- data.frame("population"=1:length(rownames(community.data.frame)), "number of drugs"=communities.sizes,community.data.frame)
    # Replace numerical row names by the examplar of each community
    row.names(final.data.frame) <- names(apcomb@exemplars)
    filename = paste(getwd(), "/Output/", "communities_combined", ".csv", sep="")
    write.csv(final.data.frame, filename, row.names=TRUE)
    dim(community.data.frame)
    
    # Iterate over the GMT_TARG file and add the drugs to a list where the names
    # of list elements are numerical instead of drug targets.
    #ll <- list()
    #for(i in 1:length(GMT_TARG)){
    #    xx <- GMT_TARG[[i]]
    #    ll[[i]] <- xx
    #}
    
    #indx <- sapply(ll, length)
    #indx <- lengths(lst) 
    #res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
    #llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
    #row.names(llx) <- unlist(names(GMT_TARG))
    #write.csv(llx," GMT_targ_chembl.csv", row.names=TRUE)
    #dim(res)
    
    # Every drug in GMT_TARG2 has drug target info
    GMT_TARG2 <- c(as.character(unlist(GMT_TARG)))
    
    # Iterate over the communities previously created by apcluster() and
    # determine which ones have at least 2 drugs in them for which
    # there is drug target info.
    clust <- apcomb@clusters
    cluster.list <- lapply(clust, function(x) {
        if (length(intersect(GMT_TARG2 ,names(x[]))) >= 2) {names(x[])} 
        else {NULL}
    })
    
    # Keep only the communities that have at least 2 drugs in them for which
    # there is drug target info.
    cluster.by.condition <- sapply(cluster.list, function(x) length(x) > 1) 
    cluster.list.refined <- cluster.list[cluster.by.condition]
    
    
    # Create the refined communities file which contains communities with 
    # with a minimum of two drugs that have drug target info available. 
    # Compute the size of each community
    community.sizes <- sapply(cluster.list.refined, length)
    # Create a data frame of communities from the cluser.list.refined list
    community.data.frame <- as.data.frame(do.call(rbind,lapply(cluster.list.refined, `length<-`,max(community.sizes))))
    # Add a column indicating community number, and another column indicating the size of each
    # community to community.data.frame
    final.data.frame <- data.frame("population"=1:length(rownames(community.data.frame)), "number of drugs"=community.sizes,community.data.frame)
    filename = paste(getwd(), "/Output/", "clusterListRefined_combined", ".csv", sep="")
    write.csv(final.data.frame, filename, row.names=TRUE)
}

