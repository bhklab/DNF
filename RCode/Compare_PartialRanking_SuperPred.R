## this script integrates a partial rnaking method to compute drug-drug similarities for drugE-Rank results ...
## 

rm(list=ls())


args <- commandArgs(trailingOnly = TRUE)
print(length(args))
args


benchNam <- as.character(args[1])



if (benchNam=="ctrp") {
  load("atc-ctrpv-new.RData")
  dataBechATC <- dataBench3 ## 51
  load("superPredResults_ctrp.Rdata")
  
} else {
  load("atc-NCI60-new.RData")
  dataBechATC <- dataBench3 ##  72 for nci60
  load("superPredResults_nci60.Rdata") ## -> drgTargs
}


atcList <- unique(uniAtc$ATC_CODE) # --> 107 atc 
m <-length(atcList) ## 107 unique atc

## NOTE:  superPred's drugSet is not complete!
length(unique(uniAtc$MOLECULE_NAME)) ## 37 for ctrp (out of 51) and 62 for nci60 (out of 72)

superPredSimil1 <- matrix(nrow = length(colnames(dataBechATC)), ncol=length(colnames(dataBechATC)), 0)
rownames(superPredSimil1) <- colnames(superPredSimil1) <- colnames(dataBechATC)

superPredSimil2 <- matrix(nrow = length(colnames(dataBechATC)), ncol=length(colnames(dataBechATC)), 0)
rownames(superPredSimil2) <- colnames(superPredSimil2) <- colnames(dataBechATC)


for (i in 1:length(colnames(dataBechATC))) {
  items1 <- vector(mode="numeric", length = m) 
  names(items1) <- atcList
  x <- uniAtc[uniAtc$MOLECULE_NAME == colnames(dataBechATC)[i],]
  ix <- match(x$ATC_CODE, atcList)
  items1[ix] <- x$SCORE
  items1 <-  1- as.numeric(items1)
  r1 <- rank(items1)
  
  for (j in ((i+1):length(colnames(dataBechATC)))) {
    if (j<=length(colnames(dataBechATC))) {
      items2 <- vector(mode="numeric", length = m) 
      names(items2) <- atcList
      xx <- uniAtc[uniAtc$MOLECULE_NAME == colnames(dataBechATC)[j],]
      ixx <- match(xx$ATC_CODE, atcList)
      items2[ixx] <- xx$SCORE
      items2 <- 1- as.numeric(items2)
      r2 <- rank(items2)
      
      ## generate the ranked lists, order, and then rank
      #superPredSimil1[i,j] <-  superPredSimil1[j,i] <- cor(r1,r2,method="spearman")
      superPredSimil2[i,j] <-  superPredSimil2[j,i] <- cor(r1,r2,method= "kendall")
    }
  }
}

## set na values to zero!
superPredSimil1[is.na(superPredSimil1)] <- 0
superPredSimil2[is.na(superPredSimil2)] <- 0


if (benchNam=="ctrp") {
    save(superPredSimil2, file="superPredSimil-ctrp-kendall.Rdata")
}
else {
    #save(superPredSimil1, file="superPredSimil-nci60-spearman.Rdata")
    save(superPredSimil2, file="superPredSimil-nci60-kendall.Rdata")
}

