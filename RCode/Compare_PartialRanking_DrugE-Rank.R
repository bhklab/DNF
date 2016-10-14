

rm(list=ls())

## this script integrates a partial rnaking method to compute drug-drug similarities for drugE-Rank results ...
## 


args <- commandArgs(trailingOnly = TRUE)
print(length(args))
args


benchNam <- as.character(args[1])


if (benchNam=="ctrp") {
  load("targetSET_CTRPV2.RData") 
  dataBenchTRG <- dataBench ## 141
  load("drugERankResults_ctrp.RData") 
} else {
  load("target-uniProt-nci60.RData")
  dataBenchTRG <- x ##  86 for nci60
  load("drugERankResults_nci60.RData")

}


targetList <- unique(drgTargs$TARGET_NAME)
m <-length(unique(drgTargs$TARGET_NAME)) ## 373 unique target ids
drgTargs$SCORE <- 21-(as.numeric(as.character(drgTargs$SCORE))*20)


drugERankSimil1 <- matrix(nrow = length(colnames(dataBenchTRG)), ncol=length(colnames(dataBenchTRG)), 0)
rownames(drugERankSimil1) <- colnames(drugERankSimil1) <- colnames(dataBenchTRG)

drugERankSimil2 <- matrix(nrow = length(colnames(dataBenchTRG)), ncol=length(colnames(dataBenchTRG)), 0)
rownames(drugERankSimil2) <- colnames(drugERankSimil2) <- colnames(dataBenchTRG)


for (i in 1:length(colnames(dataBenchTRG))) {
  items1 <- vector(mode="numeric", length = m) 
  names(items1) <- unique(drgTargs$TARGET_NAME)
  
  x <- drgTargs[drgTargs$MOLECULE_NAME == colnames(dataBenchTRG)[i],]
  ix <- match(targetList, x$TARGET_NAME)
  ix[is.na(ix)] <- 21
  items1 <- ix
  for (j in ((i+1):length(colnames(dataBenchTRG)))) {
    if (j<=length(colnames(dataBenchTRG))) {
         items2 <- vector(mode="numeric", length = m) 
         names(items2) <- unique(drgTargs$TARGET_NAME)
      
         xx <- drgTargs[drgTargs$MOLECULE_NAME == colnames(dataBenchTRG)[j],]
         ixx <- match(targetList, xx$TARGET_NAME)
         ixx[is.na(ixx)] <- 21
         items2 <- ixx
         ## generate the ranked lists, order, and then rank
         #drugERankSimil1[i,j] <-  drugERankSimil1[j,i] <- cor(items1,items2,method="spearman")
         drugERankSimil2[i,j] <-  drugERankSimil2[j,i] <- cor(items1,items2,method= "kendall")
    }
  }
}


if (benchNam=="ctrp") {
    save(drugERankSimil2, file="drugERankSimil-CTRPV2-kendall.Rdata")
}
else {
    
    #save(drugERankSimil1, file="drugERankSimil-nci60-spearman.Rdata")
    save(drugERankSimil2, file="drugERankSimil-nci60-kendall.Rdata")

}


