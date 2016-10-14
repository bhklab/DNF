

rm(list=ls())

## this script integrates a partial rnaking method to compute drug-drug similarities for drugE-Rank results ...
## 

benchNam <- "nci60"

if (benchNam=="ctrp") {
  load("targetSET_CTRPV2.RData") 
  dataBenchTRG <- dataBench ## 141
} else {
  load("target-uniProt-nci60.RData")
  dataBenchTRG <- x ##  86 for nci60
}


# badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
# 
# ## readin the uniprot dataset to have the mapping between drug names 
# chemblDrTargs <- read.csv("Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
# #drgTargets <- chemblDrTargs[,c("Name", "UniProt.Name", "UniProt.ID")]
# #colnames(drgTargets) <-  c("MOLECULE_NAME","TARGET_NAME", "ID")
# ## remove badchar from Chembl and common drug file + capitalize 
# chemblDrTargs[,"Name"] <- gsub(badchars, "",  chemblDrTargs[,"Name"])
# chemblDrTargs[,"Name"] <- toupper(chemblDrTargs[,"Name"])

## 
#load("drugERankResults_ctrp.RData") ## -> drgTargs
load("drugERankResults_nci60.RData")
targetList <- unique(drgTargs$TARGET_NAME)
m <-length(unique(drgTargs$TARGET_NAME)) ## 373 unique target ids
drgTargs$SCORE <- 21-(as.numeric(as.character(drgTargs$SCORE))*20)
#rownames(drgTargs) <- drgTargs$MOLECULE_NAME
#score <- c(20:1)/20

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
         drugERankSimil1[i,j] <-  drugERankSimil1[j,i] <- cor(items1,items2,method="spearman")
         drugERankSimil2[i,j] <-  drugERankSimil2[j,i] <- cor(items1,items2,method= "kendall")
    }
  }
}

save(drugERankSimil1, file="drugERankSimil-nci60-spearman.Rdata")
save(drugERankSimil2, file="drugERankSimil-nci60-kendall.Rdata")



