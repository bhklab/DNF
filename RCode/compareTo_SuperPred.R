

## scripts that parses the superpred results and calculates dru-drug similarity matrix for ctrv2 and nci60
## input to the script is the benchmark name "ctrp" or "nci60"

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
print(length(args))
args


benchNam <- as.character(args[1])


library(xlsx)
x <- read.xlsx("Data/superPredResultsFinal.xlsx", 1)
x <- x[!duplicated(x$drugName),] ## 85 unique drug names 

if (benchNam=="ctrp") {
   load("Data/atc-ctrpv-new.RData")
   dataBechATC <- dataBench3 ## 51
} else {
  load("Data/atc-NCI60-new.RData")
  dataBechATC <- dataBench3 ##  72 for nci60
}

xx <- x[x$drugName %in% colnames(dataBechATC),] ##51 --> 72

## remove the NAs
xx <- xx[!is.na(xx$prediction),] ## 37 --> 62
xx <- xx[xx$prediction!="",] ## 37


uniAtc <- data.frame(SYNONYMS = "xx", ATC_CODE = "xxx", Score = 0.00) #11544 rows
uniAtc$SYNONYMS <- as.character(uniAtc$SYNONYMS)
uniAtc$ATC_CODE <- as.character(uniAtc$ATC_CODE)
uniAtc$Score <- as.numeric(uniAtc$Score)

for (i in 1:dim(xx)[1]) {
  p <- unlist(strsplit(as.character(xx[i,2]), ','))
  #score <- p[length(p)]
  for (k in 1:length(p)) {
     p2 <- unlist(strsplit(as.character(p[k]), '[()]'))
     p2 <- p2[(p2 !="")] ## length(p2) -1 is the no. of codes... the last item is the score..
     for (j in 1:((length(p2))-1)) {
        uniAtc <-  rbind(uniAtc, c(as.character(xx$drugName[i]), p2[j], as.numeric(p2[length(p2)]))) #11544 rows
     }
  }
}
uniAtc <- uniAtc[-1,]

colnames(uniAtc) <- c("MOLECULE_NAME", "ATC_CODE", "SCORE")
length(unique(uniAtc$ATC_CODE))
#[1] 107 --> 155

## 
listofCodes <- list()
for(atcCode in unique(uniAtc$ATC_CODE)) {
  ind <- which(uniAtc$ATC_CODE == atcCode)
  atcGmt <- uniAtc$MOLECULE_NAME[ind]
  #atcGmt <- with(uniAtc, uniAtc$MOLECULE_NAME[uniAtc$ATC_CODE==atcCode]) 
  listofCodes[[atcCode]] <- atcGmt
}


## keep categories with length >= 2 drugs AND create final GMT file
commonChemb <- sapply(listofCodes, function(x) length(x) >= 2)
GMT_ATC <- listofCodes[commonChemb] ## only 21 out of 107
length(GMT_ATC) ## 21 --> 62


################################### 
## data frame with 0/1
dataWide <- reshape2::dcast(uniAtc, ATC_CODE ~ MOLECULE_NAME)
dataWide <- dataWide[dataWide$ATC_CODE %in% names(GMT_ATC),,drop=F] ## 21 X 34 --> 63 X 62
rownames(dataWide) <- dataWide[,1]
dataWide <- dataWide[,-1]

## only drugs with known target
dataWide <- dataWide[, colnames(dataWide) %in% unique(as.character(unlist(GMT_ATC))),drop=F] 

# ## for nci60 
ss2 <- setdiff(colnames(dataBench3), colnames(dataWide)) ## 18 --> 12

## NOW generate the similarity matrix from the dataWide matrix ...
SuperPredsimil <- matrix(nrow = length(colnames(dataBench3)), ncol=length(colnames(dataBench3)), 0)
rownames(SuperPredsimil) <- colnames(SuperPredsimil) <- c(colnames(dataWide), ss2) 

## for each row, look at the ones ...
for (i in 1:dim(dataWide)[1]) {
  ## 
  code <- rownames(dataWide)[i]
  ind <- which(dataWide[i,]!=0)  
  for (j in 1:length(ind)) {
    for(k in ((j+1):length(ind))) {
       if (k <= length(ind)) {
        ii <- intersect(which(uniAtc$MOLECULE_NAME==colnames(dataWide)[ind[j]]) , (which(uniAtc$ATC_CODE==code)))
        try(if(length(ii)!=1) stop("ii issue"))
        pi <- as.numeric(uniAtc$SCORE[ii])
        jj <- intersect(which(uniAtc$MOLECULE_NAME==colnames(dataWide)[ind[k]]) , (which(uniAtc$ATC_CODE==code)))
        try(if(length(jj)!=1) stop("jj issue"))
        pj <- as.numeric(uniAtc$SCORE[jj])
        
        distSuprPred <- 1-pi+ 1-pj+ abs(pi-pj) + 0.001
        
        ## 1/dist --> similarity
        cc <- which(colnames(SuperPredsimil)==colnames(dataWide)[ind[j]])
        dd <- which(colnames(SuperPredsimil)==colnames(dataWide)[ind[k]])
        
        SuperPredsimil[cc,dd] <- SuperPredsimil[dd,cc] <- 1/distSuprPred 
       }
    }
  }
    }
 
if (benchNam=="ctrp") { 
    save(SuperPredsimil, file="Data/SuperPredsimil-CTRPV2.Rdata")
} else {
    save(SuperPredsimil, file="Data/SuperPredsimil-NCI60.Rdata")
}