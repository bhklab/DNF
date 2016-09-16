
##NOTE this code resumes with iskar's method after mean-centering...
## does the average within cell lines , skips the detection
## calculates the pairwise dips scores within cell lines...

## NOTE: not all drugs have experiments in all cell lines ...
## so, their sigs are gonna be NA...
## if a score is calculated, it must be ignored!


# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")

rm(list=ls())

##source("connectivityScore.R")

args <- commandArgs(trailingOnly = TRUE)
print(length(args))
args

## n shows the number of up and down regulated genes. can be 30, 50 or 70
n <- as.integer(args[1])

## 
library("PharmacoGx")
library(PharmacoGx)
library(parallel)
## 


badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## get ALL the required drug names (only benchmarks) and the perturbation ids (note: a subset of DNF_pert_all.RData)
## 
## drug names of ctrp
load("Data/atc-ctrpv-new.RData") # 51
ctrpATC <- dataBench3 
load("Data/targetSET_CTRPV2.RData") # 141
ctrpTRG <- dataBench
## drug names of nci
load("Data/atc-NCI60-new.RData")
nci60ATC <- dataBench3  # 72
load("Data/target-uniProt-nci60.RData")
nci60TRG <- x   # 86

allDrugNames <- union(rownames(ctrpATC), rownames(ctrpTRG))
allDrugNames <- union(allDrugNames, rownames(nci60ATC)) 
allDrugNames <- union(allDrugNames, rownames(nci60TRG)) # 218

lincs <- read.csv("LINCS.csv", stringsAsFactors = FALSE)
lincs$pert_iname <- toupper(lincs$pert_iname)
lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)

length(intersect(lincs$pert_iname, allDrugNames))
ind <- which(lincs$pert_iname %in% allDrugNames) ## length 315 
all(lincs$pert_iname[ind] %in% allDrugNames)
lincsSub <- lincs[ind,] ## subset lincs according to the drugNames

load("Data/DNF_pert_all.RData") ## has only one pert id for each drug
indx <- which(lincsSub$pert_id %in% names(l1000.drug.signatures[1,,1])) ## get pert ids
length(intersect(lincsSub$pert_id[indx] , names(l1000.drug.signatures[1,,1])))

target <- cbind(name= lincsSub[indx,]$pert_iname, id =lincsSub[indx,]$pert_id)
target <- as.data.frame(target)
length(unique(target$name))

## now subset the perturbation object
xx <- which(names(l1000.drug.signatures[1,,1]) %in% target[,2])
pertLincsSub <- l1000.drug.signatures[,xx,]
dim(pertLincsSub)
#[1] 978 219   4

## loop must be n*n -n --> 1,2 not same as 2,1
#n=30
m = dim(target)[1]
scoreM <- matrix(nrow = m, ncol = m)

##get thedrug names
ii <- match(names(pertLincsSub[1,,1]), target$id)
all(target[ii,2] == names(pertLincsSub[1,,1]))
##[1] TRUE
rownames(scoreM) <- colnames(scoreM) <- target[ii,1]

for (i in 1:m) {
  ## generate the signatures 
  ## 0- sort according to the tstat for each drug
  pp <- pertLincsSub[order(pertLincsSub[,i,"tstat"]),i,"tstat"]
  
  ## 1- pick top and bottom n genes ...
  qSigU <- rep(1,n)
  names(qSigU) <- names(pp[(length(pp)-n+1):length(pp)])
  qSigD <- rep(-1,n)
  names(qSigD) <- names(pp[1:n])   
  ## 2- generate the drug signature 
  qSig <- c(qSigU, qSigD)
  
  cl <- makeCluster(16)
  res <- parApply(pertLincsSub[,,c("tstat", "fdr")], 2, function(x, QSIG){
    return(PharmacoGx::connectivityScore(x=x,
                                         y=QSIG,
                                         method="gsea", nperm=101))
  }, cl = cl, QSIG=qSig)
  
  stopCluster(cl)
  ## store the score...
  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  if (all(rownames(res)==colnames(scoreM[i,]))) { scoreM[i,] <- res[,1] } else print(paste(i,"error!", sep=""))
}


save(res, i, scoreM, file=paste("scoreM",n ,i,"-all-iorioPGX-101.RData",sep=""))


## post-processing on the drug similarity matrix to make it symmetric
average <- scoreM
for (i in 1:dim(scoreM)[1])
for(j in 1:dim(scoreM)[1]) {
    average[i,j] <- average[j,i] <- mean(c(scoreM[i,j], scoreM[j,i]))
}

save(average, file="Data/averageIorioPGX-all.RData")



print("the end!")


