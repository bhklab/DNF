


rm(list=ls())


args <- commandArgs(trailingOnly = TRUE)
print(length(args))
args

## n shows the number of up and down regulated genes. can be 30, 50 or 70
n <- as.integer(args[1])
n <-30


## 
library("PharmacoGx")
load("Data/DNF-meanCentered.RData")
library(PharmacoGx)
library(parallel)
## 


badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## get the required drug names (benchmarks) and the perturbation ids 
## drug names of ctrpv ATC

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


## get the perturbation ids ...
lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
## capitalize + remove badchars from lincs
lincs$pert_iname <- toupper(lincs$pert_iname)
lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
lincsAllDrugNames <- lincs[lincs$pert_iname %in% allDrugNames,,drop=F] ## 315 pert ids
length(unique(lincsAllDrugNames$pert_iname))
#[1] 218
uniqLincsAllDrugNames <- unique(lincsAllDrugNames$pert_iname)


sampleinfo <- phenoInfo(DNF_l1000_compounds_subsetted, "rna")  ## sampleinfo has 77 unique cells
dim(sampleinfo)
length(unique(sampleinfo$cell_id))
uniqCellIds <- unique(sampleinfo$cell_id)

## data structure to store the signatures 
aveSigCells <- array(dim=c(length(uniqCellIds), length(uniqLincsAllDrugNames), 
            dim(pertData_dnf_meanCentered)[1])
            , dimnames = list(unique(sampleinfo$cell_id), uniqLincsAllDrugNames, rownames(pertData_dnf_meanCentered))) ## celllines X drugs X genes 

#minUniqCellIds <- vector(length = length(uniqCellIds))
for (i in 1:length(uniqLincsAllDrugNames)) {
  
  ## for each drug 
  ## get the pert_ids 
  pertIds <- lincsAllDrugNames[lincsAllDrugNames$pert_iname == uniqLincsAllDrugNames[i], ]$pert_id
  
  ## for each cellline average the replicates ...
  ## subset again according to cell lines ...
  
  exps <- sampleinfo[sampleinfo$pert_id %in% pertIds,] ## 72 exps   
  #length(unique(exps$cell_id))
  drgUniqCellIds <- unique(exps$cell_id)
  
  for (j in 1:length(drgUniqCellIds)) {
      ## average over the same cell ids
      exps2 <- exps[exps$cell_id==drgUniqCellIds[j],]
      ## get the row names of exps2 and match with the pert data
      xx <- pertData_dnf_meanCentered[ ,colnames(pertData_dnf_meanCentered) %in% rownames(exps2)]
      if (dim(exps2)[1]==1) {aveSigCells[drgUniqCellIds[j], uniqLincsAllDrugNames[i], ] <- xx}
      else {
        xx <- rowMeans(xx)
        aveSigCells[drgUniqCellIds[j], uniqLincsAllDrugNames[i], ] <- xx
      }
  }
  ## detection call ranking --> skipped 
}
save(aveSigCells, file="Output/aveSigCells-30.RData")

## pass the array to calculate DIPS within cell lines ...
## then average over all the cell lines ... put zero for the missing values (drug pairs)! 

scores <- array(dim=c(length(uniqCellIds), length(uniqLincsAllDrugNames), length(uniqLincsAllDrugNames)),
                dimnames = list(uniqCellIds, uniqLincsAllDrugNames, uniqLincsAllDrugNames))

## subset for each cell line
for (k in 1:length(uniqCellIds) ) {
  ## rank the genes for all drugs in the cell line
  
  for (l in 1:length(uniqLincsAllDrugNames)) { ## rank individual drug signatures only if there's no NA
     if (sum(is.na(aveSigCells[uniqCellIds[k], l ,]))==0) {
      aveSigCells[uniqCellIds[k], l ,] <- rank(aveSigCells[uniqCellIds[k], l ,])
      ###### on Sep 13 NOTE: had to flip the direction becasue of weird connectivity scores....
      aveSigCells[uniqCellIds[k], l ,] <- dim(aveSigCells)[3] -  aveSigCells[uniqCellIds[k], l ,] ## reverse the rank: we want descending order of the average difference
     } ## othewise, keep na s
  }
  #if (!(all(aveSigCells[uniqCellIds[k],,]>=0))) stop("ERROR! in ranking!");
  ## do dips for drug pairs ...
  sigs <- t(aveSigCells[uniqCellIds[k], ,]) ##[1] 978 218
  for (z in 1:length(uniqLincsAllDrugNames)) {
    ## generate the signatures 
    ## 0- sort according to the tstat for each drug
    ##pp <- pertLincsSub[order(pertLincsSub[,i,"tstat"]),i,"tstat"]
    pp <- sigs[,uniqLincsAllDrugNames[z]]
    if (sum(is.na(pp))==0) {
      names(pp) <- rownames(sigs)
      pp <- pp[order(pp, decreasing= FALSE)] ## note: rank is increasing
    
      ## 1- pick top and bottom n genes ...
      qSigU <- rep(-1,n)
      names(qSigU) <- names(pp[1:n])
      qSigD <- rep(1,n)
      names(qSigD) <- names(pp[(length(pp)-n+1):length(pp)])
      ## 2- generate the drug signature 
      qSig <- c(qSigU, qSigD)
    
      cl <- makeCluster(4)
      res <- parApply(sigs, 2, function(x, QSIG){
        
          if (sum(is.na(x)==0)) {
            source('./RCode/connectivityScore2.R') #####    <----- NOTE... may need to be changed
            source('./RCode/runGSA2.R')
            source('./RCode/checkLoadArg2.R')
            source('./RCode/GSCstatBatch.R')
            source('./RCode/calcGeneSetStat.R')
            source('./RCode/GSCsignificanceBatch.R')
            source('./RCode/GSCstatGenePerm.R')
            source('./RCode/pvalFromFractionGenePerm.R')
            source('./RCode/fdrGSEA.R')
            source('./RCode/combineTest.R')  
            
         return(connectivityScore2(x=x,
                                           y=QSIG,
                                           method="gsea", nperm=1))
            } else return(-2)
      }, cl = cl, QSIG=qSig)
    
      stopCluster(cl)
      ## store the score...
      if ( (class(res)=="list") & (length(res)==length(uniqLincsAllDrugNames)) ) {
         if (all(names(res) == colnames(scores[uniqCellIds[k],z,]))) {
           for (h in 1:length(uniqLincsAllDrugNames)) {
               scores[uniqCellIds[k],z,h] <- res[[h]][1]  ## store the score ... 
           }
           }else { stop("error in matching the dims!")}
    }
  }
 } 
}

save(scores, file="Ouput/scores-30.RData")
save(aveSigCells, file="Ouput/aveSigCellsNEW-30.RData")

save(list = ls(all.names = TRUE), file = "Output/iskarRun-n=30.RData")
######################################

### 1) convert the -2 scores to NA
scoresN <-  array(dim=c(length(uniqCellIds), length(uniqLincsAllDrugNames), length(uniqLincsAllDrugNames)),
                  dimnames = list(uniqCellIds, uniqLincsAllDrugNames, uniqLincsAllDrugNames))

for (k in 1:length(uniqCellIds) ) {
  scoresN[k,,] <- apply(scores[k,,], c(1,2), function(x) {if (!is.na(x) & (x == -2)) {x <- NA} else {x <- x} } )
}
### 2) make the score matrices symetric by averaging ij and ji entries 
for (k in 1:length(uniqCellIds) ) {
  for (i in 1:dim(scoresN)[2]) {
  for(j in 1:dim(scoresN)[3]) 
    scoresN[k,i,j] <- scoresN[k,j,i] <- mean(c(scoresN[k,i,j], scoresN[k,j,i]))
  }
}

### 3) average the scores for each drug pair over all the cell lines and compute one single drug-drug score matrix

finalIskarScore <- matrix(NA, nrow = dim(scoresN)[2], ncol= dim(scoresN)[3])
rownames(finalIskarScore) <- colnames(finalIskarScore) <- rownames(scoresN[1,,])

for (i in 1:dim(scoresN)[2]) {
  for(j in i:(dim(scoresN)[3])) {
      if (j<dim(scoresN)[3]) {
        ## for each drug ij average over all cell lines
        ## first remove the NA cell lines ...
        sc <- scoresN[1:length(uniqCellIds),i,j]
        sc <- sc[!is.na(sc)] ## there might be some drug pairs not existing is any common cell line ... see below for an example
        if (length(sc)!=0) {
            finalIskarScore[i,j] <- finalIskarScore[j,i] <- mean(sc)
        }
      }
  }
}

### convert the Nan to zero! 
finalIskarScore <-apply(finalIskarScore, c(1,2), function(x) {if (is.na(x)) {x<-0} else {x <- x }} )
save(finalIskarScore, file="Output/averageIskarFinal.RData")

######
### example of a drug pair not found in a single common cell line --> the score is NA 
uniqLincsAllDrugNames[1]
# "CYTARABINE"
uniqLincsAllDrugNames[92]
# "APICIDIN"

pertIds <- lincsAllDrugNames[lincsAllDrugNames$pert_iname == uniqLincsAllDrugNames[1], ]$pert_id
exps <- sampleinfo[sampleinfo$pert_id %in% pertIds,] ## 72 exps
drgUniqCellIds <- unique(exps$cell_id)
drgUniqCellIds
#[1] "A375"  "A549"  "HA1E"  "HT29"  "MCF7"  "PC3"   "VCAP"  "HS27A" "SKM1"  "U937" 
pertIds <- lincsAllDrugNames[lincsAllDrugNames$pert_iname == uniqLincsAllDrugNames[92], ]$pert_id
exps <- sampleinfo[sampleinfo$pert_id %in% pertIds,] ## 72 exps   
drgUniqCellIds <- unique(exps$cell_id)
drgUniqCellIds
#[1] "FIBRNPC" "NEU"     "NEU.KCL" "NPC"   



