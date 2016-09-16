###############################################################################################################
## Function reads in the benchmark set and filter them out according to the intersection of the input datasets
## 
## input: 
##     benchname: name ("character") of the ATC benchmark, currently must be set to "chembl"
##     cdrugs: a vector of common drugs between the input datasets
## output: 
##        a drug x drug adjacency/similarity matrix
##
## 
###############################################################################################################





ATCBench <- function(dname, cdrugs) {
  
 if (dname=="chembl") {
     chembl_ATC <- read.delim("Data/chembl_drugs-15_3_18_59.txt", stringsAsFactor=F, na.strings=c("", "NA")) #10506 drugs!
     ## remove unused strings from synonyms or drug names such as ()
     sep = " \\("
     chembl_ATC$SYNONYMS <- toupper(unlist(lapply(chembl_ATC$SYNONYMS, function(x) strsplit(x, sep)[[1]][1])))
     ## split to get unique ATCs
     uniAtc <- strsplit(chembl_ATC$ATC_CODE, split = "; ")
     ## assign a unique ATC corresponding to each drug
     uniAtc <- data.frame(SYNONYMS = rep(chembl_ATC$SYNONYMS, sapply(uniAtc, length)), ATC_CODE = unlist(uniAtc)) #11544 rows
     ## trim on ATC and keep level 4 only
     uniAtc$ATC_CODE <- sapply(uniAtc[,2], function(x) strtrim(x, 5))
     ## remove badchar from Chembl and common drug file + capitalize (nci60 and LINCS)
     uniAtc$SYNONYMS <- gsub(badchars, "",  uniAtc$SYNONYMS)
  
     ## intersect (drugs with ATCs) with (CHEMBL 354 drugs) & keep those drugs only
     ## build a target GMT file includind drugs in common with all datasets and having known ATC associations
     intersChemblCdrugs <- uniAtc[uniAtc$SYNONYMS %in% cdrugs$pert_iname,]  ## NOTE: in nehme's code: CommonDrugs   i get 310 2
     intersChemblCdrugs <- unique(intersChemblCdrugs)  ## iget 284
  
     ## filtering chembl files before creating the GMT file, get known ATC, exclude where ATCs are NA
     uniAtc <- intersChemblCdrugs[!is.na(intersChemblCdrugs[,"ATC_CODE"]),] # 250 rows    ## 243 2

     ## create a gmt file from ATC codes, each ATC class contains a number of drugs
     listofCodes <- list()
     for(atcCode in unique(uniAtc$ATC_CODE)) {
       atcGmt <- with(uniAtc, SYNONYMS[ATC_CODE==atcCode]) 
       listofCodes[[atcCode]] <- atcGmt
     }
  
     ## keep categories with length >= 2 drugs AND create final GMT file
     commonChemb <- sapply(listofCodes, function(x) length(x) >= 2)
     GMT_ATC <- listofCodes[commonChemb]
     length(commonChemb) ## 144
  
     ## data frame with 0/1
     dataWide <- reshape2::dcast(uniAtc, ATC_CODE ~ SYNONYMS)
     dataWide <- dataWide[dataWide$ATC_CODE %in% names(GMT_ATC),,drop=F]
     rownames(dataWide) <- dataWide[,1]
     dataWide <- dataWide[,-1]
  
     ## only drugs with known target
     dataWide <- dataWide[, colnames(dataWide) %in% unique(as.character(unlist(GMT_ATC))),drop=F] 
     dataBench <- apply(dataWide, 2, function(x) ifelse(!is.na(x),1,0))
     save(dataBench, file=paste(getwd(), "/Output/", "atc_bench_chembl.RData", sep=""))

     dataBench <- as.matrix(proxy::simil(t(dataBench), method="Jaccard"))
     dataBench <- apply(dataBench, 1, function(x) ifelse(x >0 | is.na(x), 1, 0))
     save(dataBench, file=paste(getwd(), "/Output/", "drug_bench_atc_binary.RData", sep=""))
  }
  
  return(dataBench)
  
}