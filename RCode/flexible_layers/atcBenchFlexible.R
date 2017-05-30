###############################################################################################################
## Function reads in the benchmark set and filter them out according to the intersection of the input datasets
## 
## input: 
##     benchname: name ("character") of the ATC benchmark, currently must be set to "chembl"
##     common.drugs: a vector of common drugs between the input datasets
## output: 
##        a drug x drug adjacency/similarity matrix
##
## 
###############################################################################################################
ATCBenchFlexible <- function(dname, common.drugs) {
    
    if (dname=="chembl") {
        chembl.atc <- read.delim("Data/chembl_drugs-15_3_18_59.txt", stringsAsFactor=F, na.strings=c("", "NA")) #10506 drugs!
        ## remove unused strings from synonyms or drug names such as ()
        sep = " \\("
        chembl.atc$SYNONYMS <- toupper(unlist(lapply(chembl.atc$SYNONYMS, function(x) strsplit(x, sep)[[1]][1])))
        ## split to get unique ATCs
        uni.atc <- strsplit(chembl.atc$ATC_CODE, split = "; ")
        ## assign a unique ATC corresponding to each drug
        uni.atc <- data.frame(SYNONYMS = rep(chembl.atc$SYNONYMS, sapply(uni.atc, length)), ATC_CODE = unlist(uni.atc)) #11544 rows
        ## trim on ATC and keep level 4 only
        uni.atc$ATC_CODE <- sapply(uni.atc[,2], function(x) strtrim(x, 5))
        ## remove badchar from Chembl and common drug file + capitalize (nci60 and LINCS)
        uni.atc$SYNONYMS <- gsub(badchars, "",  uni.atc$SYNONYMS)
        
        ## intersect (drugs with ATCs) with (CHEMBL 354 drugs) & keep those drugs only
        ## build a target GMT file includind drugs in common with all datasets and having known ATC associations
        intersc.chembl.cdrugs <- uni.atc[uni.atc$SYNONYMS %in% common.drugs,]  ## NOTE: in nehme's code: CommonDrugs   i get 310 2
        intersc.chembl.cdrugs <- unique(intersc.chembl.cdrugs)  ## iget 284
        
        ## filtering chembl files before creating the GMT file, get known ATC, exclude where ATCs are NA
        uni.atc <- intersc.chembl.cdrugs[!is.na(intersc.chembl.cdrugs[,"ATC_CODE"]),] # 250 rows    ## 243 2
        
        ## create a gmt file from ATC codes, each ATC class contains a number of drugs
        list.of.codes <- list()
        for(atcCode in unique(uni.atc$ATC_CODE)) {
            atcGmt <- with(uni.atc, SYNONYMS[ATC_CODE==atcCode]) 
            list.of.codes[[atcCode]] <- atcGmt
        }
        
        ## keep categories with length >= 2 drugs AND create final GMT file
        common.chemb <- sapply(list.of.codes, function(x) length(x) >= 2)
        GMT_ATC <- list.of.codes[common.chemb]
        length(common.chemb) ## 144
        
        ## data frame with 0/1
        data.wide <- reshape2::dcast(uni.atc, ATC_CODE ~ SYNONYMS)
        data.wide <- data.wide[data.wide$ATC_CODE %in% names(GMT_ATC),,drop=F]
        rownames(data.wide) <- data.wide[,1]
        data.wide <- data.wide[,-1]
        
        ## only drugs with known target
        data.wide <- data.wide[, colnames(data.wide) %in% unique(as.character(unlist(GMT_ATC))),drop=F] 
        data.bench <- apply(data.wide, 2, function(x) ifelse(!is.na(x),1,0))
        save(data.bench, file=paste(getwd(), "/Output/", "atc_bench_chembl.RData", sep=""))
        
        data.bench <- as.matrix(proxy::simil(t(data.bench), method="Jaccard"))
        data.bench <- apply(data.bench, 1, function(x) ifelse(x >0 | is.na(x), 1, 0))
        save(data.bench, file=paste(getwd(), "/Output/", "drug_bench_atc_binary.RData", sep=""))
    } else if (dname=="chembl-new") {
        
        chembl.atc <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #10506 drugs!
        ## remove unused strings from synonyms or drug names such as ()
        sep = " \\("
        chembl.atc$MOLECULE_NAME <- toupper(unlist(lapply(chembl.atc$MOLECULE_NAME, function(x) strsplit(x, sep)[[1]][1])))
        ## split to get unique ATCs
        uni.atc <- strsplit(chembl.atc$ATC_CODE, split = "; ")
        ## assign a unique ATC corresponding to each drug
        uni.atc <- data.frame(MOLECULE_NAME = rep(chembl.atc$MOLECULE_NAME, sapply(uni.atc, length)), ATC_CODE = unlist(uni.atc)) #11544 rows
        ## trim on ATC and keep level 4 only
        uni.atc$ATC_CODE <- sapply(uni.atc[,2], function(x) strtrim(x, 5))
        ## remove badchar from Chembl and common drug file + capitalize (nci60 and LINCS)
        uni.atc$MOLECULE_NAME <- gsub(badchars, "",  uni.atc$MOLECULE_NAME)
        
        ## intersect (drugs with ATCs) with (CHEMBL 354 drugs) & keep those drugs only
        ## build a target GMT file includind drugs in common with all datasets and having known ATC associations
        intersc.chembl.cdrugs <- uni.atc[uni.atc$MOLECULE_NAME %in% common.drugs,]  ## NOTE: in nehme's code: CommonDrugs   i get 310 2
        intersc.chembl.cdrugs <- unique(intersc.chembl.cdrugs)  ## iget 284
        
        ## filtering chembl files before creating the GMT file, get known ATC, exclude where ATCs are NA
        uni.atc <- intersc.chembl.cdrugs[!is.na(intersc.chembl.cdrugs[,"ATC_CODE"]),] # 250 rows    ## 243 2  --> 81 2
        
        ## create a gmt file from ATC codes, each ATC class contains a number of drugs
        list.of.codes <- list()
        for(atcCode in unique(uni.atc$ATC_CODE)) {
            atcGmt <- with(uni.atc, MOLECULE_NAME[ATC_CODE==atcCode]) 
            list.of.codes[[atcCode]] <- atcGmt
        }
        
        ## keep categories with length >= 2 drugs AND create final GMT file
        common.chemb <- sapply(list.of.codes, function(x) length(x) >= 2)
        GMT_ATC <- list.of.codes[common.chemb]
        length(common.chemb) ## 144, actually looks like it's 41
        
        ## data frame with 0/1
        data.wide <- reshape2::dcast(uni.atc, ATC_CODE ~ MOLECULE_NAME)
        data.wide <- data.wide[data.wide$ATC_CODE %in% names(GMT_ATC),,drop=F]
        rownames(data.wide) <- data.wide[,1]
        data.wide <- data.wide[,-1]
        
        ## only drugs with known target
        data.wide <- data.wide[, colnames(data.wide) %in% unique(as.character(unlist(GMT_ATC))),drop=F] 
        data.bench <- apply(data.wide, 2, function(x) ifelse(!is.na(x),1,0))
        save(data.bench, file=paste(getwd(), "/Output/", "atc_bench_chembl_NEW.RData", sep=""))
        
        data.bench <- as.matrix(proxy::simil(t(data.bench), method="Jaccard"))
        data.bench <- apply(data.bench, 1, function(x) ifelse(x >0 | is.na(x), 1, 0))
        save(data.bench, file=paste(getwd(), "/Output/", "drug_bench_atc_binary_NEW.RData", sep=""))
    }
    
    
    
    
    return(data.bench)
    
}