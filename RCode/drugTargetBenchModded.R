###############################################################################################################
## Function reads in the benchmark sets and filter them out according to the intersection of the input datasets
## 
## input: 
##     benchname: name ("character") of the drug-target benchmark, e.g., "ctrpv", "chembl", or "stitch"
##     cdrugs: a vector of common drugs between the input datasets
## output: 
##        a drug x drug adjacency/similarity matrix
##
## 
###############################################################################################################


drugTargetBenchModded <- function(benchname, cdrugs, gmt_file_name="") {
    if (benchname == "ctrpv") {
        # DRUG TARGETS BENCHMARKING FROM CTRPV2
        ctrpDrTargs <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
        ctrpDrTargs$compound_name <- toupper(ctrpDrTargs$compound_name)
        ctrpDrTargs$compound_name <- gsub(badchars,"",ctrpDrTargs$compound_name)
        ## 
        ctrpDrTargs <- ctrpDrTargs[ctrpDrTargs$compound_name %in% cdrugs,,drop=F] 
        ctrpDrTargs <- ctrpDrTargs[,c(1,2)] 
        colnames(ctrpDrTargs)[1:2] <- c("MOLECULE_NAME","TARGET_NAME")
        ## split to get unique targets 
        drgTargets <- strsplit(ctrpDrTargs$TARGET_NAME, split = ";")
        ## assign a unique target corresponding to each drug (a target can have multiple assigned drugs and
        ## a drug can be found in different target categories)
        drgTargets <- data.frame(MOLECULE_NAME = rep(ctrpDrTargs$MOLECULE_NAME, sapply(drgTargets, length)), TARGET_NAME = unlist(drgTargets))
        ## single target categorty 
        drgTargets <- drgTargets[,c("MOLECULE_NAME","TARGET_NAME")]
        
    }
    
    ## read CHEMBL drug file downloaded from Chembl website
    chemblDrTargs <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
    chemblDrTargs <- chemblDrTargs[,c("MOLECULE_NAME", "TARGET_NAME")]
    ## remove badchar from Chembl and common drug file + capitalize 
    chemblDrTargs[,1] <- gsub(badchars, "",  chemblDrTargs[,1])
    chemblDrTargs[,1] <- toupper(chemblDrTargs[,1])
    
    chemblDrTargs <- chemblDrTargs[chemblDrTargs$MOLECULE_NAME %in% cdrugs,]
    chemblDrTargs <- chemblDrTargs[chemblDrTargs[,1] %in% cdrugs,]
    
    drgTargets <- rbind.data.frame(drgTargets, chemblDrTargs, stringsAsFactors = FALSE)
    
    ## read CHEMBL drug file downloaded from Chembl website
    chemblDrTargs <- read.csv("./Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
    uniprotTargs <- chemblDrTargs[,c("Name", "UniProt.Name")]
    colnames(uniprotTargs) <-  c("MOLECULE_NAME","TARGET_NAME")
    ## remove badchar from Chembl and common drug file + capitalize 
    uniprotTargs[,1] <- gsub(badchars, "",  uniprotTargs[,1])
    uniprotTargs[,1] <- toupper(uniprotTargs[,1])
    uniprotTargs <- uniprotTargs[uniprotTargs[,1] %in% cdrugs,]
    
    drgTargets <- rbind.data.frame(drgTargets, uniprotTargs, stringsAsFactors = FALSE)
    
    if (TRUE) {
        dtcTargs <- read.csv("Data/dtcTargets.csv", stringsAsFactors = FALSE)
        colnames(dtcTargs) <-  c("MOLECULE_NAME","TARGET_NAME")
        
        drgTargetsAdditional <- dtcTargs[dtcTargs[,1] %in% cdrugs,]
        
        drgTargets <- rbind.data.frame(drgTargets, drgTargetsAdditional, stringsAsFactors = FALSE)
    }
    
    #Number of drugs with TARGETS: length(unique(chembl_Targets.common.all$MOLECULE_NAME))
    drgTargets <- unique(drgTargets)
    drgTargets <- drgTargets[!is.na(drgTargets[,2]),]
    
    ## create a gmt file from Target codes, each Target class contains a number of drugs
    listofTargs <- list()
    for(targName in unique(drgTargets$TARGET_NAME)){
        targGmt <- with(drgTargets, MOLECULE_NAME[TARGET_NAME==targName]) 
        listofTargs[[targName]] <- targGmt
    }
    ## filter out the targets common between more than 2 drugs
    commonTargs <- sapply(listofTargs, function(x) length(x) >= 2)
    GMT_TARG<- listofTargs[commonTargs]
    save(GMT_TARG, file = paste(getwd(), "/Output/", gmt_file_name, sep=""))
    
    ## Build an adjacency matrix target x drugs and keep only filtered targets in GMT file
    ## data frame with 0/1
    dataWide <- reshape2::dcast(drgTargets, TARGET_NAME ~ MOLECULE_NAME)
    dataWide <- dataWide[dataWide$TARGET_NAME %in% names(GMT_TARG),, drop=F]
    rownames(dataWide) <- dataWide[,1]
    dataWide <- dataWide[,-1]
    ## keep only filtered drugs in GMT file
    ## only drugs with known target
    dataWide <- dataWide[, colnames(dataWide) %in% unique(as.character(unlist(GMT_TARG))), drop=F] 
    dataBench <- apply(dataWide, 2, function(x) ifelse(!is.na(x),1,0))
    ## create a drug x drug adjacency matrix
    dataBench <- as.matrix(proxy::simil(t(dataBench), method="Jaccard"))
    dataBench <- apply(dataBench, 1, function(x) ifelse(x>0 | is.na(x),1,0))
    
    return(dataBench)
}


