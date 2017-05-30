###############################################################################################################
## Function reads in the benchmark sets and filter them out according to the intersection of the input datasets
## 
## input: 
##     cdrugs: a vector of common drugs between the input datasets
## output: 
##        a drug x drug adjacency/similarity matrix
##
## 
###############################################################################################################


DrugTargetBenchFlexible <- function(cdrugs, gmt_file_name="", use.ctrpv2=FALSE,
                                  use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE) {
    drug.targets = data.frame(MOLECULE_NAME=character(0), TARGET_NAME=character(0))
    
    if (use.ctrpv2) {
        # DRUG TARGETS BENCHMARKING FROM CTRPV2
        ctrp.drug.targs <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
        ctrp.drug.targs$compound_name <- toupper(ctrp.drug.targs$compound_name)
        ctrp.drug.targs$compound_name <- gsub(badchars,"",ctrp.drug.targs$compound_name)
        ## 
        ctrp.drug.targs <- ctrp.drug.targs[ctrp.drug.targs$compound_name %in% cdrugs,,drop=F] 
        ctrp.drug.targs <- ctrp.drug.targs[,c(1,2)] 
        colnames(ctrp.drug.targs)[1:2] <- c("MOLECULE_NAME","TARGET_NAME")
        ## split to get unique targets 
        drug.targets <- strsplit(ctrp.drug.targs$TARGET_NAME, split = ";")
        ## assign a unique target corresponding to each drug (a target can have multiple assigned drugs and
        ## a drug can be found in different target categories)
        drug.targets <- data.frame(MOLECULE_NAME = rep(ctrp.drug.targs$MOLECULE_NAME, sapply(drug.targets, length)), TARGET_NAME = unlist(drug.targets))
        ## single target categorty 
        drug.targets <- drug.targets[,c("MOLECULE_NAME","TARGET_NAME")]
    }
    
    if (use.chembl) {
        ## read CHEMBL drug file downloaded from Chembl website
        chembl.drug.targs <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
        chembl.drug.targs <- chembl.drug.targs[,c("MOLECULE_NAME", "TARGET_NAME")]
        ## remove badchar from Chembl and common drug file + capitalize 
        chembl.drug.targs[,1] <- gsub(badchars, "",  chembl.drug.targs[,1])
        chembl.drug.targs[,1] <- toupper(chembl.drug.targs[,1])
        
        chembl.drug.targs <- chembl.drug.targs[chembl.drug.targs$MOLECULE_NAME %in% cdrugs,]
        chembl.drug.targs <- chembl.drug.targs[chembl.drug.targs[,1] %in% cdrugs,]
        
        drug.targets <- rbind.data.frame(drug.targets, chembl.drug.targs, stringsAsFactors = FALSE)
    }
    
    if (use.dbank) {
        ## Drug bank targets
        dbank.drug.targs <- read.csv("./Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
        uniprot.targs <- dbank.drug.targs[,c("Name", "UniProt.Name")]
        colnames(uniprot.targs) <-  c("MOLECULE_NAME","TARGET_NAME")
        ## remove badchar from Chembl and common drug file + capitalize 
        uniprot.targs[,1] <- toupper(uniprot.targs[,1])
        uniprot.targs[,1] <- gsub(badchars, "",  uniprot.targs[,1])
        uniprot.targs <- uniprot.targs[uniprot.targs[,1] %in% cdrugs,]
        
        drug.targets <- rbind.data.frame(drug.targets, uniprot.targs, stringsAsFactors = FALSE)
    }
    
    if (use.clue) {
        #####
        ### Get drug targets from clue.io
        clue.io.targets <- read.delim("Data/repurposing_drugs_20170327.txt", stringsAsFactors = FALSE)
        
        # Ignore first several rows which are just meta about the dataset
        clue.io.targets <- clue.io.targets[which(clue.io.targets$X.1 == 'target'):nrow(clue.io.targets), ]
        # By this point the first row has the relevant column names, so we set the column names based 
        # on this row
        colnames(clue.io.targets) <- clue.io.targets[1, ]
        
        # Get rid of that first row, since now the colnames() attribute has been set to it
        clue.io.targets <- clue.io.targets[-1, ]
        # Take only the two relevant columns of the drug name and target
        clue.io.targets <- clue.io.targets[, c("pert_iname", "target")]
        # Standardize column names like all the other column names in the benchmarks 
        colnames(clue.io.targets) <- c("MOLECULE_NAME", "TARGET_NAME")
        
        clue.io.targets$MOLECULE_NAME <- toupper(clue.io.targets$MOLECULE_NAME)
        clue.io.targets$MOLECULE_NAME <- gsub(badchars, "", clue.io.targets$MOLECULE_NAME)
        
        clue.io.targets <- clue.io.targets[clue.io.targets$MOLECULE_NAME %in% cdrugs, ]
        
        clue.targets <- strsplit(clue.io.targets$TARGET_NAME, split="|", fixed=TRUE)
        clue.targets <- data.frame(MOLECULE_NAME = rep(clue.io.targets$MOLECULE_NAME,
                                                       sapply(clue.targets, length)), TARGET_NAME = unlist(clue.targets))
        drug.targets <- rbind.data.frame(drug.targets, clue.targets, stringsAsFactors = FALSE)
        ###
        #####    
    }
    
    if (use.dtc) {
        dtc.targs <- read.csv("Data/dtcTargets.csv", stringsAsFactors = FALSE)
        colnames(dtc.targs) <-  c("MOLECULE_NAME","TARGET_NAME")
        
        drug.targetsAdditional <- dtc.targs[dtc.targs[,1] %in% cdrugs,]
        
        drug.targets <- rbind.data.frame(drug.targets, drug.targetsAdditional, stringsAsFactors = FALSE)
    }
    
    #Number of drugs with TARGETS: length(unique(chembl_Targets.common.all$MOLECULE_NAME))
    drug.targets <- unique(drug.targets)
    drug.targets <- drug.targets[!is.na(drug.targets[,2]),]
    
    ## create a gmt file from Target codes, each Target class contains a number of drugs
    list.of.targs <- list()
    for(targName in unique(drug.targets$TARGET_NAME)){
        targGmt <- with(drug.targets, MOLECULE_NAME[TARGET_NAME==targName]) 
        list.of.targs[[targName]] <- targGmt
    }
    ## filter out the targets common between more than 2 drugs
    common.targs <- sapply(list.of.targs, function(x) length(x) >= 2)
    GMT_TARG<- list.of.targs[common.targs]
    save(GMT_TARG, file = paste(getwd(), "/Output/", gmt_file_name, sep=""))
    
    ## Build an adjacency matrix target x drugs and keep only filtered targets in GMT file
    ## data frame with 0/1
    data.wide <- reshape2::dcast(drug.targets, TARGET_NAME ~ MOLECULE_NAME)
    data.wide <- data.wide[data.wide$TARGET_NAME %in% names(GMT_TARG),, drop=F]
    rownames(data.wide) <- data.wide[,1]
    data.wide <- data.wide[,-1]
    ## keep only filtered drugs in GMT file
    ## only drugs with known target
    data.wide <- data.wide[, colnames(data.wide) %in% unique(as.character(unlist(GMT_TARG))), drop=F] 
    data.bench <- apply(data.wide, 2, function(x) ifelse(!is.na(x),1,0))
    ## create a drug x drug adjacency matrix
    data.bench <- as.matrix(proxy::simil(t(data.bench), method="Jaccard"))
    data.bench <- apply(data.bench, 1, function(x) ifelse(x>0 | is.na(x),1,0))
    
    return(data.bench)
}


