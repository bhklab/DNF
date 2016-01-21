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


drugTargetBench <- function(benchname, cdrugs) {

  
  if (benchname == "ctrpv") {
       # DRUG TARGETS BENCHMARKING FROM CTRPV2
       ctrpDrTargs <- read.csv("Data/ctrpv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
       ctrpDrTargs$compound_name <- toupper(ctrpDrTargs$compound_name)
       ctrpDrTargs$compound_name <- gsub(badchars,"",ctrpDrTargs$compound_name)
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

  } else if (benchname == "chembl") {
       ## read CHEMBL drug file downloaded from Chembl website
       chemblDrTargs <- read.delim("Data/chembl_drugtargets-15_3_46_00.txt", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
       drgTargets <- chemblDrTargs[,c("MOLECULE_NAME", "TARGET_NAME")]
       ## remove badchar from Chembl and common drug file + capitalize 
       drgTargets[,1] <- gsub(badchars, "",  drgTargets[,1])
       drgTargets[,1] <- toupper(drgTargets[,1])
       drgTargets <- drgTargets[drgTargets[,1] %in% cdrugs,]
       

  } else if (benchname == "stitch") {
       ## DRUG TARGETS BENCHMARKING
       stitchDrTargs <- read.delim("Data/stitch.tsv", stringsAsFactor=F, na.strings=c("", "NA")) ## read protein-chemical from stitch embl http://stitch.embl.de/
       stitchDrTargs$chemical <- gsub("CID1","",stitchDrTargs$chemical) ## the CID1 is a standard id for merged isomers and stereoismers
       stitchDrTargs$chemical <- gsub("CID0","",stitchDrTargs$chemical) ## the CID0 is a standard id for unique cpds
       stitchDrTargs$chemical <- gsub("(^|[^0-9])0+", "\\1", stitchDrTargs$chemical, perl = TRUE) ## remove 0s to keep the pubchem CID
       stitchDrTargs <- stitchDrTargs[stitchDrTargs$experimental!=0 & stitchDrTargs$database!=0 ,,drop=F] ### select all interactions
       stitchDrTargs$protein <- gsub("9606.", "", stitchDrTargs$protein) ## remove the 9606 (homo sapiens taxon and keep ENSEMBL PROT ID)
       interscChem <- intersect(stitchDrTargs$chemical, cdrugs$pubchem_cid) # intersect pubchem id between the common 345 drugs in LINCS/NCI60 and STITCH
       stitchDrTargs <- stitchDrTargs[stitchDrTargs$chemical %in% interscChem,,drop=F] ## Keep data from STITCH (drugs in common)
       drugsCommonStitch <- cdrugs[cdrugs$pubchem_cid %in% interscChem,,drop=F] ## Keep data from LINCS/NCI60 ((drugs in common))
       colnames(stitchDrTargs)[1] <- "pubchem_cid" ### change 1st column in STITCH to pubchem_cid so we can merge between stitch and LINCS/NCI60 pheno data file
       ## merge STITCH database and drugsCommonStitch
       stitchDrTargs <- merge(stitchDrTargs, drugsCommonStitch, by="pubchem_cid") ## merge data
       ## Take only the 2 columns ("protein","pert_iname")
       stitchDrTargs <- stitchDrTargs[,c("protein","pert_iname")]
       stitchDrTargs <- stitchDrTargs[,c(2,1)]
       colnames(stitchDrTargs)[1:2] <- c("MOLECULE_NAME","TARGET_NAME") 
       ## single target categorty 
       drgTargets <- stitchDrTargs[,c("MOLECULE_NAME","TARGET_NAME")] ### rename to starg, maybe to keep stitch_targets in case of error!
       #dim(drgTargets)
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
    # (ie, remove drugs with only one target)
    commonTargs <- sapply(listofTargs, function(x) length(x) >= 2)
    GMT_TARG<- listofTargs[commonTargs]
    save(GMT_TARG, file = paste(getwd(), "/Output/", "gmt_targ_", benchname , ".RData", sep=""))
   
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
  

