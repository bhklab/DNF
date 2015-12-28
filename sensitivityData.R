###############################################################################################################
## Function reads in the sensitivity data, processes/filters according to the intersection subset provided in the preprocess stage
##
## input: 
##     dname: name ("character") of the dataset, e.g., "nci60" , "ctrpv2"   
##     intersc: a vector of "character" (for "ctrpv2") or dataframe (for "nci60") containting intersection info.
## output: 
##     sensitivity data ("dataframe") 	
##
## 
###############################################################################################################



sensitivityData <- function(dname, intersc) {

   if (dname == "ctrpv2") {
        ## Read drug sensitivity data from ctrpv2
        ctrpv2Sens <- read.csv("Data/ctrpv2_celldrug.csv",stringsAsFactors = FALSE) # 
        ctrpv2Sens <- ctrpv2Sens[!duplicated(ctrpv2Sens[,1]),,drop=F]
        rownames(ctrpv2Sens) <- ctrpv2Sens[,1]
        ctrpv2Sens <- ctrpv2Sens[,-1]
        ## Remove the columns with all NAs
        ctrpv2Sens <- ctrpv2Sens[colSums(is.na(ctrpv2Sens)) < nrow(ctrpv2Sens)]
        ## remove the rows with all NAs 
        ctrpv2Sens <- ctrpv2Sens[rowSums(is.na(ctrpv2Sens)) < ncol(ctrpv2Sens),]
        ## capitalize + remove badchars from colnames
        colnames(ctrpv2Sens) <- toupper(colnames(ctrpv2Sens))
        colnames(ctrpv2Sens) <- gsub(badchars,"",colnames(ctrpv2Sens))
        ## keep only drugs intersecting with LINCS
        ctrpv2Sens <- ctrpv2Sens[,colnames(ctrpv2Sens) %in% intersc,drop=F]
        sens <- ctrpv2Sens
        
    } else if (dname == "nci60") {
        ## Read drug sensitivity data from nci60 & reduce NCI60 to the drugs common between L1000/NCI60
        load("Data/NCI60.RData") #from PharmacoGx
        NCIdrugs <- PharmacoGx::drugInfo(NCI60)[PharmacoGx::drugInfo(NCI60)[,2] %in% intersc[,1],] #353 drugs
        NCIcells <- PharmacoGx::cellNames(NCI60) #60 celllines
        all(sort(as.numeric(NCIdrugs$`NSC #`)) == sort(intersc$NSC)) #Sanity check
        subx <- PharmacoGx::subsetTo(NCI60, cells = NCIcells, drugs = NCIdrugs[,1]) #Subset the PharmacoSet to 353 compounds
        ## Summarize sensitivity. Get recomputed AUC using PharmacoGx pkg
        NCI60Auc <- PharmacoGx::summarizeSensitivityPhenotype(subx,sensitivity.measure="auc_recomputed", summaryStat="median") 
        ## Remove the columns with all NAs
        NCI60Auc <- NCI60Auc[,colSums(is.na(NCI60Auc)) < nrow(NCI60Auc)]
        ## remove the rows with all NAs 
        NCI60Auc <- NCI60Auc[rowSums(is.na(NCI60Auc)) < ncol(NCI60Auc),]
        colnames(NCI60Auc) <- gsub("drugid_","",colnames(NCI60Auc))
        colnames(NCI60Auc) <- intersc[match(colnames(NCI60Auc),intersc$NSC),]$pert_iname
        class(NCI60Auc) <- "numeric"
        dim(NCI60Auc) 
        sens <- NCI60Auc
    }
  
  
    return(sens)
  
  }