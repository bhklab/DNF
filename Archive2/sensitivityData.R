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
        ## NOTE jan 20
        ctrpv2Sens <- ctrpv2Sens[,colnames(ctrpv2Sens) %in% intersc$pert_iname,drop=F]
        sens <- ctrpv2Sens
        
    } else if (dname == "nci60") {
        NCI60Auc <- intersc[,c(2,7:ncol(intersc)),drop=F]
        ## 1st column as the rownames
        rownames(NCI60Auc) <- NCI60Auc[,1]
        ## drop the 1st column
        NCI60Auc <- NCI60Auc[,-1]
        ## transpose the matrix
        NCI60Auc <- t(NCI60Auc)
        dim(NCI60Auc)
        sens <- NCI60Auc
    }
  
  
    return(sens)
  
  }