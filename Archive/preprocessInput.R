###############################################################################################################
## Function reads in the input datasets, preprocess and finds the set of drugs intersecting between the two sets
##
## input: 
##     dname: name ("character") of the first dataset, e.g., "nci60" , "ctrpv2"   
##     d2: currently, must be set to "lincs"  
## output: 
##     a list ("character") or "dataframe" (if dname==nci60) containing the intersection of the two datasets 			
##
## 
###############################################################################################################



preprocessInput <- function(dname , d2="lincs") {

   ## Read LINCS metadata file
   if (d2 == "lincs") {  
      lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
   } 

   if (dname == "ctrpv2") {
      ## Read ctrpv2 metadata file
      ctrpv2 <- read.csv("Data/CTRPv2_drugtarget.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
      ## capitalize + remove badchars and intersect LINCS chemicals from LINCS with ctrpv2 compounds
      lincs$pert_iname <- toupper(lincs$pert_iname)
      lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
      ctrpv2$compound_name <- toupper(ctrpv2$compound_name)
      ctrpv2$compound_name <- gsub(badchars,"",ctrpv2$compound_name)
      intrsct <- intersect(lincs$pert_iname, ctrpv2$compound_name) # 239 drugs/481 total
    } else {
        if (dname == "nci60") { 
             ## load data including smiles from NCI60
             load("Data/nci60_with_smiles.RData") # 49938 drugs x 5 column descriptions
             ## intersect LINCS cpds with NCI60 compounds
             inters <- intersect(lincs[,"pubchem_cid"], nci60_with_smiles[,"PubChem CID"]) # 360, all unique by default 
             ## keep only the intersecting drugs
             lincsInters <-lincs[match(inters,lincs$pubchem_cid),]
             nciInters <- nci60_with_smiles[match(inters,nci60_with_smiles$`PubChem CID`),]
             ## ordering of pubchem cids for easy matching
             lincsInters <- lincsInters[order(lincsInters[,"pubchem_cid"]),]
             nciInters <- nciInters[order(as.numeric(as.character(nciInters[,4]))),]
             ## keep unique
             lincsUnique <- lincsInters[!duplicated(lincsInters[,"pubchem_cid"]),]
             nciUnique <- nciInters[!duplicated(nciInters[,"PubChem CID"]),]
             all(nciUnique$`PubChem CID` == lincsUnique$pubchem_cid) #Sanity check
             intrsct <- data.frame("NSC"=nciUnique[,1], lincsUnique[,c("canonical_smiles","pert_iname","pert_id","pubchem_cid")])
             intrsct$pert_iname <- gsub(badchars, "",  intrsct[,"pert_iname"])
             intrsct$pert_iname <- toupper(intrsct$pert_iname)
             intrsct <- intrsct[!duplicated(intrsct[,"pert_iname"]),] #353 unique drugs
        }
 }
  
  
  return(intrsct) 

}