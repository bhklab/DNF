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
      ## capitalize + remove badchars from lincs
      lincs$pert_iname <- toupper(lincs$pert_iname)
      lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
   } else { stop("ERROR!")}

   if (dname == "ctrpv2") {
      ## Read ctrpv2 metadata file
      ctrpv2 <- read.csv("Data/CTRPv2_drugtarget.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
      ctrpv2$compound_name <- toupper(ctrpv2$compound_name)
      ctrpv2$compound_name <- gsub(badchars,"",ctrpv2$compound_name)
      #intersect LINCS chemicals from LINCS with ctrpv2 compounds
      intrsctLincsCtrpv2 <- intersect(lincs$pert_iname, ctrpv2$compound_name) # 239 drugs/481 total
      lincsInters <- lincs[lincs$pert_iname %in% intrsctLincsCtrpv2,,drop=F]
      lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F] ## 239 drugs
      intrsct <- lincsInters
    } else {
        if (dname == "nci60") { 
          nci <- read.csv("Data/DTP_NCI60_ZSCORE.csv",stringsAsFactors = FALSE, na.strings="na") # 20861 drugs x 66 descriptions
          nci <- nci[!nci$Drug.name=="-",,drop=F]
          nci$Drug.name <- toupper(nci$Drug.name)
          nci$Drug.name  <- gsub(badchars,"", nci$Drug.name)
          
          ## intersect LINCS cpds with NCI60 compounds by names or smiles to maximize number of drugs
          intrscLincsNciDrugNamLst <- intersect(lincs$pert_iname, nci$Drug.name) # 192, all unique by default
          intrscLincsNciSmilLst <- intersect(lincs$canonical_smiles, nci$SMILES..d.) # 106
          
          ## keep only the intersecting drugs by names and match the names
          lincsIntrsNciDrugNam <-lincs[match(intrscLincsNciDrugNamLst,lincs$pert_iname),] ## 192 X 28
          nciIntrsLincsDrugNam <- nci[match(intrscLincsNciDrugNamLst,nci$Drug.name),] ## 192 X 66
          
          ## keep only the intersecting drugs by smiles and match smiles
          lincsIntrsNciSmil <- lincs[match(intrscLincsNciSmilLst, lincs$canonical_smiles),] ## 106 X 28 
          nciIntrsLincsSmil <- nci[match(intrscLincsNciSmilLst,nci$SMILES..d.),] ## 106 X 66
          
          ### match between nci and lincs and assign same name
          nciIntrsLincsSmil <- nciIntrsLincsSmil[match(lincsIntrsNciSmil$canonical_smiles, nciIntrsLincsSmil$SMILES..d.), ] # 106 X 66
          nciIntrsLincsSmil$Drug.name <- lincsIntrsNciSmil$pert_iname
          
          ######## rbind lincs files/rbind nci60 files which keeps drugs in common = 238
          lincsboth <- rbind(lincsIntrsNciDrugNam, lincsIntrsNciSmil)
          lincsboth <- lincsboth[!duplicated(lincsboth$pert_iname),,drop=F]
          
          nciboth <- rbind(nciIntrsLincsDrugNam, nciIntrsLincsSmil)
          nciboth <- nciboth[!duplicated(nciboth$Drug.name),,drop=F]
          
          intrsct <- list(lincsboth=lincsboth, nciboth=nciboth)
        }
 }
  
  
  return(intrsct) 

}