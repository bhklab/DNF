compareCTRPv2 <- function(badchars) {
    ctrpv2_csv <- read.csv("Data/CTRPv2_drugtarget.csv",stringsAsFactors = FALSE) # 481 drugs x 28 descriptions
    ctrpv2_csv_drug_names <- toupper(ctrpv2_csv$compound_name)
    ctrpv2_csv_drug_names <- gsub(badchars,"",ctrpv2_csv_drug_names)

    load("./PSets/CTRPv2.RData")
    ctrpv2_pset_drug_names <- drugNames(CTRPv2)
    ctrpv2_pset_drug_names <- toupper(ctrpv2_pset_drug_names)
    ctrpv2_pset_drug_names <- gsub(badchars, "", ctrpv2_pset_drug_names)
    
    cat("Number of drugs in csv:")
    cat(paste(length(ctrpv2_csv_drug_names), "\n"))
    
    cat("Number of drugs in Pset:")
    cat(paste(length(ctrpv2_pset_drug_names), "\n"))
    
    cat("Number of drugs in intersection:")
    cat(paste(length(intersect(ctrpv2_csv_drug_names, ctrpv2_pset_drug_names)), "\n"))
    
    cat("Drugs in csv version that are not in Pset version\n")
    sort(ctrpv2_csv_drug_names[-which(ctrpv2_csv_drug_names %in% ctrpv2_pset_drug_names)])
}


compareIntersctionWithLINCS <- function(badchars) {
    lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
    ## capitalize + remove badchars from lincs
    lincs$pert_iname <- toupper(lincs$pert_iname)
    lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
    
    ctrpv2_csv <- read.csv("Data/CTRPv2_drugtarget.csv",stringsAsFactors = FALSE) # 481 drugs x 28 descriptions
    ctrpv2_csv_drug_names <- toupper(ctrpv2_csv$compound_name)
    ctrpv2_csv_drug_names <- gsub(badchars,"",ctrpv2_csv_drug_names)
    
    load("./PSets/CTRPv2.RData")
    ctrpv2_pset_drug_names <- drugNames(CTRPv2)
    ctrpv2_pset_drug_names <- toupper(ctrpv2_pset_drug_names)
    ctrpv2_pset_drug_names <- gsub(badchars, "", ctrpv2_pset_drug_names)
    
    intrsctLincsCSV <- intersect(lincs$pert_iname, ctrpv2_csv_drug_names)
    lincsInters <- lincs[lincs$pert_iname %in% intrsctLincsCSV,,drop=F]
    print(length(lincsInters))
    lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F]
    
    intrsctLincsCSV <- intersect(lincs$pert_iname, ctrpv2_pset_drug_names)
    lincsInters <- lincs[lincs$pert_iname %in% intrsctLincsCSV,,drop=F]
    length(lincsInters)
    lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F]
    
    
}
