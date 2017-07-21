rm(list=ls())

library(PharmacoGx)

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

drugs <- unique(unlist(sapply(datasets, 
              function(x) {
                  drug.names <- drugNames(x)
                  combos.removed <- drug.names[!grepl(":", drug.names)]
                  return (gsub(badchars, "",toupper(combos.removed))) 
              })))

mapping.table <- read.csv("Data/chantal_chembl_ids.csv")

mapping.table$drug.names <- toupper(mapping.table$drug.names)
mapping.table$drug.names <- gsub(badchars, "", mapping.table$drug.names)

dtc.targets <- read.csv("Data/DTC_data.csv", stringsAsFactors=F)
dtc.targets$compound_name <- toupper(dtc.targets$compound_name)
dtc.targets$compound_name <- gsub(badchars, "", dtc.targets$compound_name)

dtc.targets.subsetted <- dtc.targets[dtc.targets$standard_inchi_key %in% mapping.table$inchikey, ]
dtc.targets.subsetted <- rbind(dtc.targets.subsetted, 
                               dtc.targets[dtc.targets$compound_name %in% mapping.table$drug.names, ])
dtc.targets.subsetted <- rbind(dtc.targets.subsetted, 
                               dtc.targets[dtc.targets$compound_id %in% mapping.table$chembl_ids, ])


dtc.targets.subsetted <- dtc.targets.subsetted[dtc.targets.subsetted$gene_names != "", ]
dtc.targets.subsetted <- dtc.targets.subsetted[dtc.targets.subsetted$compound_name != "", ]

dtc.targets.subsetted <- dtc.targets.subsetted[,c("compound_name", "gene_names")]
dtc.targets.subsetted <- dtc.targets.subsetted[!duplicated(dtc.targets.subsetted), ]
colnames(dtc.targets.subsetted) <- c("MOLECULE_NAME", "TARGET_NAME")

dtc.target.names <- strsplit(dtc.targets.subsetted$TARGET_NAME, split=",", fixed=TRUE)
dtc.targets.subsetted <- data.frame(MOLECULE_NAME = rep(dtc.targets.subsetted$MOLECULE_NAME,
                                               sapply(dtc.target.names, length)), TARGET_NAME = unlist(dtc.target.names),
                                    stringsAsFactors = FALSE)

dtc.target.names <- strsplit(dtc.targets.subsetted$TARGET_NAME, split="-", fixed=TRUE)
dtc.targets.subsetted <- data.frame(MOLECULE_NAME = rep(dtc.targets.subsetted$MOLECULE_NAME,
                                                        sapply(dtc.target.names, length)), TARGET_NAME = unlist(dtc.target.names),
                                    stringsAsFactors = FALSE)

write.csv(dtc.targets.subsetted, "Data/dtcTargets.csv", row.names = F)
