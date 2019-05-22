# This script creates a barplot showing disagreements in terms of targets 
# between the various drug target datasets being used. Each bar corresponds
# to the number of drug that have targets which disagree between a pair of datasets.

rm(list=ls())
library(jsonlite)
library(ggplot2)
library(UniProt.ws)
library(extrafont)
library("org.Hs.eg.db")
library(reactome.db)
source("RCode/exploratory_analysis/targetDiffHelpers.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
up <- UniProt.ws(taxId=9606)

cdrugs <- read.csv("Data/drugsToAnnotate.csv", stringsAsFactors = FALSE)[,1]
cdrugs <- toupper(cdrugs)
cdrugs <- gsub(badchars, "", cdrugs)

### Get drug targets from CTRPv2
ctrp.drug.targets <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
ctrp.drug.targets$compound_name <- toupper(ctrp.drug.targets$compound_name)
ctrp.drug.targets$compound_name <- gsub(badchars,"",ctrp.drug.targets$compound_name)

ctrp.drug.targets <- ctrp.drug.targets[ctrp.drug.targets$compound_name %in% cdrugs,,drop=F] 
ctrp.drug.targets <- ctrp.drug.targets[,c(1,2)] 
colnames(ctrp.drug.targets)[1:2] <- c("MOLECULE_NAME","TARGET_NAME")
# split to get unique targets 
drug.targets <- strsplit(ctrp.drug.targets$TARGET_NAME, split = ";")
# assign a unique target corresponding to each drug (a target can have multiple assigned drugs and
# a drug can be found in different target categories)
ctrp.drug.targets <- data.frame(MOLECULE_NAME = rep(ctrp.drug.targets$MOLECULE_NAME, sapply(drug.targets, length)), TARGET_NAME = unlist(drug.targets),
                                stringsAsFactors = FALSE)
# single target categorty 
ctrp.drug.targets <- ctrp.drug.targets[,c("MOLECULE_NAME","TARGET_NAME")]


### Get drug targets from the downloaded CHEMBL file
chembl.drug.targets <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
# remove badchar from Chembl and common drug file + capitalize 
chembl.drug.targets[,"MOLECULE_NAME"] <- gsub(badchars, "",  chembl.drug.targets[,"MOLECULE_NAME"])
chembl.drug.targets[,"MOLECULE_NAME"] <- toupper(chembl.drug.targets[,"MOLECULE_NAME"])

chembl.drug.targets <- chembl.drug.targets[chembl.drug.targets$MOLECULE_NAME %in% cdrugs,]
chembl.drug.targets$gene.symbol <- NA

chembl.drug.targets$gene.symbol <- GetCHEMBLGeneSymbols(chembl.drug.targets$TARGET_CHEMBL_ID)
chembl.drug.targets <- chembl.drug.targets[, c("MOLECULE_NAME", "gene.symbol")]
colnames(chembl.drug.targets) <- c("MOLECULE_NAME", "TARGET_NAME")

### Get drug targets from Drug Bank
dbank.drug.targets <- read.csv("./Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
dbank.drug.targets$gene.symbol <- NA

# Replace the UniProt.ID with gene symbol
dbank.drug.targets$gene.symbol <- GetDrugBankGeneSymbols(up, dbank.drug.targets$UniProt.ID)
dbank.drug.targets$gene.symbol <- as.character(dbank.drug.targets$gene.symbol)

uniprot.targets <- dbank.drug.targets[,c("Name", "gene.symbol")]
colnames(uniprot.targets) <-  c("MOLECULE_NAME","TARGET_NAME")
# remove badchar from Chembl and common drug file + capitalize 
uniprot.targets[,1] <- toupper(uniprot.targets[,1])
uniprot.targets[,1] <- gsub(badchars, "",  uniprot.targets[,1])
uniprot.targets <- uniprot.targets[uniprot.targets[,1] %in% cdrugs,]

drug.bank.targets <- uniprot.targets

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
                                               sapply(clue.targets, length)), TARGET_NAME = unlist(clue.targets),
                           stringsAsFactors = FALSE)

chembl.drug.targets <- chembl.drug.targets[!(chembl.drug.targets$TARGET_NAME == ""), ]
clue.targets <- clue.targets[!(clue.targets$TARGET_NAME == ""), ]
ctrp.drug.targets <- ctrp.drug.targets[!(ctrp.drug.targets$TARGET_NAME == ""), ]
drug.bank.targets <- drug.bank.targets[!(drug.bank.targets$TARGET_NAME == ""), ]

chembl.drug.targets <- chembl.drug.targets[!is.na(chembl.drug.targets$TARGET_NAME), ]
clue.targets <- clue.targets[!is.na(clue.targets$TARGET_NAME), ]
ctrp.drug.targets <- ctrp.drug.targets[!is.na(ctrp.drug.targets$TARGET_NAME), ]
drug.bank.targets <- drug.bank.targets[!is.na(drug.bank.targets$TARGET_NAME), ]


res <- list(chembl.targets=chembl.drug.targets, clue.io.targets=clue.targets, ctrp.targets=ctrp.drug.targets,
            drug.bank.targets=drug.bank.targets)

saveRDS(res, "../../pharmacodb/targetsForDesiredDrugs.rds")