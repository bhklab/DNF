# This script creates a barplot showing disagreements in terms of targets 
# between the various drug target datasets being used. Each bar corresponds
# to the number of drug that have targets which disagree between a pair of datasets.

rm(list=ls())
library(ggplot2)
library(UniProt.ws)
library(extrafont)
library("org.Hs.eg.db")
library(reactome.db)
source("RCode/exploratory_analysis/targetDiffHelpers.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
up <- UniProt.ws(taxId=9606)
sensitivity.file.name <- "Data/combined_sens_adjusted_diag_inamedatasets_with.RData"
sensitivity.data <- readRDS(sensitivity.file.name)
cdrugs <- rownames(sensitivity.data)

### Get drug targets from CTRPv2
ctrp.drug.targs <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
ctrp.drug.targs$compound_name <- toupper(ctrp.drug.targs$compound_name)
ctrp.drug.targs$compound_name <- gsub(badchars,"",ctrp.drug.targs$compound_name)

ctrp.drug.targs <- ctrp.drug.targs[ctrp.drug.targs$compound_name %in% cdrugs,,drop=F] 
ctrp.drug.targs <- ctrp.drug.targs[,c(1,2)] 
colnames(ctrp.drug.targs)[1:2] <- c("MOLECULE_NAME","TARGET_NAME")
# split to get unique targets 
drug.targets <- strsplit(ctrp.drug.targs$TARGET_NAME, split = ";")
# assign a unique target corresponding to each drug (a target can have multiple assigned drugs and
# a drug can be found in different target categories)
ctrp.drug.targs <- data.frame(MOLECULE_NAME = rep(ctrp.drug.targs$MOLECULE_NAME, sapply(drug.targets, length)), TARGET_NAME = unlist(drug.targets),
                              stringsAsFactors = FALSE)
# single target categorty 
ctrp.drug.targs <- ctrp.drug.targs[,c("MOLECULE_NAME","TARGET_NAME")]


### Get drug targets from the downloaded CHEMBL file
chembl.drug.targs <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
chembl.drug.targs <- chembl.drug.targs[,c("MOLECULE_NAME", "TARGET_NAME")]
# remove badchar from Chembl and common drug file + capitalize 
chembl.drug.targs[,1] <- gsub(badchars, "",  chembl.drug.targs[,1])
chembl.drug.targs[,1] <- toupper(chembl.drug.targs[,1])

chembl.drug.targs <- chembl.drug.targs[chembl.drug.targs$MOLECULE_NAME %in% cdrugs,]
chembl.drug.targs <- chembl.drug.targs[chembl.drug.targs[,1] %in% cdrugs,]

### Get drug targets from Drug Bank
dbank.drug.targs <- read.csv("./Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
dbank.drug.targs$gene.symbol <- NA

# Replace the UniProt.ID with gene symbol
res <- select(up, keys = dbank.drug.targs$UniProt.ID, "GENES", "UNIPROTKB")    
splitted = strsplit(res$GENES, split=" ")
first.gene.names <- sapply(splitted, function(x) {x[1]})
dbank.drug.targs$gene.symbol <- first.gene.names
dbank.drug.targs$gene.symbol <- as.character(dbank.drug.targs$gene.symbol)

uniprot.targs <- dbank.drug.targs[,c("Name", "gene.symbol")]
colnames(uniprot.targs) <-  c("MOLECULE_NAME","TARGET_NAME")
# remove badchar from Chembl and common drug file + capitalize 
uniprot.targs[,1] <- toupper(uniprot.targs[,1])
uniprot.targs[,1] <- gsub(badchars, "",  uniprot.targs[,1])
uniprot.targs <- uniprot.targs[uniprot.targs[,1] %in% cdrugs,]

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

### Get targets from Drug Target Commons
dtc.targs <- read.csv("Data/dtcTargets.csv", stringsAsFactors = FALSE)
colnames(dtc.targs) <-  c("MOLECULE_NAME","TARGET_NAME")

dtc.targs <- dtc.targs[dtc.targs[,1] %in% cdrugs,]

### Determine the number of drugs between all possible dataset pairs
### which have disagreeing targets. I.e. dataset 1 says the target for 
### some Drug X is JAK1, while dataset 2 says that the target for 
### Drug X is MAP3K4.
all.datasets <- list(ctrp.drug.targs=ctrp.drug.targs,
                     clue.targets=clue.targets, uniprot.targs=uniprot.targs,
                     dtc.targs=dtc.targs)
discrepancy.counts <- list()
xx <- as.list(reactomeEXTID2PATHID)

for (i in 1:length(all.datasets)) {
    for (j in 1:length(all.datasets)) {
        if (i != j & j > i) {
            d1 = all.datasets[[i]]
            d2 = all.datasets[[j]]
            
            common.drugs <- intersect(d1$MOLECULE_NAME, d2$MOLECULE_NAME)
            d1 <- d1[d1$MOLECULE_NAME %in% common.drugs, , drop=FALSE]
            d2 <- d2[d2$MOLECULE_NAME %in% common.drugs, , drop=FALSE]
            
            pair.name <- paste(names(all.datasets)[i], " - ", names(all.datasets)[j], sep="")
            
            discrepancy.counts[[pair.name]] <- c()
            
            discrepancy.counts <- GetTargetDiscrepanciesLenient(d1, d2, common.drugs, discrepancy.counts, pair.name)
        }
    }    
}

### Create a dataframe that is required by ggplot to create the plot
actual.counts <- sapply(discrepancy.counts, length)
plot.data <- data.frame(names(discrepancy.counts), actual.counts)
colnames(plot.data) <- c("pair.name", "pair.counts")

### Create the bar plot
g <- ggplot(plot.data, aes(pair.name, pair.counts)) 
g <- g + geom_bar(stat = "identity") + ggtitle("Drug Target Dataset # of Target Disagreements") 
g <- g + labs(x="\nDataset Pair", y="Number of Disagreeing Drugs\n")
g + theme(axis.text=element_text(size=11),
          axis.title=element_text(size=16, family="Gadugi"),
          title=element_text(size=20, family="Century Gothic"))
