# The purpose of this script is to add drug structure data to an existing
# structure network. Once this works properly, we can think about making the script
# accept JSON data later on.

rm(list=ls())

library(rcdk)
source("RCode/flexible_layers/constStructureLayerFlexible.R")
source("RCode/flexible_layers/structureDataFlexible.R")

# Load lincs metadata and clean the pert_inames in it
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

# Load sensitivity data in order to subset the lincs metadata so that we are dealing with a network
# of reasonable size.
sens.data <- readRDS("Data/combined_sensitivity/combined_sens_iname_replaced.RData")
drugs <- rownames(sens.data)

# Subset the lincs metadata based on the drugs that are available in the combined
# sensitivity dataset
lincs.meta.subsetted <- lincs.meta[lincs.meta$pert_iname %in% drugs, ]
lincs.meta.subsetted <- lincs.meta.subsetted[!duplicated(lincs.meta.subsetted$pert_iname),]
common.drugs <- lincs.meta.subsetted$pert_iname

# Create fingerprints and tanimoto similarity matrix for the drugs that 
# are common between the sensitvity data and lincs metadata
fingerprints <- StructureDataFlexible(lincs.meta.subsetted)
structure.data.old <- ConstStructureLayerFlexible(fingerprints)

# Take a random compound from the lincs metadata that we didn't previously compute a fingerprint for,
# compute a fingerprint for it, 
new.compound <- lincs.meta[!(lincs.meta$pert_iname %in% drugs),][1, c("pert_iname", "canonical_smiles")]
new.fingerprint <- rcdk::parse.smiles(new.compound$canonical_smiles)
new.fingerprint <- rcdk::get.fingerprint(new.fingerprint[[1]])

# Add the new fingerprint to the list of fingerprints and compute a new
# tanimoto similarity matrix. 
fingerprints[[new.compound$pert_iname]] <- new.fingerprint
structure.data.new <- ConstStructureLayerFlexible(fingerprints)

# As a sanity check, if we consider only the drugs that are common between the old tanimoto
# similarity matrix and the new one, we should find that they are identical.
all(structure.data.old[common.drugs, common.drugs] == structure.data.new[common.drugs, common.drugs])
