.libPaths('/mnt/work1/users/bhklab/Rlib/')
library(PharmacoGx)
load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/L1000_compounds.RData")

common.drugs <- readRDS("common_drugs.RData")
cell.lines <- readRDS("cell_lines.RData")

myx <- phenoInfo(L1000_compounds, mDataType="rna")$duration %in% c(0,6)
L1000_compounds@molecularProfiles$rna <- L1000_compounds@molecularProfiles$rna[,myx]

## do a loop through cell lines
signatures <- list()

for (cell.line in cell.lines) {
    subsetted.drugs <- subsetTo(L1000_compounds, drugs=common.drugs, cells=cell.line)
    
    signatures[[cell.line]] <- drugPerturbationSig(subsetted.drugs, mDataType="rna", nthread=20)
}

saveRDS(signatures, "signatures.RData")