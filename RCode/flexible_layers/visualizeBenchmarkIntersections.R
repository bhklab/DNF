library(VennDiagram)
library(org.Hs.eg.db)
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/drugTargetBenchModded.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

pert.file.name <- "Data/L1000_compound_signatures.RData"
sensitivity.file.name <- "Data/combined_sens_adjusted_diag_inamedatasets_with.RData"
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

sensData <- NULL
pertData <- NULL
strcData <- NULL
luminexData <- NULL
imagingData <- NULL

sensData <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
dim(sensData) # 309 x 309 for combined data    
rownames(sensData) <- toupper(rownames(sensData))
rownames(sensData) <- gsub(badchars, "", rownames(sensData))

colnames(sensData) <- toupper(colnames(sensData))
colnames(sensData) <- gsub(badchars, "", colnames(sensData))

pertData <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 239
print(dim(pertData)) # 978 x 237 for         

colnames(pertData) <- toupper(colnames(pertData))
colnames(pertData) <- gsub(badchars, "", colnames(pertData))

pertNames <- colnames(pertData)

#pertNames <- lincs.meta$pert_iname
#pertNames <- colnames(sensData)

layers <- list(sensNames = sort(colnames(sensData)), pertNames=pertNames)
commonDrugs <- Reduce(intersect, Filter(Negate(is.null),layers))

ctrpv2 <- DrugTargetBenchModded(commonDrugs, "temp.RData", use.ctrpv2=TRUE,
                                use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE)
ctrpv2 <- unique(colnames(ctrpv2))

clue <- DrugTargetBenchModded(commonDrugs, "temp.RData", use.ctrpv2=FALSE,
                              use.clue=TRUE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE)
clue <- unique(colnames(clue))

chembl <- DrugTargetBenchModded(commonDrugs, "temp.RData", use.ctrpv2=FALSE,
                                use.clue=FALSE, use.chembl=TRUE, use.dbank=FALSE, use.dtc=FALSE)
chembl <- unique(colnames(chembl))

dbank <- DrugTargetBenchModded(commonDrugs, "temp.RData", use.ctrpv2=FALSE,
                               use.clue=FALSE, use.chembl=FALSE, use.dbank=TRUE, use.dtc=FALSE)
dbank <- unique(colnames(dbank))

dtc <- DrugTargetBenchModded(commonDrugs, "temp.RData", use.ctrpv2=FALSE,
                             use.clue=FALSE, use.chembl=FALSE, use.dbank=FALSE, use.dtc=TRUE)
dtc <- unique(colnames(dtc))

target.list <- list(ctrpv2=ctrpv2, clue=clue, dbank=dbank)

names.with.count <- toupper(paste(names(target.list), sapply(target.list, length)))


union.size <- length(unique(Reduce(union, target.list)))
main.title = paste("3 Layers Overlap Between CTRPv2-Clue-Drug Bank Benchmarks \n ", 
                   "Total Unique Drugs:",union.size)

vp <- venn.diagram(target.list,
                   fill = c(5,6,7), alpha = 0.3, filename = NULL, imagetype = "png",
                   category.names = names.with.count, cat.cex=2, cat.dist=0.1, margin=0.1, 
                   main=main.title, main.cex=2.2,
                   main.fontfamily="Arial")

old.par <- par(mar = c(0, 0, 0, 0))
plot.new()
frame()
par(old.par)
grid.draw(vp)
