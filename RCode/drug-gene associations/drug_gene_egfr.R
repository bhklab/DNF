########## R code to test gene-drug associations in CTRPv2 (EGFR inhibitors) ################
library(PharmacoGx)
library(affy)
setwd("/Users/hachemn/Desktop/TEST1")

#### load sensitivity data (AUCs from CTRPv2 study for the 9 drugs) file provided on github
load("mat_ctrp2.RData")

# download the Pset using pharmacoGx to get basal gene expression profiles
ccle <- downloadPSet("CCLE") 

### get the gene expression data for CCLE cell lines (same for the CTRPv2 study) 
ccle.rna <- summarizeMolecularProfiles(ccle,
                                       mDataType = "rna",
                                       summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
ccle.rna <- exprs(ccle.rna)
ccle.rna <- t(ccle.rna)

### intersect with ctrpv2 sensitivity data to get common cell lines
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
rownames(ccle.rna) <- gsub(badchars,"",rownames(ccle.rna))
rownames(ccle.rna) <- toupper(rownames(ccle.rna))
intcell <- intersect(rownames(ccle.rna), rownames(mat.ctrp2))

mat.ctrp2 <- mat.ctrp2[intcell,]
ccle.rna <- ccle.rna[intcell,]

#### convert ENSEMBL to GENE SYMBOL
library('biomaRt')
colnames(ccle.rna) <- gsub("_at","",colnames(ccle.rna))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- colnames(ccle.rna)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list2 <- G_list[!G_list$hgnc_symbol=="",,drop=F]
G_list2 <- G_list2[!duplicated(G_list2$ensembl_gene_id),,drop=F]
intgx <- intersect(colnames(ccle.rna), G_list2[,1])
ccle.rna2 <- ccle.rna[,intgx]
colnames(ccle.rna2) <- G_list2[,2]

##### correlate gene expression with sensitivity data for the 9 drugs
corx <- cor(mat.ctrp2, ccle.rna2, method="pearson", use="pairwise.complete.obs")

##### transpose the correlation matrix and change the sign since in CTRPv2 high AUC values mean resistance
corx <- -t(corx)

#### load a gmt file for GSVA, included in github
load("gmt_canonical_msigdb.RData")

### run GSVA with the ssGSEA setting
library(GSVA)
ssgsea.sens <- gsva(corx, gmtlistclass,min.sz=15, max.sz=250,method="ssgsea")

### get the median score across all drugs for a given pathway
med.pathway <- apply(ssgsea.sens, 1, median, na.rm=TRUE)
med.pathway <- matrix(med.pathway)
med.pathway <- data.frame("median_path"=med.pathway,ssgsea.sens)

##### plot the heatmap with the 10 top +/- scores

#### rank by ssgsea median score
med.pathway <- med.pathway[order(med.pathway[,"median_path"], decreasing = T),,drop=F]
top.pathway <- rbind(head(med.pathway,10), tail(med.pathway,10))
top.pathway <- top.pathway[,-1]

### plot the heatmap
library("gplots")
pdf("ssgsea_pathways_ctrpv2_sensitivity.pdf")

heatmap.2(as.matrix(top.pathway), col=c("#31a354","#de2d26"), scale="none", dendrogram = "col", margins =c(7,20),
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.3, cexCol = 0.7,
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=1:ncol(top.pathway),
          rowsep=1:nrow(top.pathway))
dev.off()