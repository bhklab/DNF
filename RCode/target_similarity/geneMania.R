library(org.Hs.eg.db)
library(biomaRt)

temp <- read.delim("Data/target_similarity/Shared_protein_domains.PFAM.txt", header = T, sep="", stringsAsFactors = F)

gene.names <- unique(union(temp$Gene_A, temp$Gene_B))
targets <- unique(data.bench$TARGET_NAME)
targets <- trimws(targets)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"),
       filters = "hgnc_symbol", values = targets, mart = human)
unique.mapping <- mapping[!duplicated(mapping$hgnc_symbol), ]

unique.mapping <- rbind(unique.mapping, c("AIK", "ENSG00000087586"))
unique.mapping <- rbind(unique.mapping, c("NR1B1", "ENSG00000131759"))
unique.mapping <- rbind(unique.mapping, c("TUBA2", "ENSG00000198033"))
unique.mapping <- rbind(unique.mapping, c("PDGFR", "ENSG00000113721"))
unique.mapping <- rbind(unique.mapping, c("FKBP1", "ENSG00000088832"))
unique.mapping <- rbind(unique.mapping, c("FLK1", "ENSG00000128052"))
unique.mapping <- rbind(unique.mapping, c("CDC2L4", "ENSG00000136807"))
unique.mapping <- rbind(unique.mapping, c("BFGFR", "ENSG00000077782"))
unique.mapping <- rbind(unique.mapping, c("CDHF12", "ENSG00000165731"))
unique.mapping <- rbind(unique.mapping, c("POLA", "ENSG00000101868"))
unique.mapping <- rbind(unique.mapping, c("TOP2", "ENSG00000131747"))
unique.mapping <- rbind(unique.mapping, c("CAK", "ENSG00000134058"))
unique.mapping <- rbind(unique.mapping, c("CSBP", "ENSG00000112062"))
unique.mapping <- rbind(unique.mapping, c("ABL", "ENSG00000097007"))
unique.mapping <- rbind(unique.mapping, c("CD135", "ENSG00000122025"))
unique.mapping <- rbind(unique.mapping, c("CDC2", "ENSG00000170312"))
unique.mapping <- rbind(unique.mapping, c("AIK2", "ENSG00000178999"))
unique.mapping <- rbind(unique.mapping, c("PKC2", "ENSG00000067606"))
unique.mapping <- rbind(unique.mapping, c("P53R2", "ENSG00000048392"))
unique.mapping <- rbind(unique.mapping, c("RAF", "ENSG00000132155"))
unique.mapping <- rbind(unique.mapping, c("PLK", "ENSG00000166851"))
unique.mapping <- rbind(unique.mapping, c("FRAP", "ENSG00000198793"))

shared.protein.domain <- temp
shared.protein.domain <- shared.protein.domain[shared.protein.domain$Gene_A %in% unique.mapping$ensembl_gene_id & shared.protein.domain$Gene_B %in% unique.mapping$ensembl_gene_id, ]

all.ensembl <- unique(unique.mapping$ensembl_gene_id)
all.genes <- unique(unique.mapping$hgnc_symbol)

network <- matrix(0, nrow=length(all.genes), ncol=length(all.genes))
rownames(network) <- all.genes
colnames(network) <- all.genes

for (i in 1:nrow(shared.protein.domain)) {
    e1 <- shared.protein.domain[i, "Gene_A"]
    e2 <- shared.protein.domain[i, "Gene_B"]
    
    target.similarity <- shared.protein.domain[i, "Weight"]
    
    genes.with.ensembl.1 <- unique.mapping[unique.mapping$ensembl_gene_id == e1, "hgnc_symbol"]
    genes.with.ensembl.2 <- unique.mapping[unique.mapping$ensembl_gene_id == e2, "hgnc_symbol"]
    
    
    network[genes.with.ensembl.1, genes.with.ensembl.2] <- target.similarity
    network[genes.with.ensembl.2, genes.with.ensembl.1] <- target.similarity
}

network[network == 0] <- -Inf
diag(network) <- 1
