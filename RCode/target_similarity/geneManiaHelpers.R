LoadGeneManiaSimilarityData <- function(file.path) {
    similarity.data <- read.delim(file.path, sep="", stringsAsFactors = F)
    
    return(similarity.data)
}

CreateSymbolToEnsemblMapping <- function(targets) {
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mapping <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"),
                     filters = "hgnc_symbol", values = targets, mart = human)
    mapping <- rbind(mapping, c("AIK", "ENSG00000087586"))
    mapping <- rbind(mapping, c("NR1B1", "ENSG00000131759"))
    mapping <- rbind(mapping, c("TUBA2", "ENSG00000198033"))
    mapping <- rbind(mapping, c("PDGFR", "ENSG00000113721"))
    mapping <- rbind(mapping, c("FKBP1", "ENSG00000088832"))
    mapping <- rbind(mapping, c("FLK1", "ENSG00000128052"))
    mapping <- rbind(mapping, c("CDC2L4", "ENSG00000136807"))
    mapping <- rbind(mapping, c("BFGFR", "ENSG00000077782"))
    mapping <- rbind(mapping, c("CDHF12", "ENSG00000165731"))
    mapping <- rbind(mapping, c("POLA", "ENSG00000101868"))
    mapping <- rbind(mapping, c("TOP2", "ENSG00000131747"))
    mapping <- rbind(mapping, c("CAK", "ENSG00000134058"))
    mapping <- rbind(mapping, c("CSBP", "ENSG00000112062"))
    mapping <- rbind(mapping, c("ABL", "ENSG00000097007"))
    mapping <- rbind(mapping, c("CD135", "ENSG00000122025"))
    mapping <- rbind(mapping, c("CDC2", "ENSG00000170312"))
    mapping <- rbind(mapping, c("AIK2", "ENSG00000178999"))
    mapping <- rbind(mapping, c("PKC2", "ENSG00000067606"))
    mapping <- rbind(mapping, c("P53R2", "ENSG00000048392"))
    mapping <- rbind(mapping, c("RAF", "ENSG00000132155"))
    mapping <- rbind(mapping, c("PLK", "ENSG00000166851"))
    mapping <- rbind(mapping, c("FRAP", "ENSG00000198793"))
    
    return(mapping)
}

SubsetSimilarities <- function(target.similarities, mapping) {
    target.similarities <- target.similarities[target.similarities$Gene_A %in% mapping$ensembl_gene_id & target.similarities$Gene_B %in% mapping$ensembl_gene_id, ]
    
    return(target.similarities)
}

CreateTargetSimilarityNetwork <- function(target.similarities, mapping) {
    all.ensembl <- unique(mapping$ensembl_gene_id)
    all.genes <- unique(mapping$hgnc_symbol)
    
    network <- matrix(0, nrow=length(all.genes), ncol=length(all.genes))
    rownames(network) <- all.genes
    colnames(network) <- all.genes
    
    for (i in 1:nrow(target.similarities)) {
        e1 <- target.similarities[i, "Gene_A"]
        e2 <- target.similarities[i, "Gene_B"]
        
        target.similarity <- target.similarities[i, "Weight"]
        
        genes.with.ensembl.1 <- mapping[mapping$ensembl_gene_id == e1, "hgnc_symbol"]
        genes.with.ensembl.2 <- mapping[mapping$ensembl_gene_id == e2, "hgnc_symbol"]
        
        
        network[genes.with.ensembl.1, genes.with.ensembl.2] <- target.similarity
        network[genes.with.ensembl.2, genes.with.ensembl.1] <- target.similarity
    }
    
    diag(network) <- 1
    
    return(network)
}

