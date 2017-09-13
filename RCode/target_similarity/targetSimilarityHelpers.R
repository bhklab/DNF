GetKeggIdentifierForGene <- function(gene) {
    tryCatch({
        res <- KEGGREST::keggFind("genes", gene)
    }, error = function(e) {
        return(NA)
    })
    
    for (i in 1:length(res)) {
        temp <- res[i]
        
        if (startsWith(temp, gene)) {
            return(names(res)[i])
        }
    }
    
    return(names(res)[1])
}

GetKeggIdentifiersForGenes <- function(kegg.identifiers, genes) {
    identifiers <- rep(NA, length(genes))
    names(identifiers) <- genes
    
    for (gene in genes) {
        relevant.identifiers <- grep(gene, kegg.identifiers)
        
        if (length(relevant.identifiers) == 1) {
            identifiers[gene] <- names(kegg.identifiers[relevant.identifiers])
            next()
        }
        
        for (i in 1:length(relevant.identifiers)) {
            temp <- kegg.identifiers[relevant.identifiers[i]]
            
            
            tryCatch({
                if (startsWith(temp, paste(gene, ",", sep="")) || startsWith(temp, paste(gene, ";", sep="")) || 
                    length(grep(paste(" ", gene, ",", sep=""), temp)) > 0 ||
                    length(grep(paste(" ", gene, ";", sep=""), temp)) > 0) {
                    identifiers[gene] <- names(kegg.identifiers[relevant.identifiers[i]])
                    break
                }                
            }, error = function(e) {
                print(gene)
                print(relevant.identifiers)
            })

        }
    }
    
    return(identifiers)
}

GetAminoAcidSequences <- function(kegg.identifiers) {
    tryCatch({
        aa.sequences <- list()
        buckets <- split(kegg.identifiers, ceiling(seq_along(kegg.identifiers) / 10))
        
        for (i in 1:length(buckets)) {
            bucket <- buckets[[i]]
            
            aa.sequences[[i]] <- KEGGREST::keggGet(bucket, "aaseq")
        }
    }, error = function(e) {
        return(NA)
    })
    
    return(list(aa.sequences=aa.sequences, buckets=buckets))
}

GetAminoAcidSequencesSlow <- function(kegg.identifiers) {
    aa.sequences <- list()
    
    for (i in 1:length(kegg.identifiers)) {
        identifier <- kegg.identifiers[i]
        
        tryCatch({
            aa.sequences[[identifier]] <- KEGGREST::keggGet(identifier, "aaseq")
        }, error = function(e) {
            aa.sequences[[identifier]] <- NA
            print(names(kegg.identifiers)[i])
            print(identifier)
        })     
        
        print(i)
    }
    
    return(aa.sequences)
}

NormalizedSW <- function(s1, s2) {
    numerator <- pairwiseAlignment(s1, s2, type="local")@score
    
    denominator <- sqrt(pairwiseAlignment(s1, s1, type="local")@score) * sqrt(pairwiseAlignment(s2, s2, type="local")@score)
    
    return(numerator / denominator)
}