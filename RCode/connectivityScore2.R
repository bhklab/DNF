

connectivityScore2 <- function(x, y, method=c("gsea", "gwc"), nperm=1e4, nthread=1, gwc.method=c("spearman", "pearson"), ...) {
  
  method <- match.arg(method)
  if (class(x) != "matrix") {
    x <- as.matrix(x)
  }
  if (class(y) != "matrix") {
    y <- as.matrix(y)
  }
  if ((ncol(x) != 2 || ncol(y) != 2) && method=="gwc") {
    stop ("x and y should have 2 columns: effect size and corresponding p-values")
  }
  
  if (method == "gsea" && nrow(y) >= nrow(x)) {
    warning("GSEA method: query gene set (y) larger than signature (x)")
  }
  
  if (is.null(rownames(x)) || is.null(rownames(y)) || !length(intersect(rownames(x), rownames(y)))) {
    stop ("Row names of x and y are either missing or have no intersection")
  }
 # if (nperm < 100){
  #  stop ("The minimum number of permutations for permutation testing is 100")
  #}
  switch (method,
          "gsea" = {
            ## remove missing values
            y <- y[!is.na(y[ ,1]), , drop=FALSE]
            x <- x[!is.na(x[ ,1]), , drop=FALSE]
            ## create gene set
            gset <- cbind("gene"=rownames(y), "set"=ifelse(as.numeric(y[ , 1]) >= 0, "UP", "DOWN")) 
            gset <- piano::loadGSC(gset)
            ## run enrichment analysis
            nes <- runGSA2(geneLevelStats=x[ , 1], geneSetStat="gsea", gsc=gset, nPerm=nperm + (nperm %% nthread), ncpus=nthread, verbose=FALSE, ...)
            ## merge p-values for negative and positive enrichment scores
            nes$pDistinctDir <- nes$pDistinctDirUp
            nes$pDistinctDir[is.na(nes$pDistinctDirUp), 1] <- nes$pDistinctDirDn[is.na(nes$pDistinctDirUp), 1]
            nes.up <- c(nes$statDistinctDir[which(names(nes$gsc) == "UP"), 1], nes$pDistinctDir[which(names(nes$gsc) == "UP"), 1])
            nes.down <- c(nes$statDistinctDir[which(names(nes$gsc) == "DOWN"), 1], nes$pDistinctDir[which(names(nes$gsc) == "DOWN"), 1])
            ## combine UP and DOWN
            if (length(nes.up) == 0){
              score = c("es" = -nes.down[1], "p" = nes.down[2])
            } else if (length(nes.down) == 0){
              score = c("es" = nes.up[1], "p" = nes.up[2])
            } else if (complete.cases(cbind(nes.up, nes.down)) && sign(nes.up[1]) != sign(nes.down[1])) {
              score <- c("es"=(nes.up[1] - nes.down[1]) / 2, "p"= combineTest(p=c(nes.up[2], nes.down[2]), method="fisher", na.rm=TRUE))
            } else {
              score <- c("score"=0, "p"=1)
            }
          },
          "gwc" = {
            ## intersection between x and y
            ii <- intersect(rownames(x), rownames(y))
            if(length(ii) < 10) {
              stop ("Less than 10 probes/genes in common between x and y")
            }
            score <- gwc(x1=x[ii, 1], p1=x[ii, 2], x2=y[ii, 1], p2=y[ii, 2], method.cor=gwc.method, nperm=nperm, ...)
            names(score) <- c("score", "p")
          }
  )
  return (score)
}