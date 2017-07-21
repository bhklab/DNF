###############################################################################################################
## Function reads in the structure data, and generates "affinity matrix" for the structure layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     targ.fps: target finger prints generated in "structureData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


ConstStructureLayerFlexible <-  function(targ.fps) {
    ## Correlation for Structure (Tanimoto metric)
    fp.sim <- fingerprint::fp.sim.matrix(targ.fps, method = "tanimoto")
    rownames(fp.sim) <- names(targ.fps)
    colnames(fp.sim) <- names(targ.fps)
    
    saveRDS(fp.sim, "Data/uploading_features/structure/structure_similarities.RData")
    ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
    fp.aff <- SNFtool::affinityMatrix(1-fp.sim, 20, 0.5)
    return(fp.aff)
}