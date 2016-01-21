###############################################################################################################
## Function reads in the structure data, and generates "affinity matrix" for the structure layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     targFps: target finger prints generated in "structureData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


constStructureLayer <-  function(targFps) {
     ## Correlation for Structure (Tanimoto metric)
     fpSim <- fingerprint::fp.sim.matrix(targFps, method = "tanimoto")
     rownames(fpSim) <- names(targFps)
     colnames(fpSim) <- names(targFps)
     ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
     fpAff <- SNFtool::affinityMatrix(1-fpSim, 20, 0.5)
     return(fpAff)
}