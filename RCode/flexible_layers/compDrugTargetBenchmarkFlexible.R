CompDrugTargetBenchmarkFlexible <- function(common.drugs, strc.aff.mat,
                                    sens.aff.mat, pert.aff.mat, integration, luminex.aff.mat, imaging.aff.mat,
                                    target.roc.file.name, use.ctrpv2,
                                    use.clue, use.chembl, use.dbank, use.dtc) {
    # 1. Drug target benchmarking

    dataBench <- drugTargetBenchModded(common.drugs, 
                                       "temp.RData", use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl,
                                       use.dbank=use.dbank, use.dtc=use.dtc) # 141 x 141 drug-drug adjacency matrix --> 141
    
    pairs.target <- GenerateDrugPairsFlexible(dataBench, strcAff=strc.aff.mat, sensAff=sens.aff.mat,
                                              pertAff=pert.aff.mat, integration=integration,
                                              luminexAff=luminex.aff.mat, imagingAff=imaging.aff.mat)
    
    res.target <- CompConcordIndxFlexible(pairs.target)
    
    PrintCIndices(res.target$c.index.list)
    PrintPVals(res.target$p.vals.list)
    
    GenerateROCPlotFlexible(pairs.target, target.roc.file.name, nrow(dataBench))
    
}