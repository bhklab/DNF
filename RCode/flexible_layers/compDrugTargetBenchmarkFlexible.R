CompDrugTargetBenchmarkFlexible <- function(common.drugs, strc.aff.mat,
                                    sens.aff.mat, pert.aff.mat, integration, luminex.aff.mat, imaging.aff.mat,
                                    target.roc.file.name, use.ctrpv2,
                                    use.clue, use.chembl, use.dbank, use.dtc) {
    # 1. Drug target benchmarking

    data.bench <- DrugTargetBenchFlexible(common.drugs, 
                                       "temp.RData", use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl,
                                       use.dbank=use.dbank, use.dtc=use.dtc) # 141 x 141 drug-drug adjacency matrix --> 141
    
    pairs.target <- GenerateDrugPairsFlexible(data.bench, strc.aff=strc.aff.mat, sens.aff=sens.aff.mat,
                                              pert.aff=pert.aff.mat, integration=integration,
                                              luminex.aff=luminex.aff.mat, imaging.aff=imaging.aff.mat)
    
    res.target <- CompConcordIndxFlexible(pairs.target)
    
    PrintCIndices(res.target$c.index.list)
    PrintPVals(res.target$p.vals.list)
    
    GenerateROCPlotFlexible(pairs.target, target.roc.file.name, nrow(data.bench))
    
}