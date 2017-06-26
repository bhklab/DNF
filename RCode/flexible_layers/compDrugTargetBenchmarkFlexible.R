CompDrugTargetBenchmarkFlexible <- function(common.drugs, gmt.file.name, strc.aff.mat,
                                    sens.aff.mat, pert.aff.mat, integration, luminex.aff.mat, imaging.aff.mat,
                                    target.roc.file.name, use.ctrpv2,
                                    use.clue, use.chembl, use.dbank, use.dtc) {
    # 1. Drug target benchmarking
    
    data.bench <- DrugTargetBenchFlexible(common.drugs, 
                                       gmt.file.name, use.ctrpv2=use.ctrpv2,
                                       use.clue=use.clue, use.chembl=use.chembl,
                                       use.dbank=use.dbank, use.dtc=use.dtc) # 141 x 141 drug-drug adjacency matrix --> 141
    
    print("Benchmark obtained")
    pairs.target <- GenerateDrugPairsFlexible(data.bench, strc.aff=strc.aff.mat, sens.aff=sens.aff.mat,
                                              pert.aff=pert.aff.mat, integration=integration,
                                              luminex.aff=luminex.aff.mat, imaging.aff=imaging.aff.mat)
    
    print("Drug pairs obtained")
    res.target <- CompConcordIndxFlexible(pairs.target)
    
    print("Concordance Index computed")
    PrintCIndices(res.target$c.index.list)
    PrintPVals(res.target$p.vals.list)
    saveRDS(res.target$bad.performers, "bad_performers.RData")
    saveRDS(pairs.target, "pairs.RData")
    
    GenerateROCPlotFlexible(pairs.target, target.roc.file.name, nrow(data.bench))
    
}