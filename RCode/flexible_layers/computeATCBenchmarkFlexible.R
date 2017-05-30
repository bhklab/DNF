ComputeATCBenchmarkFlexible <- function(atc.benchmark.name, common.drugs, strc.aff.mat, sens.aff.mat,
                                        pert.aff.mat, integration, luminex.aff.mat, imaging.aff.mat, atc.roc.file.name) {
    # 2. ATC benchmarking
    
    atc.bench <- ATCBenchFlexible("chembl-new", common.drugs)
    
    pairs.atc <- GenerateDrugPairsFlexible(atc.bench, strcAff=strc.aff.mat, sensAff=sens.aff.mat,
                                           pertAff=pert.aff.mat, integration=integration,
                                           luminexAff=luminex.aff.mat, imagingAff=imaging.aff.mat)
    
    res.atc <- CompConcordIndxFlexible(pairs.atc)
    
    PrintCIndices(res.atc$c.index.list)
    PrintPVals(res.atc$p.vals.list)
    
    GenerateROCPlotFlexible(pairs.atc, atc.roc.file.name, nrow(atc.bench))
}