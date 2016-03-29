# DNF (Drug Network Fusion)
===========================

This project assess the utility of data integration from multiple drug information layers.

Hypothesis: Better characterization of compound mechanism of action (MoA) from integrative pharmacogenomics.

Impact: Identifying mechanism for compounds with unknown MoA in early stages of drug development without the need 
of sophisticated info (side effects or pharmacological indications...) 

# Based on the method described in http://www.ncbi.nlm.nih.gov/pubmed/24464287

Please cite: coming soon

# Reproducibility of the Analysis 

The following steps will reproduce the output files (figure, tables...) mentioned in the main manuscript.
The script will be using data files such as:

## Files to run the scripts 

- Drug-target (benchmark) from CHEMBL and CTRP
- Sensitivity measurements (drugs x cell lines) from CTRPv2 and NCI-60
- Transcriptomic data from the LINCS database http://lincs.hms.harvard.edu/ created using PharmacoGx package https://cran.r-project.org/web/packages/PharmacoGx/index.html

## Run the R scripts 

* main-ctrpv-lincs.R and 
* main-nci60-lincs.R

(this will execute the following R codes):

```
# process the raw data, find common drugs...
preprocessInput.R 	

# remove problematic drugs and missing data
sensitivityData.R 	

# get the structural fingerprints (extended-connectivity descriptors from RCDK)
structureData.R 

# remove problematic drugs and missing data and find drug names from LINCS metadata
perturbationData.R 	

# Get the similarity matrix from structure (tanimoto metric)
constStructureLayer.R 	

# Get the similarity matrix from sensitivity (pearson metric), rescale to 0-1 with function in SNF
constSensitivityLayer.R	

# Get the similarity matrix from perturbation (pearson metric), rescale to 0-1 with SNF affinitymatrix
constPerturbationLayer.R

# Integrate all 3 layers using SNFtools
integrateStrctSensPert.R	

# ATC code gold standard benchmark from CHEMBL
ATCbench.R	

# Drug-target benchmarks from CHEMBL and CTRPv2
drugTargetBench.R	

# get the pairs of drugs sharing the same target (0 or 1)
generateDrugPairs.R	

# Compute concordance index between benchmark data and similarity network data
compConcordIndx.R	

# Get the p-value between concordance indices
predPerf.R	

# Generate ROC plots and AUC values
generateRocPlot.R	

# Community clustering from apcluster package in R
communityGen.R

# Network based on exemplar drugs from apcluster
main-network-generation.R	

```

# Set up the software environment

```
R version 3.2.1 (2015-06-18)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.9.5 (Mavericks)

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] proxy_0.4-15         reshape2_1.4.1       minerva_1.4.5        survcomp_1.18.0      prodlim_1.5.7       
 [6] survival_2.38-3      ROCR_1.0-7           gplots_2.17.0        SNFtool_2.2          org.Hs.eg.db_3.1.2  
[11] RSQLite_1.0.0        DBI_0.3.1            annotate_1.46.1      XML_3.98-1.3         AnnotationDbi_1.30.1
[16] GenomeInfoDb_1.4.3   IRanges_2.2.9        S4Vectors_0.6.6      Biobase_2.28.0       BiocGenerics_0.14.0 
[21] rcdk_3.3.2           fingerprint_3.5.2    apcluster_1.4.2      PharmacoGx_1.0.6    

loaded via a namespace (and not attached):
 [1] gtools_3.5.0       slam_0.1-32        sets_1.0-16        splines_3.2.1      rJava_0.9-7        lattice_0.20-33   
 [7] marray_1.46.0      sm_2.2-5.4         piano_1.8.2        magicaxis_1.9.4    RColorBrewer_1.1-2 heatmap.plus_1.3  
[13] plyr_1.8.3         stringr_1.0.0      lava_1.4.1         survivalROC_1.0.3  rcdklibs_1.5.8.4   caTools_1.17.1    
[19] Rcpp_0.12.3        KernSmooth_2.23-15 xtable_1.8-0       relations_0.6-6    limma_3.24.15      gdata_2.17.0      
[25] rmeta_2.16         plotrix_3.6-1      bootstrap_2015.2   png_0.1-7          digest_0.6.9       stringi_1.0-1     
[31] SuppDists_1.1-9.1  grid_3.2.1         tools_3.2.1        bitops_1.0-6       magrittr_1.5       cluster_2.0.3     
[37] MASS_7.3-45        Matrix_1.2-3       downloader_0.4     iterators_1.0.8    igraph_1.0.1 

```