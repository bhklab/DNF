# Integrative pharmacogenomics to infer large-scale drug taxonomy 

Please cite: Integrative Cancer Pharmacogenomics to Infer Large-Scale Drug Taxonomy.
El-Hachem N, Gendoo DMA, Ghoraie LS, Safikhani Z, Smirnov P, Chung C, Deng K, Fang A, Birkwood E, Ho C, Isserlin R, Bader GD, Goldenberg A, Haibe-Kains B.
Cancer Res. 2017 Jun 1;77(11):3057-3069. doi: 10.1158/0008-5472.CAN-17-0096. Epub 2017 Mar 17.
PMID: 28314784

This project assess the utility of data integration from multiple drug information layers.

Hypothesis: Better characterization of compound mechanism of action (MoA) from integrative pharmacogenomics.

Impact: Identifying mechanism for compounds with unknown MoA in early stages of drug development without the need 
of sophisticated info (side effects or pharmacological indications...) 

# Based on the method described in 
Similarity network fusion for aggregating data types on a genomic scale.
Wang B, Mezlini AM, Demir F, Fiume M, Tu Z, Brudno M, Haibe-Kains B, Goldenberg A.
Nat Methods. 2014 Mar;11(3):333-7. doi: 10.1038/nmeth.2810. Epub 2014 Jan 26.
PMID: 24464287

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

A Docker Terminal Environment has been setup with all dependencies installed. Access it by installing [Docker Engine](https://www.docker.com/community-edition) and running the following lines of code:

```
#clone image
docker pull gosuzombie/dnf

#run image in interactive command line 
docker run -it gosuzombie/dnf
```

Alternatively, one can setup the Docker environment by using the provided Dockerfile. Instructions can be found [here](https://docs.docker.com/engine/getstarted/step_four/) 
