source("https://bioconductor.org/biocLite.R")

list.of.CRAN.packages <- c("reshape2", "apcluster", "rcdk", "fingerprint", "SNFtool", "ROCR", "proxy", "PRROC", "Hmisc")
new.CRAN.packages <- list.of.CRAN.packages[!(list.of.CRAN.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.CRAN.packages, repos='http://cran.utstat.utoronto.ca/', dependencies = TRUE, clean = TRUE, Ncpus = 2, verbose = TRUE, quiet = TRUE)


list.of.bioC.packages <- c("PharmacoGx", "annotate", "org.Hs.eg.db", "survcomp")
new.bioC.packages <- list.of.bioC.packages[!(list.of.bioC.packages %in% installed.packages()[,"Package"])]
biocLite(pkgs = new.bioC.packages, Ncpus = 2, ask = FALSE)

dependencies <- c(list.of.CRAN.packages, list.of.bioC.packages)

message("=== package install success status ===")

for(package in dependencies){
  message(paste(package, library(package, character.only = TRUE, quietly = TRUE, logical.return = TRUE, verbose = FALSE), sep = ": "))
}
