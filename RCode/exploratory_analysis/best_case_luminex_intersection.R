load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

l1000 <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)


luminex_meta <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.meta.cpd.txt", stringsAsFactors = FALSE)
luminex <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.profiles.txt", stringsAsFactors = FALSE)

datasets = list(CCLE=CCLE, cTRPv2=CTRPv2, FIMM=FIMM, gCSI=gCSI, GDSC1000=GDSC1000)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

all.names <- sapply(datasets, 
                    function(x) {
                        if (!is.null(drugInfo(x)$broad_cpd_id)) {
                            return(list(drug_name=drugNames(x)))
                        } else {
                            return(list(drug_name=drugNames(x)))
                        }})

all.names <- unlist(all.names)
all.names <- toupper(all.names)
all.names <- gsub(badchars, "", all.names)

luminex.names <- luminex_meta$name
luminex.names <- toupper(luminex.names)
luminex.names <- gsub(badchars, "", luminex.names)

luminex_meta_rows_name <- luminex_meta[luminex.names %in% all.names,]
luminex_meta_rows_broad <- luminex_meta[luminex_meta$name %in% all.names,]

drug_mappings <- read.csv("Data/drugs_with_ids.csv", stringsAsFactors = FALSE)
drug_mapping_clean_id <- drug_mappings$unique.drugid
drug_mapping_clean_id <- toupper(drug_mappings$unique.drugid)
drug_mapping_clean_id <- gsub(badchars, "", drug_mapping_clean_id)

smiles <- drug_mappings[drug_mapping_clean_id %in% all.names, "smiles"]


for (i in 1:length(datasets)) {
    View(drugInfo(datasets[[i]]))
}
