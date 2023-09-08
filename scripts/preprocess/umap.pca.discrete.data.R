# For use this script we need the environment preprocess

# we open libraries
library(umap)
library(ggplot2)
library(limma)
library(RColorBrewer)
# open functions

source("functions/preprocess/aditional.functions.R")
source("functions/main.functions.R")


# Open config file
config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

# Open data

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)


methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)


clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
sep = '\t',
header = T,
row.names = 1
)

###########################################################################

#Create paths
if (dir.exists(paste0(wkdir, "/preprocess/discrete_features")) == F) {

    dir.create(paste0(wkdir, "/preprocess/discrete_features"))

}

if (dir.exists(paste0(wkdir,
"/preprocess/discrete_features/expression")) == F) {

    dir.create(paste0(wkdir, "/preprocess/discrete_features/expression"))

}


if (dir.exists(paste0(wkdir,
"/preprocess/discrete_features/methylation")) == F) {

    dir.create(paste0(wkdir, "/preprocess/discrete_features/methylation"))

}

#Select features to test
columns <- c("BIOPSY_PRIMARY_LOCATION",
            "BIOPSY_PRIMARY_LOCATION",
            "BIOPSY_METASTASIS_LOCATION",
            "GENDER",
            "ECOG_PS",
            "PRIMARY_TUMOR",
            "GRADE",
            "DIAGNOSIS_STAGE_INITIAL",
            "FUNCTIONAL_dichotomic",
            "X5HIAA_HIGH",
            "BIOPSY_COMPLETE_LOCATION",
            "BIOPSY_LOCATION",
            "CIMP.status")

pca_exp_qual <- pca_corr_discrete(data = expression,
base = clinical.data,
columns = c("TRANSCRIPTOME_BATCH", columns),
outpath = paste0(wkdir, "/preprocess/discrete_features/expression"),
flag = "expression",
anotate1 = c(50, 90),
anotate2 = c(50, 80)

)

umap_exp_qual <- umap_corr_discrete(data = expression,
base = clinical.data,
columns = c("TRANSCRIPTOME_BATCH", columns),
outpath = paste0(wkdir, "/preprocess/discrete_features/expression"),
flag = "expression",
anotate1 = c(1.8, 4),
anotate2 = c(1.8, 3.7)
)

write.table(pca_exp_qual,
paste0(wkdir, 
"/preprocess/discrete_features/expression/Correlation.PCA.expression.tsv"),
sep = "\t",
row.names = T
,col.names = NA
)

write.table(umap_exp_qual,
paste0(wkdir,
"/preprocess/discrete_features/expression/Correlation.UMAP.expression.tsv"),
sep = "\t",
row.names = T
, col.names = NA
)
#############################################
#Methylation

pca_met_qual <- pca_corr_discrete(data = methylation,
base = clinical.data,
columns = c("METHYLOMA_BATCH", columns),
outpath = paste0(wkdir, "/preprocess/discrete_features/methylation"),
flag = "metrhylation",
anotate1 = c(350, 500),
anotate2 = c(350, 470)

)

umap_met_qual <- umap_corr_discrete(data = methylation,
base = clinical.data,
columns = c("METHYLOMA_BATCH", columns),
outpath = paste0(wkdir, "/preprocess/discrete_features/methylation"),
flag = "metrhylation",
anotate1 = c(2.5, 3),
anotate2 = c(2.5, 2.77)
)

write.table(pca_met_qual,
paste0(wkdir, 
"/preprocess/discrete_features/methylation/Correlation.PCA.methylation.tsv"),
sep = "\t",
row.names = T
,col.names = NA
)

write.table(umap_met_qual,
paste0(wkdir,
"/preprocess/discrete_features/methylation/Correlation.UMAP.methylation.tsv"),
sep = "\t",
row.names = T
, col.names = NA
)
