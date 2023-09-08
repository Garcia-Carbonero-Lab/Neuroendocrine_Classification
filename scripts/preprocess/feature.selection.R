# This script needs preprocess environment
library(tidyverse)


source("functions/preprocess/aditional.functions.R")

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

# Open expression data with quality adjusted

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

# Open methylation data with quality adjusted

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

#Filtering using the genes with more
#MAD taking into account primary tumor separately


if (dir.exists(paste0(wkdir, "/preprocess/mad")) == F) {

    dir.create(paste0(wkdir, "/preprocess/mad"))

}

setwd(paste0(wkdir, "/preprocess/mad"))

tissues <- c("COLORECTAL", "PANCREAS", "LUNG", "GASTRIC", "SMALL INTESTINE")

expression_fs <- filter_tissue(expression, clinical.data,
tissues = tissues, flag = "expression")

methylation_fs <- filter_tissue(methylation, clinical.data,
tissues = tissues, flag = "methylation")



write.table(expression_fs,
paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
sep = "\t",
row.names = T,
col.names = NA
)



write.table(methylation_fs,
paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
sep = "\t",
row.names = T,
col.names = NA
)
