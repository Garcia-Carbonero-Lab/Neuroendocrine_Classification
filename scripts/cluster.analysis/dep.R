#This script needs the environment deg


library(limma)
library(ggplot2)
library(rsample)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

#load functions

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/cluster.analysis/methylation.R")

# we read config file
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

if (dir.exists(paste0(wkdir, "/clustering/dep")) == F) {

    dir.create(paste0(wkdir, "/clustering/dep"))

}


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

base$PRIMARY_TUMOR <- as.factor(gsub("_", " ", base$PRIMARY_TUMOR))

# open methylation data
data <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

base$PRIMARY_TUMOR <- as.factor(gsub(" ", "_", base$PRIMARY_TUMOR))


# differentially methylated probes
deg.function(data,
base,
group = "Subtype",
cov = "PRIMARY_TUMOR",
outdir = paste0(wkdir, "/clustering/dep"),
flag = "dep_cov_primary")

# open results of limma
dep1 <- read.table(paste0(wkdir,
"/clustering/dep/TopTable_dep_cov_primary_Group1-Group2.txt"),
sep = "\t",
row.names = 1,
header = T
)

dep2 <- read.table(paste0(wkdir,
"/clustering/dep/TopTable_dep_cov_primary_Group1-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

dep3 <- read.table(paste0(wkdir,
"/clustering/dep/TopTable_dep_cov_primary_Group2-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

# anotation of methylation data
epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)
    

# transform methylation anotation matrix
met_anot <- create.met.anot(epic)

# anotate limma results with genes
dep1 <- anotation(dep1,
met_anot,
outdir = paste0(wkdir, "/clustering/dep"),
flag = "Group1-Group2")

dep2 <- anotation(dep2,
met_anot,
outdir = paste0(wkdir, "/clustering/dep"),
flag = "Group1-Group3")

dep3 <- anotation(dep3,
met_anot,
outdir = paste0(wkdir, "/clustering/dep"),
flag = "Group2-Group3")

# filter results
deps1 <- dep1$ID[dep1$adj.P.Val < 0.05 & abs(dep1$logFC) > 2]
deps2 <- dep2$ID[dep2$adj.P.Val < 0.05 & abs(dep2$logFC) > 2]
deps3 <- dep3$ID[dep3$adj.P.Val < 0.05 & abs(dep3$logFC) > 2]



probes <- c(deps1,deps2,deps3)

met.deps <- data[unique(probes),]

write.table(met.deps,
paste0(wkdir,
"/clustering/dep/Methylation.dep.txt"),
sep = "\t",
row.names= T,
col.names= NA)

