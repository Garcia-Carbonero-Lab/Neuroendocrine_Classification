# IMPORTANT
#For this script we need the environment clustering.ensemble.yml

library(MOFA2)
library(ComplexHeatmap)
library(diceR)
library(circlize)
library(foreach)
library(doParallel)
library(NMF)
library(ggplot2)
#load functions

source("functions/main.functions.R")
source("functions/clustering/clustering.functions.R")
source("functions/MOFA/adictional.functions.R")

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

#We create output folder
if (dir.exists(paste0(wkdir, "/clustering")) == F) {

    dir.create(paste0(wkdir, "/clustering"))

}

if (dir.exists(paste0(wkdir, "/clustering/clustering_process")) == F) {

    dir.create(paste0(wkdir, "/clustering/clustering_process"))

}

#load_models


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T)


#factors to be removed
factor.rm <- c("Factor3", "Factor6")


#algorithms  to use in ensemble clustering
algorithms <- c('hc', 'diana', 'km', 'pam',
'gmm', 'block', 'cmeans', 'hdbscan')


outdir <- paste0(wkdir, "/clustering/clustering_process")


mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

#Create the similarity matrix for each k and algorithm
clustering.ensemble(mod,
factor.rm,
algorithms = algorithms,
outdir = outdir,
scale = F)

#Read the results of each clustering in each permutation

clust_input <- readRDS(paste0(outdir, "/clust_input.rds"))


#Read the similarity matrix for each k and algorithms
comb <- readRDS(paste0(outdir, "/comb.rds"))


#Calculate the best K
consensus <- lapply(comb, consensus_matrix)
coph <- lapply(consensus,cophcor)
coph <- unlist(coph)
coph.plot(coph,outdir)
k <- as.numeric(names(coph[coph == max(coph)]))


#Obtain the class of each sample
class_scale <- clustering.ensemble.matrix(mod = mod,
factors.rm = factor.rm,
clust_input = clust_input,
k = k,
outdir = outdir)


class_scale <- read.table(paste0(outdir, "/Class.matrix_3.txt"),
sep = "\t",
header = T,
row.names = 1
)


#Change class name 
class_scale$Subtype[class_scale$Subtype == 1] <- 4
class_scale$Subtype[class_scale$Subtype == 3] <- 1
class_scale$Subtype[class_scale$Subtype == 4] <- 3

class_scale <- class_scale[base$CLINICAL_TRIAL_CODE,,drop = F]


base$Subtype <- class_scale$Subtype

write.table(base,
paste0(datadir, "/clinical_data.txt"),
sep = "\t",
row.names = F,
col.names = T)
