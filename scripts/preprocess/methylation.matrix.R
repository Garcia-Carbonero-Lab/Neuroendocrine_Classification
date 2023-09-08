# This script need enviroment champ

library(ChAMP)
library(tidyverse)
# open functions

source("functions/Preprocess/methylation.process.R")
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

clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    header = T
)

myImport <- champ.import(paste0(datadir, "/methylome/raw"),
    arraytype = "EPIC"
)


# we calculate
qual <- qual.met(myImport)

clinical.data <- merge(clinical.data,
    qual,
    by.x = "CLINICAL_TRIAL_CODE",
    by.y = "row.names"
)



myLoad <- champ.filter(
    beta = myImport$beta, M = myImport$M,
    SampleCutoff = 0.15, pd = myImport$pd,
    arraytype = "EPIC",
    detP = myImport$detP
)



if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/methylome")) == F) {

    dir.create(paste0(wkdir, "/preprocess/methylome"))

}


setwd(paste0(wkdir, "/preprocess/methylome"))

# Normalizaton don't support names with spaces



myNormB <- champ.norm(
    beta = myLoad$beta, arraytype = "EPIC",
    cores = 10,
    method = "BMIQ"
)

myNormB <- myNormB[complete.cases(myNormB), ]
myNormB <- apply(myNormB, c(1,2), as.numeric)



write.table(myNormB, "myNormB.txt", sep = "\t", row.names = T, col.names = NA)



write.table(clinical.data, paste0(datadir,"/clinical_data.txt"),
    sep = "\t", row.names = F, col.names = T
)

# TRansform beta values to m

myNormM <- myNormB
myNormM[myNormM <= 0.001] <- 0.001
myNormM[myNormM >= 0.999] <- 0.999
myNormM <- log((myNormM / (1 - myNormM)), 2)

write.table(myNormM,
    "myNormM.txt",
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
)

