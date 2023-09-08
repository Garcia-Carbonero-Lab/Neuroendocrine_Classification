#For this script we need the enviroment arrays_m

library(limma)
library(makecdfenv)
library(affy)
library(vsn)



#Open config file
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]


#load functions
source("functions/Preprocess/expression.process.R")


targetfile <- paste0(datadir, "/transcriptome/raw/targets.txt")
inputpath <- paste0(datadir, "/transcriptome/raw")



if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/transcriptome")) == F) {

    dir.create(paste0(wkdir, "/preprocess/transcriptome"))

}

outputpath <- paste0(wkdir, "/preprocess/transcriptome")

clariomshumancdf <- make.cdf.env("Clariom_S_Human.r1.Gene.CDF", inputpath)

Norm_exp(targetfile = targetfile
         , outputpath = outputpath
         , flag = "RMA"
         , inpath = inputpath)
