#This script needs MOFA.nlcor environment
library(MOFA2)
library(nlcor)
library(pheatmap)
library(ggplot2)


source("functions/MOFA/MOFA.nonlinear.R")
source("functions/MOFA/adictional.functions.R")

config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]
#Open clinical data and model        

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

#CHANGE COLUMNS CLASS
base$X5.hydroxyindoleacetic_acid..mg.24H. <- as.numeric(base$X5.hydroxyindoleacetic_acid..mg.24H.)


# We indicate the out folder
outdir <- paste0(wkdir, "/MOFA/analysis_results/")




if (dir.exists(paste0(outdir, "/continous.variables/")) == F) {

    dir.create(paste0(outdir, "/continous.variables/"))

}


# Anotate the model
mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

df <- base[mod@samples_metadata$sample, ]

Z <- obtain_Z(mod)

df <- cbind(df, Z)

mod@samples_metadata <- merge(mod@samples_metadata,
    df, by.x = 'sample', by.y = 'row.names')



continous_noise <- c(
"pos.vs.neg.auc",
"qual.met"
)



continous.features.nonlinear(mod = mod,
p.th = 0.05,
columns = continous_noise,
outdir = outdir,
flag = "noise_non_linear",
select = F)




continous_features <- c(
"Ki.67...",
"X5.hydroxyindoleacetic_acid..mg.24H.",
"AGE", "StromalScore", "ImmuneScore", "ESTIMATEScore","TumorPurity"
)




continous.features.nonlinear(mod = mod,
p.th = 0.05,
columns = continous_features,
outdir = outdir,
flag = "normal_features_non_linear")
