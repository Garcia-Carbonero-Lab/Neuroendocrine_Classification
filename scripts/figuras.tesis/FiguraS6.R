#environment figuras

library(MOFA2)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(survminer)
library(survival)
library(tidyverse)
library(rstatix)
library(patchwork)
library(tidyHeatmap)
library(figpatch)
#load functions
source("functions/main.functions.R")
source("functions/figuras/figuras.functions.R")

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



mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

factors_names(mod) <- gsub("Factor", "Factor.", factors_names(mod))


views_names(mod) <- c("Expresión", "Metilación")

groups_names(mod) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS6")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS6"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS6")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("StromalScore",
"ImmuneScore",
"ESTIMATEScore",
"TumorPurity")

base <- base[,c(features_cm, "PRIMARY_TUMOR")]

variables <- c("Score.Estromal",
"Score.Inmune",
"Score.ESTIMATE",
"Pureza.Tumoral")

colnames(base) <- c(variables, "Tumor.Primario")


df <- base[mod@samples_metadata$sample, ]

Z <- obtain_Z(mod)

df <- cbind(df, Z)

mod@samples_metadata <- merge(mod@samples_metadata,
    df, by.x = 'sample', by.y = 'row.names')

################################################################################
#Figura 2A

columns <-  c("Score.Estromal",
"Score.Inmune",
"Score.ESTIMATE",
"Pureza.Tumoral")

fs6<- heatmap.correlacion(mod,
columns,
outdir,
flag = "correlation.estimate")

