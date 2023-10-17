#This script needs the environment deg


library(limma)
library(ggplot2)
library(rsample)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(patchwork)

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/figuras/figuras.functions.R")


# we read config file
config <- read.csv("config/config.template.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS16")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS16"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS16")


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)
#############################################################

features_cm <- c("PRIMARY_TUMOR",
"Subtype")

base <- base[,c(features_cm)]

variables <- c("Tumor.Primario",
"Subtipo")

colnames(base) <- c(variables)

base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")

l_colors <- list(
    "Subtipo" = c("SN1" = "gold",
    "SN2" = "cornflowerblue",
    "SN3" = "tomato1"),
    "Tumor.Primario" = c(
        "Colorrectal" = "orange",
        "Pulmón" = "blue",
        "I.Delgado" = "yellow",
        "Páncreas" = "green",
        "Estómago" = "pink"
    ))


rownames(exp) <- gsub("-","_",rownames(exp))

genes <- c("CHGA","CHGB",
"PTPRN","BEX1",
"GSTA1", "FBP1", "FABP1","PIGR" ,"SPINK1", "GJD2",
"GRB10", "HEPACAM2")


bxp <- deg.boxplot(as.data.frame(t(exp)),
base,
"Tumor.Primario",
"Subtipo",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "degs_3_importantes_primary",
genes = genes
)

pdf(paste0(outdir, "/FiguraS16.pdf"),
width = 20,
height = 20)

bxp[[1]] + bxp[[2]] + 
bxp[[3]] + bxp[[4]] +
bxp[[5]] + bxp[[6]]  +
bxp[[7]] + bxp[[8]] + 
bxp[[9]] + bxp[[10]] +
bxp[[11]] + bxp[[12]]  +
plot_layout(ncol = 2, guide = "collect")

dev.off()



