########################################################################
#Necesitamos el ambiente preprocess

library(umap)
library(ggplot2)
library(limma)
library(cowplot)
library(patchwork)
source("functions/Preprocess/aditional.functions.R")
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
"/preprocess/transcriptome/RMA_Expression_data_genes.tsv"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)


methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/myNormM.txt"),
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

#Cmabiamos el nombre de la calidad de


colnames(clinical.data)

colnames(clinical.data) <- gsub("pos.vs.neg.auc", "Calidad.Expresión", colnames(clinical.data))
colnames(clinical.data) <- gsub("qual.met", "Calidad.Metilación", colnames(clinical.data))



if (dir.exists(paste0(wkdir, "/Figuras_tesis")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis"))

}


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS1")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS1"))

}

pca_exp <- pca_corr(data = expression,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "Calidad.Expresión"),
outpath = paste0(wkdir, "/Figuras_tesis/FiguraS1"),
flag = "antes",
anotate1 = c(70, 150),
anotate2 = c(70, 140)
)


pca_met <- pca_corr(data = methylation,
base = clinical.data,
columns = c("METHYLATION_DAYS", "Calidad.Metilación"),
outpath = paste0(wkdir, "/Figuras_tesis/FiguraS1"),
flag = "antes",
anotate1 = c(650, 700),
anotate2 = c(650, 620)
)



expression_despues <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)


methylation_despues <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)



pca_exp_despues <- pca_corr(data = expression_despues,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "Calidad.Expresión"),
outpath = paste0(wkdir, "/Figuras_tesis/FiguraS1"),
flag = "despues",
anotate1 = c(40, 90),
anotate2 = c(40, 80)
)


pca_met_despues <- pca_corr(data = methylation_despues,
base = clinical.data,
columns = c("METHYLATION_DAYS", "Calidad.Metilación"),
outpath = paste0(wkdir, "/Figuras_tesis/FiguraS1"),
flag = "despues",
anotate1 = c(350, 800),
anotate2 = c(350, 750)
)


plt <- (pca_exp | pca_exp_despues) / (pca_met | pca_met_despues)


 ggsave(plot = plt, filename = paste0(wkdir, "/Figuras_tesis/FiguraS1/Figura_S1.pdf")
, height = 15, width = 15, device = "pdf"
,dpi =  300)
