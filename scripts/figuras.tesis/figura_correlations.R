#environment figuras

library(MOFA2)
library(circlize)
library(ComplexHeatmap)
library(patchwork)
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura_correlación")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura_correlación"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura_correlación")

mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

factors_names(mod) <- gsub("Factor", "Factor.", factors_names(mod))


views_names(mod) <- c("Expresión", "Metilación")

exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)

corr_pca_lf(exp,  mod, outdir, "Expression")
corr_pca_lf(met,  mod, outdir, "Methylation")
