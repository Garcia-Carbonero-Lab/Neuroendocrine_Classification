# For use this script we need the environment preprocess

# we open libraries
library(umap)
library(ggplot2)
library(limma)

# open functions

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

#######################################################
# Check correlation between quality and data


if (dir.exists(paste0(wkdir, "/preprocess")) == F) {

    dir.create(paste0(wkdir, "/preprocess"))

}

if (dir.exists(paste0(wkdir, "/preprocess/quality.plots")) == F) {

    dir.create(paste0(wkdir, "/preprocess/quality.plots"))
    dir.create(paste0(wkdir, "/preprocess/quality.plots/transcriptome"))
    dir.create(paste0(wkdir, "/preprocess/quality.plots/methylome"))

}



pca_exp <- pca_corr(data = expression,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "pos.vs.neg.auc"),
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "before.expression"
)


umap_exp <- umap_corr(data = expression,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "pos.vs.neg.auc"),
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "before.expression",
anotate1 = c(3,4),
anotate2 = c(3,3.7)


)

write.table(pca_exp,
paste0(wkdir, "/preprocess/quality.plots/transcriptome/Correlation.PCA.features.before.expression.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

write.table(umap_exp,
paste0(wkdir, "/preprocess/quality.plots/transcriptome/Correlation.UMAP.features.before.expression.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)


pca_met <- pca_corr(data = methylation,
base = clinical.data,
columns = c("METHYLATION_DAYS", "qual.met"),
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "before.methylation",
anotate1 = c(650, 450),
anotate2 = c(650, 420)
)


umap_met <- umap_corr(data = methylation,
base = clinical.data,
columns = c("METHYLATION_DAYS", "qual.met"),
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "before.methylation",
anotate1 = c(2, 4),
anotate2 = c(2, 3.7)
)

write.table(pca_met,
paste0(wkdir, "/preprocess/quality.plots/methylome/Correlation.PCA.features.before.methylation.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

write.table(umap_met,
paste0(wkdir, "/preprocess/quality.plots/methylome/Correlation.UMAP.features.before.methylation.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

#########################################################
# eliminate PCA1 associated to quality

pcaexp <- generate_pca(expression)

pcaexp <- pcaexp[colnames(expression),,drop = F]

expression_qual <- removeBatchEffect(x = expression,
covariates = pcaexp[, "PC1"])


write.table(expression_qual,
paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = "\t",
row.names = T,
col.names = NA
)



pcamet <- generate_pca(methylation)

pcamet <- pcamet[colnames(methylation),,drop = F]

methylation_qual <- removeBatchEffect(x = methylation,
covariates = pcamet[, "PC1"])

write.table(methylation_qual,
paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = "\t",
row.names = T,
col.names = NA
)
#########################################################
#

pca_exp_qual <- pca_corr(data = expression_qual,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "pos.vs.neg.auc"),
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "after.expression",
anotate1 = c(50, 90),
anotate2 = c(50, 80)

)


umap_exp_qual <- umap_corr(data = expression_qual,
base = clinical.data,
columns = c("TRANSCRIPTOME_DAYS", "pos.vs.neg.auc"),
outpath = paste0(wkdir, "/preprocess/quality.plots/transcriptome"),
flag = "after.expression",
anotate1 = c(1.8, 4),
anotate2 = c(1.8, 3.7)
)

write.table(pca_exp_qual,
paste0(wkdir, "/preprocess/quality.plots/transcriptome/Correlation.PCA.features.expression.after.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

write.table(umap_exp_qual,
paste0(wkdir, "/preprocess/quality.plots/transcriptome/Correlation.UMAP.features.expression.after.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)
########################################################

pca_met_qual <- pca_corr(data = methylation_qual,
base = clinical.data,
columns = c("METHYLATION_DAYS", "qual.met"),
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "after.methylation",
anotate1 = c(350, 500),
anotate2 = c(350, 470)
)


umap_met_qual <- umap_corr(data = methylation_qual,
base = clinical.data,
columns = c("METHYLATION_DAYS", "qual.met"),
outpath = paste0(wkdir, "/preprocess/quality.plots/methylome"),
flag = "after.methylation",
anotate1 = c(2.5, 3),
anotate2 = c(2.5, 2.77)
)

write.table(pca_met_qual,
paste0(wkdir, "/preprocess/quality.plots/methylome/Correlation.PCA.features.after.methylation.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

write.table(umap_met_qual,
paste0(wkdir, "/preprocess/quality.plots/methylome/Correlation.UMAP.features.after.methylation.tsv"), 
sep = "\t",
row.names = T
, col.names = NA
)

