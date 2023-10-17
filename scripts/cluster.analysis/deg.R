#This script needs the environment deg


library(limma)
library(ggplot2)
library(rsample)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)

#load functions

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")

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

if (dir.exists(paste0(wkdir, "/clustering/deg")) == F) {

    dir.create(paste0(wkdir, "/clustering/deg"))

}


outdir <- paste0(wkdir, "/clustering/deg")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

base$PRIMARY_TUMOR <- as.factor(gsub(" ", "_", base$PRIMARY_TUMOR))


# open expression data
exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


# Obtain DEGS
deg.function(exp,
base,
group = "Subtype",
cov = "PRIMARY_TUMOR",
outdir = paste0(wkdir, "/clustering/deg"),
flag = "deg_cov_primary")


# write DEGS
deg1 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group1-Group2.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg2 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group1-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg3 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group2-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)


#Volcano plots

v3 <- EnhancedVolcano(deg3,
lab = rownames(deg3),
x = "logFC",
y = "adj.P.Val",
title = "Subtype2 vs Subtype3",
FCcutoff = 0.2,
pCutoff = 0.1,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(deg1[["logFC"]], na.rm = TRUE) - 0.5, max(deg1[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg1[["adj.P.Val"]]), na.rm = TRUE) + 1))

ggsave(v3,
filename = paste0(wkdir, "/clustering/deg/Volcano_S2_vs_S3.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)


v2 <- EnhancedVolcano(deg2,
lab = rownames(deg2),
x = "logFC",
y = "adj.P.Val",
title = "Subtype1 vs Subtype2",
FCcutoff = 0.5,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(deg2[["logFC"]], na.rm = TRUE) - 0.5, max(deg2[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg2[["adj.P.Val"]]), na.rm = TRUE) + 1))

ggsave(v2,
filename = paste0(wkdir, "/clustering/deg/Volcano_S1_vs_S2.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

v1 <- EnhancedVolcano(deg1,
lab = rownames(deg1),
x = "logFC",
y = "adj.P.Val",
title = "Subtype1 vs Subtype2",
FCcutoff = 0.5,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
xlim = c(min(deg3[["logFC"]], na.rm = TRUE) - 0.5, max(deg3[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg3[["adj.P.Val"]]), na.rm = TRUE) + 1))


ggsave(v1,
filename = paste0(wkdir, "/clustering/deg/Volcano_S1_vs_S2.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

#Table filter

deg3_filter <- deg3[deg3$adj.P.Val < 0.1 & abs(deg3$logFC) > 0.2,]
deg2_filter <- deg2[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 0.5,]
deg1_filter <- deg1[deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 0.5,]


# Write filtered deg table
write.table(deg1_filter,
paste0(outdir, "/Table_SN1_vs_SN2_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(deg2_filter,
paste0(outdir, "/Table_SN1_vs_SN3_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(deg3_filter,
paste0(outdir, "/Table_SN2_vs_SN3_filter.txt"),
sep = "\t",
row.names = T,
col.names = NA)


# genes selected by filters
degs3 <- rownames(deg3)[deg3$adj.P.Val < 0.1 & abs(deg3$logFC) > 0.2]

degs2 <- rownames(deg2)[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 0.5]

degs1 <- rownames(deg1)[deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 0.5]



# Plot genes selected by filters
degs <- unique(c(degs1,degs2,degs3))

deg.plot(as.data.frame(t(exp)),
base,
"PRIMARY_TUMOR",
"Subtype",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = paste0(wkdir, "/clustering/deg"),
flag = "degs_primary",
genes = degs
)



exp.degs <- exp[ unique(c(degs1,degs2,degs3)),]

write.table(exp.degs,
paste0(wkdir,
"/clustering/deg/Expression.deg.txt"),
sep = "\t",
row.names= T,
col.names= NA)

