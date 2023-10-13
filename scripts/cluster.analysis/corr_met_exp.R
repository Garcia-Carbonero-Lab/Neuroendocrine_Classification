#This script needs correlation.exp.met environment

library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

#load functions

source("functions/main.functions.R")
source("functions/cluster.analysis/corr.met.exp.functions.R")

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

if (dir.exists(paste0(wkdir, "/clustering/corr_exp_met")) == F) {

    dir.create(paste0(wkdir, "/clustering/corr_exp_met"))

}

outdir <- paste0(wkdir, "/clustering/corr_exp_met")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

exp <- read.table(paste0(wkdir,
"/clustering/deg/Expression.deg.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


met <- read.table(paste0(wkdir,
"/clustering/dep/Methylation.dep.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)

met_anot <- create.met.anot(epic)


df.corr <- corr.met.exp.unique(exp, met, met_anot)
df.corr$p.adj <- p.adjust(df.corr$p.val, method = "fdr")
df.corr <- merge(df.corr, met_anot, by.x = "CpG", by.y = "ID")
df.corr <- df.corr[order(df.corr$p.val),]

write.table(df.corr,
paste0(outdir,"/Correlation_exp_met.txt"),
sep = "\t",
row.names = T,
col.names = NA)



#corr per primary tumor

tissues <- unique(base$PRIMARY_TUMOR)

corr.tissue <- lapply(tissues, function(x) {
    df <- base[base$PRIMARY_TUMOR == x,]


    exp.df = exp[,rownames(df)]
    met.df = met[,rownames(df)]

    corr <- corr.met.exp.unique(exp.df, met.df, met_anot)
    corr$p.adj <- p.adjust(corr$p.val, method = "fdr")
    corr$tissue <- rep(x,nrow(corr))
    corr <- merge(corr, met_anot, by.x = "CpG", by.y = "ID")
    corr <- corr[order(corr$p.val),]
    return(corr)

})


corr.tissue.df <- as.data.frame(do.call(rbind, corr.tissue))

write.table(corr.tissue.df,
paste0(outdir,"/Correlation_exp_met_primary.txt"),
sep = "\t", row.names = F, col.names = T)

#We select the highest correlated cpg per gene
df.corr <- df.corr[!duplicated(df.corr$gene),]


#significant correlation
df.corr <- df.corr[df.corr$p.adj < 0.05,]
df.corr.rep <- df.corr



symbol_genes <- alias2SymbolTable(rownames(exp))
rownames(exp)[is.na(symbol_genes) == F] <- symbol_genes[is.na(symbol_genes) == F]


met.corr <- met[df.corr$CpG,]
rownames(met.corr) <- df.corr$gene

# Correlation plots of genes that correlates expression and methylation


genes <- c("PTPRN", "PTPRN2", "CHGA", "PIGR", "CHGB")




correlation.plot(exp = exp,
met = met,
anot = df.corr.rep,
base = base,
genes = genes,
outdir = outdir
)


#