#This script needs the environment deg

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)

source("functions/main.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/cluster.analysis/clinical.functions.R")
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


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)




#We create output folder
if (dir.exists(paste0(wkdir, "/clustering/cnv")) == F) {

    dir.create(paste0(wkdir, "/clustering/cnv"))

}

outdir <- paste0(wkdir, "/clustering/cnv")

######################################################################
#CNV discrete
cnvdis <- read.table(paste0(wkdir, "/CNV/preprocess/CNV.peaks.dis.txt"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

cnvdis <- as.data.frame(t(cnvdis))
cnvdis[rownames(base),"Subtype"] <- base$Subtype
cnvdis[rownames(base),"Tumor.Primario"] <- base$PRIMARY_TUMOR

#study differences in CNV among subtypes
table <- clinical.table(data = cnvdis, 
features = colnames(cnvdis)[1:ncol(cnvdis) -1],
group = "Subtype")


table_general <- clinical.table(data = base, 
features = c("amp_dis", "del_dis", "cnv_dis"),
group = "Subtype")

table_general2 <- clinical.table(data = base, 
features = c("amp_dis", "Delection_discrete"),
group = "Subtype")


table_primary <- clinical.table(data = cnvdis, 
features = colnames(cnvdis)[1:ncol(cnvdis) -1],
group = "Tumor.Primario")

write.table(table,
paste0(outdir, "/CNV_discrete_table_fisher.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)

write.table(table_general,
paste0(outdir, "/CNV_discrete_table_fisher_general.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)

write.table(table_general2,
paste0(outdir, "/CNV_discrete_table_fisher_general2.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)

write.table(table_primary,
paste0(outdir, "/CNV_discrete_table_fisher_primary.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)


#Calculate differences by primary tumor

cnvdis[rownames(base), "PRIMARY.TUMOR"] <- base$PRIMARY_TUMOR


pm <- unique(cnvdis$PRIMARY.TUMOR)

for (p in pm){

if (dir.exists(paste0(outdir, "/", p)) == F) {

    dir.create(paste0(outdir, "/", p))

}

df <- cnvdis[cnvdis$PRIMARY.TUMOR == p,]
df_general <- base[base$PRIMARY_TUMOR == p,]


table <- clinical.table(data = df, 
features = colnames(cnvdis)[1:ncol(cnvdis) -1],
group = "Subtype")


table_general <- clinical.table(data = df_general, 
features = c("amp_dis", "del_dis", "cnv_dis"),
group = "Subtype")


write.table(table_general,
paste0(outdir, "/", p, "/CNV_discrete_table_fisher_general_",p,".txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)


write.table(table,
paste0(outdir, "/", p, "/CNV_discrete_table_fisherl_",p,".txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)


}
