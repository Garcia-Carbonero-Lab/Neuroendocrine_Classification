
#This script needs the environment clinical

library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(tidyverse)
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


if (dir.exists(paste0(wkdir, "/clustering/clinical")) == F) {

    dir.create(paste0(wkdir, "/clustering/clinical"))

}


outdir <- paste0(wkdir, "/clustering/clinical")



base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

base$Subtype <- as.character(base$Subtype)
base$Subtype <- as.factor(base$Subtype)


features <- c("TRANSCRIPTOME_BATCH","AGE.median",
    "BIOPSY_SITES", "BIOPSY_METASTASIS_LOCATION"
             , "GENDER", "PRIMARY_TUMOR"
             , "GRADE", "FUNCTIONAL_dichotomic"
             , "BIOPSY_LOCATION", "X5HIAA_HIGH", "X5HIAA_HIGH2"
             , "CIMP.status", "STAGE_TNM", "STAGE_IV")


base$BIOPSY_METASTASIS_LOCATION[is.na(base$BIOPSY_METASTASIS_LOCATION)] <- "PRIMARY.TUMOR"


continous_features <- c("AGE","Ki.67...",
"StromalScore",
"ImmuneScore",
"ESTIMATEScore",
"TumorPurity")

#Create table indicating the distribution of each features among subtypes
table <- clinical.table(data = base,
features = features,
group = "Subtype")

#Create table indicating the distribution of each features among subtypes
#taking into account the NA level
table_na <- clinical.table(data = base,
features = features,
group = "Subtype",
na = T)


#Create table indicating the distribution of each features among primary tumors


table_primary <- clinical.table(data = base,
features = c(features, "Subtype"),
group = "PRIMARY_TUMOR")

#Create table indicating the distribution of each features among primary tumors
#taking into account the NA level
table_primary_na <- clinical.table(data = base,
features = c(features, "Subtype"),
group = "PRIMARY_TUMOR",
na = T)


write.table(table,
paste0(outdir, "/Clinical.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)



write.table(table_primary,
paste0(outdir, "/Clinical_primary.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)


write.table(table_primary_na,
paste0(outdir, "/Clinical_primary_na.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)

write.table(table_na,
paste0(outdir, "/Clinical_subtypes_na.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)


pm <- unique(base$PRIMARY_TUMOR)
for (p in pm){

df <- base[base$PRIMARY_TUMOR == p,]
table <- clinical.table(data = df, features = features, group = "Subtype")

write.table(table,
paste0(outdir, "/Clinical_",p,".txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F)

} 



l_colors <- list(
    "Subtype" = c("1" = "gold",
    "2" = "cornflowerblue",
    "3" = "tomato1"),
    "PRIMARY_TUMOR" = c(
        "COLORECTAL" = "orange",
        "LUNG" = "blue",
        "SMALL INTESTINE" = "yellow",
        "PANCREAS" = "green",
        "GASTRIC" = "pink"
    ),
    "BIOPSY_LOCATION" = c(
        "COLORECTAL" = "orange",
        "LUNG" = "blue",
        "SMALL INTESTINE" = "yellow",
        "PANCREAS" = "green",
        "GASTRIC" = "pink",
        "PERITONEUM" = "lightsalmon",
        "LYMPH NODE" = "lightskyblue1",
        "LIVER" = "brown",
        "LUNG" = "blue",
        "OVARY" = "lightpink",
        "SKIN" = "gold",
        "PAROTID GLAND" = "magenta"
    ),
    "AGE.median" = c("<= 56" = "seagreen", "> 56" = "chocolate1"), 
    "GENDER" = c("FEMALE" = "lightpink", "MALE" = "lightskyblue1"),
    "GRADE" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "STAGE_IV" = c("YES" = "red4",  "NO" = "forestgreen"
    ),
    "X5HIAA_HIGH" = c("YES" = "steelblue1", "NO" = "brown2"),
    "CIMP.status" = c("Low" = "blue", "Inter" = "green", "High" = "red")
)

#############################################
#Boxplot Continous


continous_features <- c(
"StromalScore",
"ImmuneScore",
"ESTIMATEScore",
"TumorPurity")

#Calculate pvalues of kruskal-Wallis test
kruskal <- kruskal(base,
"Subtype",
continous_features,
outdir = outdir,
flag = "ESTIMATE")

#Create the boxplots and study differences using Wilcoxon rank test
boxplots(base,
group.x = "PRIMARY_TUMOR",
group.fill = "Subtype",
color.fill = c("gold", "cornflowerblue", "tomato1"),
features = continous_features,
outdir = outdir,
flag = "ESTIMATE_primary")

boxplots(base,
group.x = "Subtype",
group.fill = "Subtype",
color.fill = c("gold", "cornflowerblue", "tomato1"),
features = continous_features,
outdir = outdir,
flag = "ESTIMATE_complete")


