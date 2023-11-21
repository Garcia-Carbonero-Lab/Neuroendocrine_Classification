#environment is cnv.plot
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
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



if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura10")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura10"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura10")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)



cnvdis <- read.table(paste0(wkdir, "/CNV/preprocess/CNV.peaks.dis.txt"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)


fisher <- read.table(paste0(wkdir, "/clustering/cnv/CNV_discrete_table_fisher.txt"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

fisher <- fisher[which(fisher$P.value < 0.05),]
significant <- gsub("\\.Amplification","",rownames(fisher))
significant <- gsub("\\.Deletion","",significant)

rownames(fisher) <- gsub("\\.Amplification", "", rownames(fisher))
rownames(fisher) <- gsub("\\.Deletion", "", rownames(fisher))

fisher <- fisher[order(fisher$P.value),]
cnvdis <- cnvdis[rownames(fisher),]


exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

base$Subtype <- as.factor(base$Subtype)
levels(base$Subtype) <- c("SN1", "SN2", "SN3")


base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)
levels(base$PRIMARY_TUMOR) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")


base$amp_dis <- as.factor(base$amp_dis)
levels(base$amp_dis) <- c("Sí", "No")

base$Delection_discrete <- as.factor(base$Delection_discrete)
levels(base$Delection_discrete) <- c("Sí", "No")

cnvdis[cnvdis == "Amplification"] <-  "Amplificación"
cnvdis[cnvdis == "Wild.type"] <-  "Normal"
cnvdis[cnvdis == "Deletion"] <-  "Deleción"

base <- merge(base, t(cnvdis), by = "row.names")
rownames(base) <- base$Row.names


df1 <- base %>%
    transmute(
        Subtipo = Subtype,
        Tumor.Primario = PRIMARY_TUMOR,
        Amplificación = amp_dis,
       Deleción= Delection_discrete,
      )


cnvs <- c("4p16.3",
"11q14.1",
"16q22.3",
"22q13.32",
"21q22.3",
"7p13",
"17q25.3")


df2 <- base[,c("Subtype","PRIMARY_TUMOR",cnvs)]

colnames(df2) <- c("SUBTIPO", "TUMOR.PRIMARIO", cnvs)



colors1<- list(
    "Subtipo" = c("SN1" = "gold",
    "SN2" = "cornflowerblue",
    "SN3" = "tomato1"),
    "Tumor.Primario" = c(
        "Colorrectal" = "orange",
        "Pulmón" = "blue",
        "I.Delgado" = "yellow",
        "Páncreas" = "green",
        "Estómago" = "pink"
    ),
    "Amplificación" = c("No" = "grey", "Sí" = "red"),
    "Deleción" = c("No" = "grey",
    "Sí" = "blue")   
)



colors2<- list(
    "SUBTIPO" = c("SN1" = "gold",
    "SN2" = "cornflowerblue",
    "SN3" = "tomato1"),
    "TUMOR.PRIMARIO" = c(
        "COLORRECTAL" = "orange",
        "PULMÓN" = "blue",
        "INTESTINO" = "yellow",
        "PÁNCREAS" = "green",
        "ESTÓMAGO" = "pink"
    ),
    "4p16.3" = c("Normal" = "grey", "Amplificación" = "red"),
    "11q14.1" = c("Normal" = "grey", "Delección" = "blue"),
    "16q22.3" = c("Normal" = "grey", "Delección" = "blue"),
   "22q13.32" = c("Normal" = "grey", "Delección" = "blue"),
   "21q22.3" = c("Normal" = "grey", "Delección" = "blue"),
   "7p13" = c("Normal" = "grey", "Delección" = "blue"),
   "17q25.3" = c("Normal" = "grey", "Delección" = "blue")

    
)






lgds1 <- lapply(names(colors1),function(x){

lgd <- Legend(labels = names(colors1[[x]]),
legend_gp = gpar(fill = colors1[[x]]),
title = x,
title_gp = gpar(fontsize = 80),
labels_gp = gpar(fontsize = 80),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm"),
title_position = "topcenter"
)
return(lgd)
})

lgds2 <- lapply(names(colors2),function(x){

lgd <- Legend(labels = names(colors2[[x]]),
legend_gp = gpar(fill = colors2[[x]]),
title = x,
title_gp = gpar(fontsize = 80),
labels_gp = gpar(fontsize = 80),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm")
)
return(lgd)
})

#########################################################
#MOre filteres
#Expression

deg1 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group1-Group2.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg2 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group1-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg3 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group2-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg3_filter <- rownames(deg3)[deg3$adj.P.Val < 0.1 & abs(deg3$logFC) > 1]
deg2_filter <- rownames(deg2)[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 2]
deg1_filter <- rownames(deg1)[deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 2]

genes <- unique(c(deg1_filter,deg2_filter,deg3_filter))
#METILACION

dep1 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group1-Group2.txt"),
sep = "\t",
header = T,
row.names = 1
)

dep2 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group1-Group3.txt"),
sep = "\t",
header = T,
row.names = 1
)

dep3 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group2-Group3.txt"),
sep = "\t",
header = T,
row.names = 1
)


dep1_filter <- rownames(dep1)[dep1$adj.P.Val < 0.05 & abs(dep1$logFC) > 2]
dep2_filter <- rownames(dep2)[dep2$adj.P.Val < 0.05 & abs(dep2$logFC) > 2]
dep3_filter <- rownames(dep3)[dep3$adj.P.Val < 0.05 & abs(dep3$logFC) > 2]

probes <- unique(c(dep1_filter, dep2_filter, dep3_filter))


#########################################################
source("functions/figuras/figuras.functions.R")


heatmap.molecular(exp = exp[genes,],
met = met[probes,],
df1,
outdir = outdir,
col = colors1,
size = list("width" = 50, "height" = 50),
anot_param_legend = lgds1,
flag = "Figura10.1",
min.met = NULL,
names.size = 3,
column = "Subtipo",
anot_h = c(1.5,1,1,1),
annotation_name_gp = 80,
show_row_names = F,
order = c("SN1", "SN2", "SN3"),
height_anotation = 15,
nlegend = 4)



heatmap.molecular(exp = exp[genes,],
met = met[probes,],
df2,
outdir = outdir,
col = colors2,
size = list("width" = 50, "height" = 60),
anot_param_legend = lgds2,
flag = "Figura10.2",
max.exp = 5,
min.exp = -3,
max.met = 4,
min.met = NULL,
names.size = 3,
column = "SUBTIPO",
anot_h = c(1.5,1,1,1),
annotation_name_gp = 60,
show_row_names = F,
order = c("SN1", "SN2", "SN3"))


