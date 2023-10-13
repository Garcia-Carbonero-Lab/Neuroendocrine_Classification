
#This script needs environment MOFA

library(MOFA2)
library(ComplexHeatmap)
library(circlize)


#load functions

source("functions/main.functions.R")
source("functions/clustering/clustering.functions.R")
source("functions/MOFA/adictional.functions.R")

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

#load_models

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)



factor.rm <- c("Factor3", "Factor6")



outdir <- paste0(wkdir, "/clustering/clustering_process")


base_general <- base[,c("Subtype","PRIMARY_TUMOR", "GRADE", "STAGE_TNM")]
colnames(base_general) <- c("SUBTYPE", "PRIMARY TUMOR", "GRADE", "STAGE")

base_noise <- base[,c("Subtype",
"PRIMARY_TUMOR",
"GRADE",
"STAGE_TNM",
"METHYLOMA_BATCH",
"METHYLATION_DAYS",
"CIMP.status")]

colnames(base_noise) <- c("SUBTYPE",
"PRIMARY TUMOR",
"GRADE",
"STAGE",
"METHYLOMA_BATCH",
"METHYLATION_DAYS",
"CIMP.status")


l_colors <- list(
    "SUBTYPE" = c("1" = "gold",
    "2" = "cornflowerblue",
    "3" = "tomato1"),
    "PRIMARY TUMOR" = c(
        "COLORECTAL" = "orange",
        "LUNG" = "blue",
        "SMALL INTESTINE" = "yellow",
        "PANCREAS" = "green",
        "GASTRIC" = "pink"
    ),
    "GRADE" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "STAGE" = c(
        "I" = "green", "II" = "yellow",
        "III" = "orange", "IV" = "red"
    )
)




l_colors_noise <- list(
    "SUBTYPE" = c("1" = "gold",
    "2" = "cornflowerblue",
    "3" = "tomato1"),
    "PRIMARY TUMOR" = c(
        "COLORECTAL" = "orange",
        "LUNG" = "blue",
        "SMALL INTESTINE" = "yellow",
        "PANCREAS" = "green",
        "GASTRIC" = "pink"
    ),
    "GRADE" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "STAGE" = c(
        "I" = "green", "II" = "yellow",
        "III" = "orange", "IV" = "red"
    ),
    "METHYLOMA_BATCH" = c("BATCH1" = "blue",
    "BATCH2" = "palegreen",
    "BATCH3" = "darkolivegreen",
    "BATCH4" = "purple"),
    "CIMP.status" = c("Low" = "blue", "Inter" = "green", "High" = "red"),
    "METHYLATION_DAYS" = colorRamp2(
            c(
                min(base_noise$METHYLATION_DAYS, na.rm = T),
                max(base_noise$METHYLATION_DAYS, na.rm = T)),
            c("white", "green")))



anot_param_legend <- list(
    "SUBTYPE" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    
    ),
    "PRIMARY TUMOR" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    ),
    
    "GRADE" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 1
    ),
    "STAGE" = list(
       title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    )
)

anot_param_legend_noise <- list(
    "SUBTYPE" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    
    ),
    "PRIMARY TUMOR" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    ),
    
    "GRADE" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 1
    ),
    "STAGE" = list(
       title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2),
    "METHYLOMA_BATCH" =  list(
       title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2),
        "CIMP.status" = list(
       title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2),
        
        "METHYLATION_DAYS" = list(title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                legend_height = unit(8, "cm"),
                grid_width = unit(2, "cm"))
)

mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))




base_general <- merge(mod@samples_metadata,
base_general, by.x = "sample", by.y = "row.names")


base_noise <- merge(mod@samples_metadata,
base_noise, by.x = "sample", by.y = "row.names")


base_general$SUBTYPE <- as.factor(as.character(base_general$SUBTYPE))
base_noise$SUBTYPE <- as.factor(as.character(base_noise$SUBTYPE))


Z <- obtain_Z(mod)



heatmap_classes(Z = Z, base = base_general,
outdir = outdir,
col = l_colors,
size = list("width" = 50, "height" = 45),
factors.rm = factor.rm,
anot_param_legend = anot_param_legend,
anot_h = c(40,20,20,20),
column = "SUBTYPE",
flag = "general")



heatmap_classes(Z = Z, base = base_noise,
outdir = outdir,
col = l_colors_noise,
size = list("width" = 50, "height" = 50),
factors.rm = factor.rm,
anot_param_legend = anot_param_legend_noise,
anot_h = c(40,20,20,20,20,20,20),
column = "SUBTYPE",
flag = "noise")
