
library(MOFA2)
library(survival)
library(survminer)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
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


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

factor.rm <- c("Factor3", "Factor6")



if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura8")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura8"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura8")

base <- base[,c("Subtype",
"PRIMARY_TUMOR",
"GENDER",
"GRADE",
"CIMP.status",
"STAGE_IV",
"X5HIAA_HIGH")]


colnames(base) <- c("Subtipo",
"Tumor.Primario",
"Género",
"Grado",
"CIMP",
"Metástasis",
"5HIAA")

base$Género <- as.factor(base$Género)
levels(base$Género) <- c("Mujer", "Hombre")

base[,"5HIAA"] <- as.factor(base[,"5HIAA"])
levels(base[,"5HIAA"]) <- c("Bajo", "Alto")

base$CIMP <- as.factor(base$CIMP)
levels(base$CIMP) <- c("Alto", "Bajo")

base$Metástasis <- as.factor(base$Metástasis)
levels(base$Metástasis) <- c("No", "Sí")

base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón", 
"Páncreas",
"I.Delgado")


base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")


mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))




base<- merge(mod@samples_metadata,
base, by.x = "sample", by.y = "row.names")

Z <- obtain_Z(mod)



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
    ),
    "Grado" = c("G1" = "green",
    "G2" = "yellow",
    "G3" = "red"),
    "CIMP" = c("Alto" = "red",
    "Bajo"= "blue"),
    "Metástasis" = c("No" = "forestgreen", "Sí" = "red4"),
    "5HIAA" = c("Alto"= "brown2",
    "Bajo" = "steelblue1")
)

lgds <- lapply(names(l_colors),function(x){

lgd <- Legend(labels = names(l_colors[[x]]),
legend_gp = gpar(fill = l_colors[[x]]),
title = x,
title_gp = gpar(fontsize = 120),
labels_gp = gpar(fontsize = 120),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm"),
title_position = "topcenter"

)
return(lgd)
})


source("functions/figuras/figuras.functions.R")

heatmap_classes(Z = Z,
base = base,
outdir = outdir,
col = l_colors,
size = list("width" = 60, "height" = 70),
factors.rm = factor.rm,
anot_param_legend = lgds,
anot_h = c(1.5,1,1,1,1,1),
column = "Subtipo",
flag = "Figura8")


heatmap_classes(Z = Z,
base = base,
outdir = outdir,
col = l_colors,
size = list("width" = 40, "height" = 40),
factors.rm = factor.rm,
anot_param_legend = lgds,
anot_h = c(1.5,1,1,1,1,1),
column = "Subtipo",
flag = "Figura8v2")
###################################################
mod <- load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

base <- base[,c("Subtype",
"PRIMARY_TUMOR",
"OS.time",
"EXITUS")]

colnames(base) <- c("Subtipo", 
"Tumor.Primario",
"OS.time",
"EXITUS")


base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")


source("functions/figuras/figuras.functions.R")

fig8c <- Kapplan_Meyer(cut_survival(base, 150),
"Subtipo",
"OS.time",
"EXITUS",
NULL,
outdir,
"Figura8",
c("gold",
"cornflowerblue",
"tomato1"),
width= 8,
height=8) 

#######################################333
#Figura 8C

base<- merge(mod@samples_metadata,
base, by.x = "sample", by.y = "row.names")

mod@samples_metadata <- base


pt <- plot_factors(mod, c(1,2), color_by = "Subtipo")+
            
    guides(colour = guide_legend(override.aes = list(size=20))) +

            theme(
                text = element_text(size = 25),
                axis.text.x = element_text(vjust = 0.5, hjust = 1),
                axis.ticks.length.y = unit(5,"pt"),
                axis.ticks.length.x = unit(5,"pt")

            ) +

scale_fill_manual(values = c("gold","cornflowerblue","tomato1"))




pdf(paste0(outdir, "/Figura8C.pdf"), width = 10, height = 5)
pt
dev.off()


pdf(paste0(outdir, "/Figura8Cv2.pdf"), width = 5,
height = 5)
pt
dev.off()

