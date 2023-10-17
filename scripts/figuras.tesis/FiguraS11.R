#environment figuras
library(ggplot2)
library(ggpubr)
library(circlize)
library(tidyverse)
library(rstatix)
library(patchwork)
#load functions
source("functions/main.functions.R")
source("functions/figuras/figuras.functions.R")

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


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS11")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS11"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS11")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("PRIMARY_TUMOR",
"Subtype", "ImmuneScore", "TumorPurity")

base <- base[,c(features_cm)]

variables <- c("Tumor.Primario",
"Subtipo", "Score.Inmune", "Pureza.Tumoral")

colnames(base) <- c(variables)

base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")



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
    ))

#########################################################

figS10a1 <- complete.boxplot(base,
factor = "Pureza.Tumoral",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Pureza_complete",
points = T,
add.limits = F,
title = F)

figS10b1 <- wrap.boxplot(base,
factor = "Pureza.Tumoral",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Pureza_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T,
remove.y = T)

figS10a2 <- complete.boxplot(base,
factor = "Score.Inmune",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Score_Inmune",
points = T,
add.limits = F,
title = F)

figS10b2 <- wrap.boxplot(base,
factor = "Score.Inmune",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Score_Inmune_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T,
remove.y = T)


pdf(paste0(outdir, "/FiguraS11.pdf"),
width = 30,
height = 20)


figS10a1 + figS10b1 + figS10a2 + figS10b2 +

plot_layout(widths = c(1,2)) + 
plot_annotation(tag_levels = 'A') 


dev.off()





