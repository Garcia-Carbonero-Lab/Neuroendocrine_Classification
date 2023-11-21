#environment figuras

library(MOFA2)
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
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]



mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

factors_names(mod) <- gsub("Factor", "Factor.", factors_names(mod))


views_names(mod) <- c("Expresión", "Metilación")

groups_names(mod) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS2")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS2"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS2")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("GENDER",
"GRADE",
"STAGE_IV",
"X5HIAA_HIGH",
"CIMP.status",
"OS.time",
"EXITUS")

base <- base[,c(features_cm, "PRIMARY_TUMOR")]

variables <- c("Género",
"Grado",
"Metástasis",
"X5HIAA",
"CIMP",
"Tiempo",
"Exitus")

colnames(base) <- c(variables, "Tumor.Primario")

base$Género <- as.factor(base$Género)
levels(base$Género) <- c("Mujer", "Hombre")

base$Metástasis <- as.factor(base$Metástasis)
levels(base$Metástasis) <- c("No", "Sí")

base[,"X5HIAA"] <- as.factor(base[,"X5HIAA"])
levels(base[,"X5HIAA"]) <- c("Bajo", "Alto")
base[,"X5HIAA"] <- factor(base[,"X5HIAA"], levels = c("Alto", "Bajo"))

base$CIMP <- as.factor(base$CIMP)
levels(base$CIMP) <- c("Alto", "Bajo")


base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

df <- base[mod@samples_metadata$sample, ]

Z <- obtain_Z(mod)

df <- cbind(df, Z)

mod@samples_metadata <- merge(mod@samples_metadata,
    df, by.x = 'sample', by.y = 'row.names')

###################################################################

source("functions/figuras/figuras.functions.R")

colors_fcm <- list(
    "Género" = c("Mujer" = "lightpink", "Hombre" = "lightskyblue1"),
    "Grado" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "Metástasis" = c("No" = "forestgreen", "Sí" = "red4"),
    "X5HIAA" = c("Alto" = "brown2", "Bajo" = "steelblue1"),
    "CIMP" = c("Alto" = "red","Bajo" = "blue")

)


figs2a1 <- complete.boxplot(df,
factor = "Factor.1",
column = "Grado",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_Grado_complete.pdf")


figs2a2 <- complete.boxplot(df,
factor = "Factor.1",
column = "X5HIAA",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_5HIAA_complete.pdf",
remove.y = T)

figs2a3 <- complete.boxplot(df,
factor = "Factor.1",
column = "Metástasis",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_estadio_complete.pdf",
remove.y = T)



figs2b1 <- wrap.boxplot(df,
factor = "Factor.1",
column = "Grado",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_Gradowrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)


figs2b2 <- wrap.boxplot(df,
factor = "Factor.1",
column = "X5HIAA",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_X5HIAA_wrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)


figs2b3 <- wrap.boxplot(df,
factor = "Factor.1",
column = "Metástasis",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_Estadio_wrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)

figs2b12 <- scatterplot(df,
factor = "Factor.1",
column = "Grado",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_Grado_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)


figs2b22 <- scatterplot(df,
factor = "Factor.1",
column = "X5HIAA",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_5HIAA_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)



figs2b32 <- scatterplot(df,
factor = "Factor.1",
column = "Metástasis",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_Estadio_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)



pdf(file = paste0(outdir, "/FiguraS2.pdf"), height = 30, width = 50)
figs2a1 <- figs2a1 + labs(tag = "A")
figs2b1 <- figs2b1 + labs(tag = "B")

figs2a2 <- figs2a2 + labs(tag = " ")
figs2b2 <- figs2b2 + labs(tag = " ")


figs2a3 <- figs2a3 + labs(tag = " ")
figs2b3 <- figs2b3 + labs(tag = " ")

patch1 <- figs2a1 / (figs2b1 / figs2b12 + plot_layout(guides = "auto")) +
plot_layout(guides = "auto", heights = c(1, 2))

patch2 <- figs2a2 / (figs2b2 / figs2b22 + plot_layout(guides = "auto")) +
plot_layout(guides = "auto", heights = c(1, 2))


patch3 <- figs2a3 / (figs2b3 / figs2b32 + plot_layout(guides = "auto")) +
plot_layout(guides = "auto", heights = c(1, 2))

(patch1 | patch2 | patch3)


dev.off()


