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
config <- read.csv("config/config.template.tsv",
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura3")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura3"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura3")

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




fig3a <- complete.boxplot(df,
factor = "Factor.1",
column = "CIMP",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_FMIC_complete.pdf")


fig3b1 <- wrap.boxplot(df,
factor = "Factor.1",
column = "CIMP",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_FMIC_wrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)


fig3b2 <- scatterplot(df,
factor = "Factor.1",
column = "CIMP",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor1_FMIC_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)



pdf(file = paste0(outdir, "/Figura3.pdf"), height = 30, width = 20)
fig3a <- fig3a + labs(tag = "A")
fig3b1 <- fig3b1 + labs(tag = "B")
fig3b2 <- fig3b2 + labs(tag = " ")

fig3a + (fig3b1 + fig3b2 + plot_layout(ncol = 1, guides = "collect")) +
plot_layout(ncol = 1, guides = "auto", heights = c(1,2))

dev.off()










