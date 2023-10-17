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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura7")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura7"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura7")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("BIOPSY_METASTASIS_LOCATION")




base <- base[,c(features_cm, "PRIMARY_TUMOR")]

colnames(base) <- c("Localización.Biopsia",
"Tumor.Primario")

base$Localización.Biopsia[is.na(base$Localización.Biopsia)] <- "Tumor.Primario"

base$Localización.Biopsia <- as.factor(base$Localización.Biopsia)
levels(base$Localización.Biopsia) <- c("Hígado",
"Pulmón",
"Nódulo.Linfático",
"Ovario",
"Parótida",
"Peritoneo",
"Piel",
"Tumor.Primario")


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
    "Localización.Biopsia" = c("Hígado" = "brown",
    "Pulmón" = "blue",
    "Nódulo.Linfático" = "lightskyblue1",
    "Ovario" = "lightpink",
    "Parótida" = "magenta",
    "Peritoneo"= "lightsalmon",
    "Piel" ="gold",
    "Tumor.Primario" = "grey"
    )
)


fig7 <- scatterplot(df,
factor = "Factor.6",
column = "Localización.Biopsia",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor6_Biopsia",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
legend = T)

pdf(paste0(outdir, "/Figura7.pdf"), width = 20, height = 10)
fig7
dev.off()

