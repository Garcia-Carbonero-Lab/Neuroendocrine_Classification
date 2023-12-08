library(survival)
library(survminer)
library(patchwork)
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura_supervivencia")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura_supervivencia"))

}

outdir <- paste0(wkdir, "/Figuras_tesis/Figura_supervivencia")


base <- base[,c(
"PRIMARY_TUMOR",
"OS.time",
"EXITUS")]

colnames(base) <- c( 
"Tumor.Primario",
"OS.time",
"EXITUS")


base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal", "Estómago", "Pulmón", "Páncreas", "I.Delgado")
base$Tumor.Primario <- factor(base$Tumor.Primario, levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"))


fig<- Kapplan_Meyer(cut_survival(base, 150),
"Tumor.Primario",
"OS.time",
"EXITUS",
NULL,
outdir,
"Figura8",
c("orange",
"pink",
"yellow",
"green",
"blue"),
width= 8,
height=8,
risk.table = T) 
