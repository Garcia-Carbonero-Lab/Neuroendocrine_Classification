#This script need the environment survival
library(survival)
library(survminer)
library(tidyverse)
library(aod)
library(forestmodel)
library(patchwork)
source("functions/cluster.analysis/survival.functions.R")

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


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS9")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS9"))

}

outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS9")


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

base$GENDER <- as.factor(base$GENDER)
levels(base$GENDER) <- c("Mujer", "Hombre")
base$Subtype <- as.factor(as.character(base$Subtype))
levels(base$Subtype) <- c("SN1", "SN2","SN3")


base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)

levels(base$PRIMARY_TUMOR) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

base$PRIMARY_TUMOR <- factor(base$PRIMARY_TUMOR,
levels = c("Estómago","Colorrectal", "Pulmón", "Páncreas", "I.Delgado"))

base$STAGE_IV <- as.factor(base$STAGE_IV)
levels(base$STAGE_IV) <- c("No", "Sí")
base$STAGE_IV <- factor(base$STAGE_IV, levels = c("Sí"))

df <- base %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        Subtipo = Subtype,
        Edad = AGE,
        Género = GENDER,
        Grado= GRADE,
        Metástasis = STAGE_IV,
        `Tumor Primario` = PRIMARY_TUMOR
      )




p1 <- forest_model(coxph(Surv(time, status) ~ ., df))


ggsave(p1,
file =paste0(outdir, "/Forest_todos.pdf"),
device = "pdf", width =12 , height = 10)




