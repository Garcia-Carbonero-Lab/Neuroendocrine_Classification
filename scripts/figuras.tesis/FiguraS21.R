#This script need the environment survival
library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)
library(glmnet)
library(rms)
library(fastDummies)

source("functions/cluster.analysis/survival.functions.R")
source("functions/main.functions.R")

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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS21")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS21"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS21")


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


rownames(expression) <- gsub("-", "_", rownames(expression))

#######################################################
#Cambiar nombres


base <- base %>%
    transmute(
        Edad = AGE,
        Género= GENDER,
        Grado = GRADE,
        Metástasis = STAGE_IV,
        Tumor.Primario = PRIMARY_TUMOR,
        time = OS.time,
        status = EXITUS
      )

base$Género <- as.factor(base$Género)
levels(base$Género) <- c("Mujer", "Hombre")

base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")

base$Metástasis <- as.factor(base$Metástasis)
levels(base$Metástasis) <- c("No", "Sí")

base$Tumor.Primario <- factor(base$Tumor.Primario,
levels = c("Estómago","Colorrectal", "Pulmón", "Páncreas", "I.Delgado"))


genes <- c("JCHAIN",
"PAK3",
"TCEA3",
"CAPN8",
"CXXC4",
"GIPC2",
"FAM177B",
"ANKRD30B")

genes_cov <- c(
    "RAP1GAP2",
    "FYB",
    "TCEAL2",
    "PAK3",
    "RNASE1",
    "HMGCS2",
    "HOXB_AS3",
    "CYP2C18",
    "CXXC4",
    "GIPC2",
    "FAM177B",
    "ANKRD30B",
    "HSPA4L",
    "GRB10"

)

df_genes <- merge(base, t(expression[genes,]), by = "row.names")

Hazard_Cox(base = df_genes,
class = genes[1],
OS = "time",
OS_event = "status",
covariates = genes[2:length(genes)],
outdir = outdir,
flag = "Final_sn_covariates")


df <- df_genes[,c("time","status",genes)]


p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_all_genes.pdf"),
device = "pdf", width =12 , height = 10)


Hazard_Cox(base = df_genes,
class = "Grado",
OS = "time",
OS_event = "status",
covariates = c("Metástasis", "Tumor.Primario", "Género", "Edad", genes),
outdir = outdir,
flag = "Final_sn_covariates")


df <- df_genes[,c("time","status",
"Edad",
"Género",
"Grado",
"Metástasis",
"Tumor.Primario",
genes)]

p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_all_genes_multivariante.pdf"),
device = "pdf", width =12 , height = 10)


df_genes_cov <- merge(base, t(expression[genes_cov,]), by = "row.names")

Hazard_Cox(base = df_genes_cov,
class = "Grado",
OS = "time",
OS_event = "status",
covariates = c("Metástasis", "Tumor.Primario", "Género", "Edad", genes_cov),
outdir = outdir,
flag = "Final_covariates")


df <- df_genes_cov[,c("time","status",
"Edad",
"Género",
"Grado",
"Metástasis",
"Tumor.Primario",
genes_cov)]


p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_genes_covariables.pdf"),
device = "pdf", width =8 , height = 8)



