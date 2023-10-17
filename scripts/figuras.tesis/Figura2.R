#environment figuras

library(MOFA2)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(survminer)
library(survival)
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura2")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura2"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura2")

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


colors_fcm <- list(
    "Género" = c("Mujer" = "lightpink", "Hombre" = "lightskyblue1"),
    "Grado" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "Metástasis" = c("No" = "forestgreen", "Sí" = "red4"),
    "X5HIAA" = c("Bajo" = "steelblue1", "Alto" = "brown2"),
    "CIMP" = c("Bajo" = "blue", "Alto" = "red")

)


################################################################################
#Figura 2A
columns <-  c("Género", "Grado", "Metástasis", "X5HIAA", "CIMP")

f2a <- heatmap.kruskal(mod,
columns,
outdir,
flag = "Kruskal.clinica")


##########################
#Figura 2B

f2b <- surv.model(mod,time = "Tiempo", event =  "Exitus", flag = "overall_survival")

#################################
#Figura 2C 


mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)


rownames(mcp)

rownames(mcp) <- c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Katural.killer",
"Células.B", 
"Linaje.monocítico",
"Células.dendríticas",
"Neutrófilos",
"Células.endoteliares",
"Fibroblastos")

mcp <- mcp[c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Natural.killer",
"Células.B", 
"Linaje.monocítico",
"Células.dendríticas",
"Neutrófilos",
"Células.endoteliares",
"Fibroblastos"),]


f2c <- immune_lm(as.data.frame(t(mcp)),
mod,
"mcp", outdir = outdir)



#######################################33
#Figura 2D
#######################################


expression_gsea <- list.files(paste0(wkdir,
"/MOFA/analysis_results/features/"), recursive = T)

expression_gsea <- expression_gsea[grep("GSEA_expression", expression_gsea)]
expression_gsea <- expression_gsea[grep("General", expression_gsea,invert = T)]
expression_gsea <- expression_gsea[grep("Subtypes", expression_gsea,invert = T)]
expression_gsea <- expression_gsea[grep("pathways", expression_gsea,invert = T)]




expression_gsea <- paste0(wkdir,
"/MOFA/analysis_results/features/", expression_gsea)

expression_gsea <- expression_gsea %>%
map(read_tsv) %>%
reduce(bind_rows)

f2d <- gsea.heatmp(expression_gsea, outdir)
##############################################################
#Create Panel

#############################################################
mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("StromalScore",
"ImmuneScore",
"ESTIMATEScore",
"TumorPurity",
"PRIMARY_TUMOR")

base <- base[,c(features_cm)]

variables <- c("Score Stromal",
"Score Inmune",
"Score ESTIMATE",
"Pureza Tumoral")

colnames(base) <- c(variables, "PRIMARY_TUMOR")

df <- base[mod@samples_metadata$sample, ]

Z <- obtain_Z(mod)

df <- cbind(df, Z)

mod@samples_metadata <- merge(mod@samples_metadata,
    df, by.x = 'sample', by.y = 'row.names')

continous.features(mod = mod,
p.th = 0.05,
columns = variables,
outdir = outdir,
flag = "estimate")





