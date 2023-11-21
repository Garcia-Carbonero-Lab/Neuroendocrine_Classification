#This script needs the environment figuras

library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(patchwork)
#load functions



source("functions/main.functions.R")
source("functions/cluster.analysis/corr.met.exp.functions.R")
source("functions/cluster.analysis/additional.functions.R")
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS20")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS20"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS20")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


base$Subtype <- as.factor(base$Subtype)
levels(base$Subtype) <- c("SN1", "SN2", "SN3")

base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)
levels(base$PRIMARY_TUMOR) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")

base$PRIMARY_TUMOR  <- factor(base$PRIMARY_TUMOR, levels =c("Colorrectal",
"Estómago", "I.Delgado", "Páncreas", "Pulmón") )


base <- base %>%
    transmute(
        Subtipo = Subtype,
        Tumor.Primario = PRIMARY_TUMOR,
      )

exp <- read.table(paste0(wkdir,
"/clustering/deg/Expression.deg.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)



met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)

met_anot <- create.met.anot(epic)


genes <- c("GIPC2",
"CHGA",
"PIGR",
"NOVA1",
"LGALS4",
"PTPRN",
"CHGB")

probes <- c("cg26947773",
"cg19433225",
"cg02244697",
"cg04974854",
"cg19419519",
"cg03970036",
"cg20851464")


crp <- correlation.plot(exp,
met,
met_anot,
base,
genes,
probes,
outdir,
"Tumor.Primario")


pdf(paste0(outdir, "/FiguraS20.pdf"),
width = 32,
height = 20)

crp[[1]] + crp[[2]] + 
crp[[3]] + crp[[4]] +
crp[[5]] + crp[[6]]  + plot_layout(ncol = 2, guide = "collect")

dev.off()

#############################################################3
#Table

primario <- read.table(paste0(wkdir, "/Figuras_tesis/Figura13/Correlation_exp_met_primary_filtrado.txt"),
sep = "\t",
header = T)

primario <- primario[primario$CpG %in% probes,]

write.table(primario, paste0(outdir, "/Tabla_primarios.tsv"),
sep = "\t", row.names = F, col.names = T)


subtipo <- read.table(paste0(wkdir, "/Figuras_tesis/Figura13/Correlation_exp_met_subtypes_filtrado.txt"),
sep = "\t",
header = T)

subtipo <- subtipo[subtipo$CpG %in% probes,]

write.table(subtipo, paste0(outdir, "/Tabla_subtipo.tsv"),
sep = "\t", row.names = F, col.names = T)



