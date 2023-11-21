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
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS13")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS13"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS13")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("PRIMARY_TUMOR",
"Subtype")

base <- base[,c(features_cm)]

variables <- c("Tumor.Primario",
"Subtipo")

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


mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)


rownames(mcp)

rownames(mcp) <- c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Natural.killer",
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

base <- merge(base, t(mcp), by = "row.names")


figs12.1 <- wrap.boxplot(base,
factor = "Células.T",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Células.T_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.2 <- wrap.boxplot(base,
factor = "T.CD8",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "T.CD8_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.3 <- wrap.boxplot(base,
factor = "Linfocitos.citotóxicos",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Linfocitos.citotóxicos_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.4 <- wrap.boxplot(base,
factor = "Natural.killer",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Natural.killer_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.5 <- wrap.boxplot(base,
factor = "Células.B",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Células.B_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.6 <- wrap.boxplot(base,
factor = "Linaje.monocítico",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Linaje.monocítico_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.7 <- wrap.boxplot(base,
factor = "Células.dendríticas",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Células.dendríticas_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.8 <- wrap.boxplot(base,
factor = "Neutrófilos",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Neutrófilos_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.9 <- wrap.boxplot(base,
factor = "Células.endoteliares",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Células.endoteliares_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

figs12.10 <- wrap.boxplot(base,
factor = "Fibroblastos",
column = "Subtipo",
group = "Tumor.Primario",
colors = l_colors,
outdir = outdir,
flag = "Fibroblastos_primario",
levels = c("Colorrectal", "Estómago",
"I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
points = T,
names = T)

pdf(paste0(outdir, "/FiguraS13.pdf"),
width = 30,
height = 30)

figs12.1 + figs12.2 + 
figs12.3 + figs12.4 +
figs12.5 + figs12.6 +
figs12.7 + figs12.8 +
figs12.9 + figs12.10 +
plot_layout(ncol = 2, guide = "collect")

dev.off()

