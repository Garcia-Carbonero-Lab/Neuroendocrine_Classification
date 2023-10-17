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


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS12")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS12"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS12")

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


figS11.1 <- complete.boxplot(base,
factor = "Células.T",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Células.T_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.2 <- complete.boxplot(base,
factor = "T.CD8",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "T.CD8_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.3 <- complete.boxplot(base,
factor = "Linfocitos.citotóxicos",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Linfocitos.citotóxicos_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)



figS11.4 <- complete.boxplot(base,
factor = "Natural.killer",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Natural.killer_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.5 <- complete.boxplot(base,
factor = "Células.B",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Células.B_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.6 <- complete.boxplot(base,
factor = "Linaje.monocítico",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Linaje.monocítico_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.7 <- complete.boxplot(base,
factor = "Células.dendríticas",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Células.dendríticas_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.8 <- complete.boxplot(base,
factor = "Neutrófilos",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Neutrófilos_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)


figS11.9 <- complete.boxplot(base,
factor = "Células.endoteliares",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Células.endoteliares_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

figS11.10 <- complete.boxplot(base,
factor = "Fibroblastos",
column = "Subtipo",
colors = l_colors,
outdir = outdir,
flag = "Fibroblastos_complete",
points = T,
add.limits = F,
title = F,
remove.x = T,
legend = T)

pdf(paste0(outdir, "/FiguraS12.pdf"),
width = 20,
height = 45)

figS11.1 + figS11.2 + 
figS11.3 + figS11.4 +
figS11.5 + figS11.6 +
figS11.7 + figS11.8 +
figS11.9 + figS11.10 +
plot_layout(ncol = 2, guide = "collect")

dev.off()


