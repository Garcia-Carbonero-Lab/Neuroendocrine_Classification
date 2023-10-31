#This script need the environment figuras
library(survival)
library(survminer)
library(tidyverse)
library(patchwork)
library(aod)
library(limma)
library(ggplot2)
library(ggpubr)
library(rstatix)

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS23")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS23"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS23")

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

rownames(expression) <- gsub("-","_",rownames(expression))


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

genes <- c(
    "GIPC2",
    "PAK3"
)

medians <- generate_medians(expression, genes)



base <- base %>%
    transmute(
        Subtipo = Subtype,
        Tumor.Primario = PRIMARY_TUMOR,
        time = OS.time,
        status = EXITUS
      )


base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")

base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1",
"SN2", "SN3")



base <- merge(base,medians, by = "row.names")



FiguraS23A <- Kapplan_Meyer(cut_survival(base, 150,
time = "time",
event = "status"),
"PAK3",
"time",
"status",
NULL,
outdir,
"FiguraS23A",
c("red",
"blue"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS23A$table <- FiguraS23A$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS23A$plot <- FiguraS23A$plot + ggtitle("Mediana PAK3") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/PAK3",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS23A)
dev.off()



df <- base[base$Subtipo == "SN1",]



FiguraS23B <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"PAK3",
"time",
"status",
NULL,
outdir,
"FiguraS23B",
c("gold4", "gold"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS23B$table <- FiguraS23B$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS23B$plot <- FiguraS23B$plot + ggtitle("Mediana PAK3 SN1") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/PAK3_SN1",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS23B)
dev.off()


df <- base[base$Subtipo == "SN2",]


FiguraS23C <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"PAK3",
"time",
"status",
NULL,
outdir,
"FiguraS23C",
c("skyblue4", "skyblue"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS23C$table <- FiguraS23C$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS23C$plot <- FiguraS23C$plot + ggtitle("Mediana PAK3 SN2") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/PAK3_SN2",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS23C)
dev.off()


df <- base[base$Subtipo == "SN3",]


FiguraS23D <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"PAK3",
"time",
"status",
NULL,
outdir,
"FiguraS23D",
c("tomato4", "tomato"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS23D$table <- FiguraS23D$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS23D$plot <- FiguraS23D$plot + ggtitle("Mediana PAK3 SN3") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/PAK3_SN3",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS23D)
dev.off()


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)



base <- base %>%
    transmute(
        Subtipo = Subtype,
        Tumor.Primario = PRIMARY_TUMOR,
        time = OS.time,
        status = EXITUS
      )


base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")

base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1",
"SN2", "SN3")



FiguraS23E <- deg.boxplot(as.data.frame(t(expression)),
base,
"Tumor.Primario",
"Subtipo",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "PAK3",
genes = "PAK3"
)


FiguraS23E2 <- deg.boxplot(as.data.frame(t(expression)),
base,
"Subtipo",
"Subtipo",
cov = "Tumor.Primario",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "PAK3_complete",
genes = "PAK3"
)

pdf(paste0(outdir, "/FiguraS23.1.pdf"),
width = 15,
height = 25)

((FiguraS23A$plot / FiguraS23A$table + plot_layout(heights = c(2,1))) | (FiguraS23B$plot / FiguraS23B$table + plot_layout(heights = c(2,1)))) /
((FiguraS23C$plot / FiguraS23C$table + plot_layout(heights = c(2,1))) | (FiguraS23D$plot / FiguraS23D$table + plot_layout(heights = c(2,1)))) /
 FiguraS23E 

dev.off()



pdf(paste0(outdir, "/FiguraS23.2.pdf"),
width = 15,
height = 25)

((FiguraS23A$plot / FiguraS23A$table + plot_layout(heights = c(2,1))) | (FiguraS23B$plot / FiguraS23B$table + plot_layout(heights = c(2,1)))) /
((FiguraS23C$plot / FiguraS23C$table + plot_layout(heights = c(2,1))) | (FiguraS23D$plot / FiguraS23D$table + plot_layout(heights = c(2,1)))) /
 FiguraS23E2 

dev.off()

