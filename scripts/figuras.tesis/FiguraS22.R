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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS22")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS22"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS22")

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



FiguraS22BA <- Kapplan_Meyer(cut_survival(base, 150,
time = "time",
event = "status"),
"GIPC2",
"time",
"status",
NULL,
outdir,
"FiguraS22BA",
c("red",
"blue"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS22BA$table <- FiguraS22BA$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS22BA$plot <- FiguraS22BA$plot + ggtitle("Mediana GIPC2") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/GIPC2",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS22BA)
dev.off()



df <- base[base$Subtipo == "SN1",]



FiguraS22BB <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"GIPC2",
"time",
"status",
NULL,
outdir,
"FiguraS22BB",
c("gold4", "gold"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS22BB$table <- FiguraS22BB$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS22BB$plot <- FiguraS22BB$plot + ggtitle("Mediana GIPC2 SN1") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/GIPC2_SN1",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS22BB)
dev.off()


df <- base[base$Subtipo == "SN2",]


FiguraS22BC <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"GIPC2",
"time",
"status",
NULL,
outdir,
"FiguraS22BC",
c("skyblue4", "skyblue"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS22BC$table <- FiguraS22BC$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS22BC$plot <- FiguraS22BC$plot + ggtitle("Mediana GIPC2 SN2") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/GIPC2_SN2",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS22BC)
dev.off()


df <- base[base$Subtipo == "SN3",]


FiguraS22BD <- Kapplan_Meyer(cut_survival(df, 150,
time = "time",
event = "status"),
"GIPC2",
"time",
"status",
NULL,
outdir,
"FiguraS22BD",
c("tomato4", "tomato"),
width= 8,
height=8,
risk.table = T,
font = 20,
fontsize = 7) 

FiguraS22BD$table <- FiguraS22BD$table  +
scale_y_discrete(labels = c("Bajo", "Alto")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
FiguraS22BD$plot <- FiguraS22BD$plot + ggtitle("Mediana GIPC2 SN3") + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/GIPC2_SN3",
"_Kapplan_Meyer.pdf"), width = 12, height = 12)
print(FiguraS22BD)
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



FiguraS22BE <- deg.boxplot(as.data.frame(t(expression)),
base,
"Tumor.Primario",
"Subtipo",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "GIPC2",
genes = "GIPC2"
)


FiguraS22BE2 <- deg.boxplot(as.data.frame(t(expression)),
base,
"Subtipo",
"Subtipo",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "GIPC2_complete",
genes = "GIPC2"
)

pdf(paste0(outdir, "/FiguraS22.1.pdf"),
width = 15,
height = 25)

((FiguraS22BA$plot / FiguraS22BA$table + plot_layout(heights = c(2,1))) | (FiguraS22BB$plot / FiguraS22BB$table + plot_layout(heights = c(2,1)))) /
((FiguraS22BC$plot / FiguraS22BC$table + plot_layout(heights = c(2,1))) | (FiguraS22BD$plot / FiguraS22BD$table + plot_layout(heights = c(2,1)))) /
 FiguraS22BE
dev.off()



pdf(paste0(outdir, "/FiguraS22.2.pdf"),
width = 15,
height = 25)

((FiguraS22BA$plot / FiguraS22BA$table + plot_layout(heights = c(2,1))) | (FiguraS22BB$plot / FiguraS22BB$table + plot_layout(heights = c(2,1)))) /
((FiguraS22BC$plot / FiguraS22BC$table + plot_layout(heights = c(2,1))) | (FiguraS22BD$plot / FiguraS22BD$table + plot_layout(heights = c(2,1)))) /
 FiguraS22BE2 

dev.off()
