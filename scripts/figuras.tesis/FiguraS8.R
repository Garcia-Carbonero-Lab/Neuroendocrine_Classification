
library(MOFA2)
library(survival)
library(survminer)
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


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS8")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS8"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS8")

colnames(base)


base <- base[,c("Subtype",
"PRIMARY_TUMOR",
"OS.time",
"EXITUS",
"STAGE_TNM")]

colnames(base) <- c("Subtipo", 
"Tumor.Primario",
"OS.time",
"EXITUS",
"Estadio")



base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón", 
"Páncreas",
"I.Delgado")
base$Estadio <- as.factor(base$Estadio)

base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")


plots <- list()


for(tissue in levels(base$Tumor.Primario)){
print(tissue)
df <- base[base$Tumor.Primario == tissue,]

figS8 <- Kapplan_Meyer(cut_survival(df, 150),
"Subtipo",
"OS.time",
"EXITUS",
NULL,
outdir,
"FiguraS8",
c("gold",
"cornflowerblue",
"tomato1"),
width= 8,
height=8,
risk.table = T,
font = 20) 

figS8$table <- figS8$table  +
scale_y_discrete(labels = c("SN3", "SN2", "SN1")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
figS8$plot <- figS8$plot + ggtitle(tissue) + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))


plots[[tissue]] <- figS8


pdf(paste0(outdir, "/",tissue,
"_Kapplan_Meyer.pdf"), width = 8, height = 8)
print(figS8)
dev.off()

}

pdf(paste0(outdir, "/FiguraS8.pdf"), width = 15, height = 27)
arrange_ggsurvplots(plots, print = T, ncol = 2, nrow = 3)
dev.off()
##############################################
#stage
for(stage in levels(base$Estadio)){
print(stage)

df <- base[base$Estadio == stage,]
df <- df[complete.cases(df$Estadio),]


figS8.2 <- Kapplan_Meyer(cut_survival(df, 150),
"Subtipo",
"OS.time",
"EXITUS",
NULL,
outdir,
"FiguraS8.2",
c("gold",
"cornflowerblue",
"tomato1"),
width= 8,
height=8,
risk.table = T,
font = 20) 

figS8.2$table <- figS8.2$table  +
scale_y_discrete(labels = c("SN3", "SN2", "SN1")) +
theme(axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.y = ggtext::element_markdown(size = 20),
    plot.margin = unit(c(1,1,1,1), "cm")
        ) 
figS8.2$plot <- figS8.2$plot + ggtitle(stage) + theme(
    plot.title = element_text(hjust = 0.5,
    vjust= 1, size = 30),
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(margin = margin(20, 10, 20, 20)),
    plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir, "/",stage,
"_Kapplan_Meyer.pdf"), width = 8, height = 8)
print(figS8.2)
dev.off()

}
