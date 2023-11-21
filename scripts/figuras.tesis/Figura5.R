#this script need the environment figures

library(ggplot2)
library(tidyverse)
library(MOFA2)
library(ggpubr)



#load functions
source("functions/MOFA/MOFA.functions.R")
source("functions/MOFA/adictional.functions.R")
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


mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

factors_names(mod) <- gsub("Factor", "Factor.", factors_names(mod))


views_names(mod) <- c("Expresión", "Metilación")

groups_names(mod) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")


if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura5")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura5"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura5")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

features_cm <- c("pos.vs.neg.auc")

base <- base[,c(features_cm), drop = F]


colnames(base) <- c("Calidad.ARN")

df <- base[mod@samples_metadata$sample, , drop = F]

Z <- obtain_Z(mod)

df <- cbind(df, Z)
##########################################################

plt <- ggplot(df, aes_string(x = "Calidad.ARN", y = "Factor.3", group = 1)) +
    theme_classic() +
    geom_point(size = 3) +

    stat_cor(size = 10) +

    stat_smooth(color = "black") +
    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(margin = margin(10, 10, 10, 20),
        colour = "black"),
        axis.ticks.length.y = unit(10,"pt"),
        axis.ticks.length.x = unit(10,"pt")

    )



    ggsave(
        plot = plt, filename = paste0(
            outdir, "/Figura5.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 10,
        height = 10
    )

