#this script needs the environment deg

library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/cluster.analysis/clinical.functions.R")

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

if (dir.exists(paste0(wkdir, "/clustering/immune")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune"))

}

outdir <- paste0(wkdir, "/clustering/immune")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)


base$Subtype <- as.factor(as.character(base$Subtype))


base <- merge(base, t(mcp), by = "row.names")


# Study differences of immune cells using MCP counter data
kruskal <- kruskal(base,
"Subtype",
rownames(mcp),
outdir = outdir,
flag = "mcp")


if (dir.exists(paste0(wkdir, "/clustering/immune/complete.plots")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune/complete.plots"))

}
rownames(base) <- base$Row.names


# Differences of immune cells in each subtype
wilcox_plots(data = as.data.frame(t(mcp)),
base = base,
group.x = "Subtype",
group.fill =  "Subtype",
color.fill = c("gold", "cornflowerblue", "tomato1"),
outdir= paste0(outdir, "/complete.plots"),
flag = "MCP_Complete")

if (dir.exists(paste0(wkdir, "/clustering/immune/primary")) == F) {

    dir.create(paste0(wkdir, "/clustering/immune/primary"))

}
# Differences of immune cells in each subtype by primary tumors

wilcox_plots(data = as.data.frame(t(mcp)),
base = base,
group.x = "PRIMARY_TUMOR",
group.fill =  "Subtype",
color.fill = c("gold", "cornflowerblue","tomato1"),
outdir= paste0(outdir, "/primary"),
flag = "MCP_primary")

l_colors <- list(
    "Subtype" = c("1" = "gold",
    "2" = "cornflowerblue",
    "3" = "tomato1"),
    "PRIMARY_TUMOR" = c(
        "COLORECTAL" = "orange",
        "LUNG" = "blue",
        "SMALL INTESTINE" = "yellow",
        "PANCREAS" = "green",
        "GASTRIC" = "pink"
    ))



anot_param_legend <- list(
    "Subtype" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    
    ),
    "PRIMARY_TUMOR" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    )
)


# Heatmap of immune cells
heatmap(mcp,
base = base,
outdir = outdir,
col = l_colors,
size = list("width" = 50, "height" = 45),
anot_param_legend = anot_param_legend,
flag = "mcp",
anot_h = c(40,20),
column = "Subtype",
annotation_name_gp = 90,
names.size = 20
)
