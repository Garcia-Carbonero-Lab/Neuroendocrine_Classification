
library(MOFA2)
library(ggplot2)
library(ggpubr)
library(circlize)
library(tidyverse)
library(rstatix)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura6")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura6"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura6")

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
base[,"X5HIAA"] <- factor(base[,"X5HIAA"], levels = c("Alto", "Bajo"))


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

###################################################################



###########36A
source("functions/figuras/figuras.functions.R")

colors_fcm <- list(
    "Género" = c("Mujer" = "lightpink", "Hombre" = "lightskyblue1"),
    "Grado" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "Metástasis" = c("No" = "forestgreen", "Sí" = "red4"),
    "X5HIAA" = c("Alto" = "brown2", "Bajo" = "steelblue1"),
    "CIMP" = c("Alto" = "red","Bajo" = "blue")

)



fig6a1 <- complete.boxplot(df,
factor = "Factor.4",
column = "Género",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_Género_complete.pdf")


fig6a2 <- complete.boxplot(df,
factor = "Factor.4",
column = "X5HIAA",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_X5HIAA_complete.pdf",
remove.y = T)

fig6b1 <- wrap.boxplot(df,
factor = "Factor.4",
column = "Género",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_Género_wrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)


fig6b2 <- wrap.boxplot(df,
factor = "Factor.4",
column = "X5HIAA",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_X5HIAA_wrapper",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)



fig6b12 <- scatterplot(df,
factor = "Factor.4",
column = "Género",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_Género_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F)


fig6b22 <- scatterplot(df,
factor = "Factor.4",
column = "X5HIAA",
group = "Tumor.Primario",
colors = colors_fcm,
outdir = outdir,
flag = "Factor4_X5HIAA_scatter",
levels = c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"),
wrap = F,
remove.y = T)

pdf(file = paste0(outdir, "/Figura6A.pdf"), height = 30, width = 35)
fig6a1 <- fig6a1 + labs(tag = "A")
fig6b1 <- fig6b1 + labs(tag = "B")

fig6a2 <- fig6a2 + labs(tag = " ")
fig6b2 <- fig6b2 + labs(tag = " ")



patch1 <- fig6a1 / (fig6b1 / fig6b12 + plot_layout(guides = "auto")) +
plot_layout(guides = "auto", heights = c(1, 2))

patch2 <- fig6a2 / (fig6b2 / fig6b22 + plot_layout(guides = "auto")) +
plot_layout(guides = "auto", heights = c(1, 2))



(patch1 | patch2)


dev.off()
#############################6B

#Fig6B



genes <- c("ASCL1", "DLL3", "ROBO1", "SLIT1", "ANGPTL3", "ERBB4",
"UGT2A3", "UGT2B15", "UGT2B7", "UGT2B11", "UGT2B4", "UGT2B17", "OTP", "HNF1A")

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

expgenes <- expression[rownames(expression) %in% genes,]


mod@samples_metadata <- merge(mod@samples_metadata, t(expgenes),
by.x = "sample", by.y = "row.names")

samples <- mod@samples_metadata$sample[mod@samples_metadata$PRIMARY_TUMOR == "LUNG"]


modlung <- subset_samples(mod, rownames(base)[base$Tumor.Primario == "Pulmón"])


Z <- Reduce(rbind, modlung@expectations$Z)
Z <- Z[modlung@samples_metadata$sample, ]

dbase <- modlung@samples_metadata
Z <- Z[order(Z[, "Factor.4"], decreasing = T), ]

df <- dbase[match(rownames(Z), dbase$sample), ]
data <- df[, genes]
rownames(data) <- df$sample
data <- t(zscore.rows2(t(data)))

df$Tumor.Primario <- droplevels(df$Tumor.Primario)

heatmap_continous(mod = modlung,
columns = rownames(expgenes),
flag = "genes_alcala",
outdir = outdir,
colorv = col_immune,
size = list("width" = 24, "height" = 20))

col <- list("Tumor.Primario" = c(
"Pulmón" = "blue"),
"Factor.4" =  colorRamp2(
            c(
                min(df[, "Factor.4"], na.rm = T),
                0,
                max(df[, "Factor.4"], na.rm = T)
            ),
            c("blue", "white", "red"))
)



anot <- df[, names(col)]
rownames(anot) <- rownames(df)


ha_row <- rowAnnotation(
            df = df[,c("Tumor.Primario","Factor.4")],
            show_legend = F,
            show_annotation_name = c(Tumor.Primario = F,
            Factor.4 = F),
            col = col,
            annotation_name_gp = gpar(fontsize = 50),
            border = unit(50, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 50, col = "black"),
            annotation_name_side = "top")

lgd1 = Legend(labels = levels(df[,"Tumor.Primario"]),
legend_gp = gpar(fill = c("blue")),
title = "Tumor Primario",
title_gp = gpar(fontsize = 50),
labels_gp = gpar(fontsize = 50),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1
)

lgd2 = Legend(
col_fun = col[["Factor.4"]],
title = "Factor 4",
title_gp = gpar(fontsize = 50),
labels_gp = gpar(fontsize = 50),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1
)

hm <- Heatmap(as.matrix(data),
            left_annotation =  ha_row,
            show_row_names = F,
            show_column_names = T,
            col = colorRamp2(
                c(max(data), 0, min(data)),
                c("red", "white", "blue")
            ),
            name = "Heatmap",
            gap = unit(5, "mm"),
            cluster_rows = F,
            show_column_dend = F,
            cluster_columns = T,
            clustering_distance_columns = "spearman",
            heatmap_legend_param = list(
                title = "Z score(expresión)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 50),
                labels_gp = gpar(fontsize = 50,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm")
            ),
            column_names_side = "bottom",
            column_names_gp = gpar(fontsize = 50),
            column_names_rot = 45
        )


pd = packLegend(list = list(lgd1, lgd2) , 
row_gap = unit(3, "cm"))



        p <- draw(hm, heatmap_legend_side = "bottom",
        padding = unit(c(2, 120, 2, 2), "mm"),
        annotation_legend_side = "right",
        annotation_legend_list = pd)

        pdf(paste0(
            outdir, "/Figura6B.pdf"
        ),
        width = 45, 15
        )
        print(p)
        dev.off()


