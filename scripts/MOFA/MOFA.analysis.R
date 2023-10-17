# IMPORTANT NOTES
#Script to analyse MOFA model
#For this script we need MOFA.yml enviroment

#LOAD packages
library(MOFA2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(fgsea)
library(reshape2)
library(rstatix)

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

#We create output folder
if (dir.exists(paste0(wkdir, "/MOFA/analysis_results")) == F) {

    dir.create(paste0(wkdir, "/MOFA/analysis_results"))

}



#Prepare anotation files

#We open epic hg38 manifest 
#download it here:
#https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/archive/202209/EPIC.hg38.manifest.tsv.gz

epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)
    

met_anot <- create.met.anot(epic)


#Open clinical data and model        

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


#CHANGE COLUMNS CLASS
base$X5.hydroxyindoleacetic_acid..mg.24H. <- as.numeric(base$X5.hydroxyindoleacetic_acid..mg.24H.)


# We indicate the out folder
outdir <- paste0(wkdir, "/MOFA/analysis_results/")
# Introduce clinical data and Factor value in samples_metadata
mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))

df <- base[mod@samples_metadata$sample, ]

Z <- obtain_Z(mod)

df <- cbind(df, Z)

mod@samples_metadata <- merge(mod@samples_metadata,
    df, by.x = 'sample', by.y = 'row.names')


#anotate metilation
mod <- anotate.metilation(mod, met_anot)


#############################################
#General Results

# Create folder of general results
if (dir.exists(paste0(outdir, "/general_results")) == F) {

    dir.create(paste0(outdir, "/general_results"))

}


#write samples_metadata
write.table(df,
paste0(outdir, "/general_results/Clinica_table_factors.txt"),
sep = "\t",
row.names = T,
col.names = NA
)

# Plot variance

pve <- plot_variance_explained(mod, x = 'group', y = 'factor')
pve <- theme_general(pve)

ggsave(filename = paste0(outdir, "/general_results/Variance_groups.pdf"),
device = "pdf"
           ,dpi =  300, width = 29, height = 15)


calculate.mean.variance(mod, paste0(outdir, "/general_results"))
################################################


################################################
#Associated discrete variables

if (dir.exists(paste0(outdir, "/discrete.variables")) == F) {

    dir.create(paste0(outdir, "/discrete.variables"))

}


features_cm <- c("GENDER",
"GRADE",
"STAGE_IV",
"FUNCTIONAL_dichotomic",
"X5HIAA_HIGH",
"CIMP.status")

colors_fcm <- list(
    "GENDER" = c("FEMALE" = "lightpink", "MALE" = "lightskyblue1"),
    "GRADE" = c("G1" = "green", "G2" = "yellow", "G3" = "red"),
    "STAGE_IV" = c("YES" = "red4", "NO" = "forestgreen"),
    "FUNCTIONAL_dichotomic" = c("YES" = "steelblue1", "NO" = "brown2"),
    "X5HIAA_HIGH" = c("YES" = "steelblue1", "NO" = "brown2"),
    "CIMP.status" = c("Low" = "blue", "Inter" = "green", "High" = "red")

)


##########################################################
#Associated discrete features

discrete.features(mod = mod,
p.th = 0.05,
columns = features_cm,
outdir = outdir,
flag = "clinical",
colors = colors_fcm,
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))


##################################################################
#Associate continous variables


if (dir.exists(paste0(outdir, "/continous.variables/")) == F) {

    dir.create(paste0(outdir, "/continous.variables/"))

}

continous_features <- c(
"Ki.67...",
"X5.hydroxyindoleacetic_acid..mg.24H.",
"AGE", "StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"
)

continous.features(mod = mod,
p.th = 0.05,
columns = continous_features,
outdir = outdir,
flag = "normal_features",
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))


###############################################
#plots correlation genes Alcala Factor 4


genes <- c("ASCL1", "DLL3", "ROBO1", "SLIT1", "ANGPTL3", "ERBB4",
"UGT2A3", "UGT2B15", "UGT2B7", "UGT2B11", "UGT2B4", "UGT2B17", "OTP", "HNF1A")

if (dir.exists(paste0(outdir, "/continous.variables/genes_alcala")) == F) {

    dir.create(paste0(outdir, "/continous.variables/genes_alcala"))

}


# open expression data
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

# select genes
expgenes <- expression[rownames(expression) %in% genes,]


# add genes information to model metadata
mod@samples_metadata <- merge(mod@samples_metadata, t(expgenes),
by.x = "sample", by.y = "row.names")

samples <- mod@samples_metadata$sample[mod@samples_metadata$PRIMARY_TUMOR == "LUNG"]

# Calculate correlatuon between genes and factors
continous.features(mod = subset_samples(mod, samples),
p.th = 0.05,
columns = rownames(expgenes),
outdir = outdir,
flag = "genes_alcala",
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))



for (gene in rownames(expgenes)){


if (dir.exists(paste0(outdir, "/continous.variables/genes_alcala/",gene)) == F) {

    dir.create(paste0(outdir, "/continous.variables/genes_alcala/", gene))

}

# To obtain the correlation plot
correlation.plot(mod, gene,
factor = "Factor4",
"genes_alcala",
outdir =  outdir,
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))


plt <- plot_factor(mod,
                factors = "Factor4",
                color_by = gene
            ) +
                theme(
                    text = element_text(size = 18),
                    axis.text.x = element_text(angle = 90,
                    vjust = 0.5, hjust = 1)
                )
ggsave(
    plot = plt, filename = paste0(
    outdir, "/continous.variables/genes_alcala/", gene,
            "/Factors_Scatterplot_genes_alcala_", gene, ".pdf"
    ),
    device = "pdf",
    dpi = 500, width = 10,
     height = 5
            )
}


col_alcala <- list("PRIMARY_TUMOR" = c("SMALL INTESTINE" = "yellow",
"LUNG" = "blue",
"COLORECTAL" = "orange",
"GASTRIC" = "pink",
"PANCREAS" = "green"),
"BIOPSY_METASTASIS_LOCATION" = c("PERITONEUM" = "lightsalmon", 
    "LYMPH NODE" = "lightskyblue1",
    "LIVER" = "brown",
    "LUNG" = "blue",
    "OVARY" = "lightpink",
    "SKIN" = "gold",
    "PAROTID GLAND" = "magenta")
)



anot_legend_param <-  list(
                "PRIMARY_TUMOR" = list(
                
                title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                grid_height = unit(2, "cm"),
                grid_width = unit(2, "cm"),
                gap = unit(40, "mm"),
                ncol = 2

            ),

            "BIOPSY_METASTASIS_LOCATION" = list(
                
                title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                grid_height = unit(2, "cm"),
                grid_width = unit(2, "cm"),
                gap = unit(5, "cm"),
                ncol = 2

            ),

                factor = list(title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40),
                legend_height = unit(8, "cm"),
                grid_width = unit(2, "cm"),
                title_position= "leftcenter-rot")
        )






modlung <- subset_samples(mod, rownames(base)[base$PRIMARY_TUMOR == "LUNG"])

heatmap_continous(mod = modlung,
columns = rownames(expgenes),
flag = "genes_alcala",
outdir = outdir,
colorv = col_alcala,
size = list("width" = 24, "height" = 20),
anot_legend_param = anot_legend_param)

###################################
#Analisis immune population

# load mcp counter immune infiltration data
mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)



mod@samples_metadata <- merge(mod@samples_metadata, t(mcp),
by.x = "sample", by.y = "row.names")



continous.features(mod = mod,
p.th = 0.05,
columns = rownames(mcp),
outdir = outdir,
flag = "mcp",
group = "PRIMARY_TUMOR",
levels = c("COLORECTAL", "GASTRIC", "LUNG", "PANCREAS", "SMALL INTESTINE"))


col_immune <- list("PRIMARY_TUMOR" = c("SMALL INTESTINE" = "yellow",
"LUNG" = "blue",
"COLORECTAL" = "orange",
"GASTRIC" = "pink",
"PANCREAS" = "green"),
"BIOPSY_METASTASIS_LOCATION" = c("PERITONEUM" = "lightsalmon", 
    "LYMPH NODE" = "lightskyblue1",
    "LIVER" = "brown",
    "LUNG" = "blue",
    "OVARY" = "lightpink",
    "SKIN" = "gold",
    "PAROTID GLAND" = "magenta")
)


heatmap_continous(mod = mod,
columns =  c("StromalScore", "ImmuneScore", "ESTIMATEScore","TumorPurity"
),
flag = "normal_features",
outdir = outdir,
colorv = col_immune,
size = list("width" = 2400, "height" = 1100),
anot_legend_param)


heatmap_continous(mod = mod,
columns = rownames(mcp),
flag = "mcp",
outdir = outdir,
colorv = col_immune,
size = list("width" = 2400, "height" = 1900),
anot_legend_param)


immune_lm(as.data.frame(t(mcp)),
mod,
"mcp", outdir)

########################################################3
#Survival
if (dir.exists(paste0(outdir, "/survival/")) == F) {

    dir.create(paste0(outdir, "/survival/"))

}



surv.model(mod,time = "OS.time", "EXITUS", flag = "overall_survival", outdir)
surv.model(mod,time = "PFS.time", "RELAPSE_DICHOTOMIC", flag = "pfs_survival", outdir)





if (dir.exists(paste0(outdir, "/features/")) == F) {

    dir.create(paste0(outdir, "/features/"))

}


#############################################################
#Methylation importance estatus

pathways <- gmtPathways(paste0(wkdir, "/genesets/General_genesets.gmt"))

#Obtain the important features in each omics
important.features(mod,
outdir = outdir,
gmt = pathways,
met_anot = met_anot)


importance_met <- read.table(
paste0(outdir, "/features/Methylation_Importance_features.txt"),
sep = "\t",
header = T
)


importance_met <- importance_met[!duplicated(importance_met$id),]
rownames(importance_met) <- importance_met$id




if (dir.exists(paste0(outdir, "/features/methylation_statis")) == F) {

    dir.create(paste0(outdir, "/features/methylation_statis"))

}


#########################################################################
#GSEA

#open gmts
gmts <- list.files(paste0(wkdir,"/genesets"))
gmts <- gmts[grep(".gmt", gmts)]
gmts <- paste0(wkdir, "/genesets/", gmts)


#Prepare methylation data


methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F)



gsea_factors(expression, mod, gmts = gmts, outdir,
type = "expression", met_anot = NULL, flag = "expression" )


gsea_factors(methylation, mod, gmts = gmts, outdir,
type = "methylation", met_anot = met_anot, flag = "methylation")

###################################################################
