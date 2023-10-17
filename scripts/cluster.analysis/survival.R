#This script need the environment survival
library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)
source("functions/cluster.analysis/survival.functions.R")

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


if (dir.exists(paste0(wkdir, "/clustering/survival")) == F) {

    dir.create(paste0(wkdir, "/clustering/survival"))

}


outdir <- paste0(wkdir, "/clustering/survival")



base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


base$Subtype <- as.character(base$Subtype)
base$PRIMARY_TUMOR <- factor(base$PRIMARY_TUMOR,
levels = c("GASTRIC","COLORECTAL", "LUNG", "PANCREAS", "SMALL INTESTINE"))


# Create univariate cox model of features known to be related to OS
for (feature in c("Subtype",
"GRADE",
"STAGE_IV",
"AGE",
"PRIMARY_TUMOR",
"GENDER")) {

Hazard_Cox(base,
class = feature,
OS = "OS.time",
OS_event = "EXITUS",
covariates = NULL,
outdir = outdir,
flag = paste0("Univariate_",feature))
}

Hazard_Cox(base,
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("GRADE", "STAGE_IV", "AGE", "PRIMARY_TUMOR", "GENDER"),
outdir = outdir,
flag = "Multivariate_stage_IV")

#Create Forest plot
df <- base %>%
    transmute(
        time = OS.time,
        status = EXITUS,
        GENDER,
        GRADE,
        AGE,
        STAGE.IV = STAGE_IV,
        PRIMARY.TUMOR = PRIMARY_TUMOR,
        SUBTYPE = Subtype
    )

p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_all_stage_IV_.pdf"),
device = "pdf", width =12 , height = 10)


#Create Kapplan Meyer
Kapplan_Meyer(cut_survival(base, 150),
class = "Subtype",
OS = "OS.time",
OS_event = "EXITUS",
covariates = NULL,
outdir = outdir,
flag = "cut_150",
custcolor = c("gold", "cornflowerblue", "tomato1")
)

#Kapplan Meyer 150 for each tissue

for (tissue in unique(base$PRIMARY_TUMOR)){

    df <- base[base$PRIMARY_TUMOR == tissue,]

    Kapplan_Meyer(cut_survival(df, 150),
    class = "Subtype",
    OS = "OS.time",
    OS_event = "EXITUS",
    covariates = NULL,
    outdir = outdir,
    flag = paste0("cut_150_", tissue),
    custcolor = c("gold", "cornflowerblue", "tomato1")
    )


}

#Kapplan Meyer 150 for stage IV and stage I-III

for (meta in unique(base$STAGE_IV)){

    df <- base[base$STAGE_IV == meta,]

    Kapplan_Meyer(cut_survival(df, 150),
    class = "Subtype",
    OS = "OS.time",
    OS_event = "EXITUS",
    covariates = NULL,
    outdir = outdir,
    flag = paste0("cut_150_Metastasis", meta),
    custcolor = c("gold", "cornflowerblue", "tomato1")
    )


}



data <- read.table(paste0(outdir, "/Multivariate_Hazard_ratio_factors.txt"),
sep = "\t",
header = T)
