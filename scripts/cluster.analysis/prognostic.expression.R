#This script need the environment survival
library(survival)
library(survminer)
library(forestplot)
library(tidyverse)
library(aod)
library(forestmodel)
library(glmnet)
library(rms)
library(fastDummies)

source("functions/cluster.analysis/survival.functions.R")
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


if (dir.exists(paste0(wkdir, "/clustering/prognostic.features")) == F) {

    dir.create(paste0(wkdir, "/clustering/prognostic.features"))

}


outdir <- paste0(wkdir, "/clustering/prognostic.features")



base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

#Select genes


deg1 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group1-Group2.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg2 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group1-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg3 <- read.table(paste0(wkdir,
"/clustering/deg/TopTable_deg_cov_primary_Group2-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)


deg3_filter <- deg3[deg3$adj.P.Val < 0.1 & abs(deg3$logFC) > 0.2,]
deg2_filter <- deg2[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 0.5,]
deg1_filter <- deg1[deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 0.5,]

genes <- unique(c(rownames(deg1_filter), rownames(deg2_filter), rownames(deg3_filter)))
rownames(expression) <- gsub("-", "_", rownames(expression))
genes <- gsub("-","_",genes)

expression <- expression[genes,]
################################################################################
#1
#Uses only LASSO as variable selector


genes_pen <- 
cox.lasso(base = base,
data = expression,
OS = "OS.time",
OS_event = "EXITUS",
covariates = NULL,
outdir = outdir,
flag = "Cox_Lasso_sn_covaraites",
alpha = 1)


genes_pen <- data.frame("names" = names(genes_pen),
"coeficcient" = genes_pen)

write.table(genes_pen,
file = paste0(outdir, "/Lasso_cox_Coeficcients_lambda.min.txt"),
sep = "\t",
row.names = F,
col.names = T)

base <- base[, !colnames(base) %in% "STAGE_IV_UKNOWN"]

base_dummy <- base[complete.cases(base$STAGE_IV),]

base_dummy <- dummy_cols(base_dummy,
select_columns = c("GRADE", "STAGE_IV", "PRIMARY_TUMOR", "GENDER"),
           remove_selected_columns = TRUE)

rownames(base_dummy) <- rownames(base[complete.cases(base$STAGE_IV),])


covariates <- c("GRADE", "STAGE_IV", "PRIMARY_TUMOR", "GENDER")

covariates <- unlist(lapply(covariates,function(x){
    return(colnames(base_dummy)[grep(x,colnames(base_dummy))])
}))


gene_selection_cov_forced <- 
cox.lasso(base = base_dummy,
data = expression,
OS = "OS.time",
OS_event = "EXITUS",
covariates =  c(covariates, "AGE"),
outdir = outdir,
flag = "Cox_Lasso_covariates_forced",
alpha = 1,
forced = T)



genes_pen_cov_forced <- data.frame("names" = names(gene_selection_cov_forced),
"coeficcient" = gene_selection_cov_forced)

write.table(genes_pen_cov_forced,
file = paste0(outdir, "/Lasso_cox_Coeficcients_lambda.min_cov_forced.txt"),
sep = "\t",
row.names = F,
col.names = T)


gene_selection_cov <- 
cox.lasso(base = base_dummy,
data = expression,
OS = "OS.time",
OS_event = "EXITUS",
covariates =  c(covariates, "AGE"),
outdir = outdir,
flag = "Cox_Lasso_covariates",
alpha = 1,
forced = F)


gene_selection_cov <- data.frame("names" = names(gene_selection_cov),
"coeficcient" = gene_selection_cov)


write.table(gene_selection_cov,
file = paste0(outdir, "/Lasso_cox_Coeficcients_lambda.min_cov.txt"),
sep = "\t",
row.names = F,
col.names = T)
################################################################################

################################################################################
#Final models

genes <- c("JCHAIN",
"PAK3",
"TCEA3",
"CAPN8",
"CXXC4",
"GIPC2",
"FAM177B",
"ANKRD30B")

genes_cov <- c(
    "RAP1GAP2",
    "FYB",
    "TCEAL2",
    "PAK3",
    "RNASE1",
    "HMGCS2",
    "HOXB_AS3",
    "CYP2C18",
    "CXXC4",
    "GIPC2",
    "FAM177B",
    "ANKRD30B",
    "HSPA4L",
    "GRB10"

)


df_genes <- merge(base, t(expression[genes,]), by = "row.names")

Hazard_Cox(base = df_genes,
class = genes[1],
OS = "OS.time",
OS_event = "EXITUS",
covariates = genes[2:length(genes)],
outdir = outdir,
flag = "Final_sn_covariates")



df <- df_genes[,c("OS.time","EXITUS",genes)]
colnames(df) <- c("time","status",genes)


p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_all_genes.pdf"),
device = "pdf", width =12 , height = 10)


df_genes_cov <- merge(base, t(expression[genes_cov,]), by = "row.names")

Hazard_Cox(base = df_genes_cov,
class = "GRADE",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("STAGE_IV", "PRIMARY_TUMOR", "GENDER", "AGE", genes_cov),
outdir = outdir,
flag = "Final_covariates")



df <- df_genes_cov[,c("OS.time","EXITUS",
"STAGE_IV", "PRIMARY_TUMOR", "GENDER", "AGE", "GRADE", genes_cov)]
colnames(df) <- c("time","status","STAGE.TNM", "PRIMARY.TUMOR", "GENDER", "AGE", "GRADE", genes_cov)


p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_genes_covariables.pdf"),
device = "pdf", width =12 , height = 10)
############################################
genes_multivariate <- c("GSTA1",
"GIPC2")

df_genes_multi<- merge(base, t(expression[genes_multivariate,]), by = "row.names")

Hazard_Cox(base = df_genes_multi,
class = "GRADE",
OS = "OS.time",
OS_event = "EXITUS",
covariates = c("STAGE_TNM", "PRIMARY_TUMOR", "GENDER", "AGE", genes_multivariate),
outdir = outdir,
flag = "Final_covariates_genes_multivariate")




df <- df_genes_multi[,c("OS.time","EXITUS",
"STAGE_TNM", "PRIMARY_TUMOR", "GENDER", "AGE", "GRADE", genes_multivariate)]
colnames(df) <- c("time","status","STAGE.TNM", "PRIMARY.TUMOR", "GENDER", "AGE", "GRADE", genes_multivariate)


p <- forest_model(coxph(Surv(time, status) ~ ., df))

ggsave(p,
file =paste0(outdir, "/Forest_plot_genes_covariables_multivariate.pdf"),
device = "pdf", width =12 , height = 10)
################################################################################
