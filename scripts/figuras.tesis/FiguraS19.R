#This script needs the environment deg


library(limma)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(patchwork)

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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS19")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS19"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS19")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

master <- read.table(paste0(wkdir,
"/master_regulators/viper/viper.results.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)




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


rownames(master) <- gsub("-","_",rownames(master))


genes <- c(
"PTPRN",
"RGS11",
"MAPK10",
"MC2R",
"STMN3",
"SSTR2",
"INSM1",
"PTPRN2",
"KCNA5",
"TNFSF10"
)



bxp <- deg.boxplot(as.data.frame(t(master)),
base,
"Tumor.Primario",
"Subtipo",
color.fill = c("gold", "cornflowerblue", "tomato1") ,
outdir = outdir,
flag = "master",
genes = genes
)

pdf(paste0(outdir, "/FiguraS19.pdf"),
width = 20,
height = 20)

bxp[[1]] + bxp[[2]] + 
bxp[[3]] + bxp[[4]] +
bxp[[5]] + bxp[[6]]  + 
bxp[[7]] + bxp[[8]]  + 
bxp[[9]] + bxp[[10]]  + plot_layout(ncol = 2, guide = "collect")

dev.off()






#######################################################################


deg1 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group1-Group2.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg2 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group1-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)

deg3 <- read.table(paste0(wkdir,
"/clustering/master/TopTable_master_Group2-Group3.txt"),
sep = "\t",
row.names = 1,
header = T
)
#Table filter

deg1_filter <- deg1[deg1$adj.P.Val < 0.1 & abs(deg1$logFC) > 1,]
deg2_filter <- deg2[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 2,]
deg3_filter <- deg3[deg3$adj.P.Val < 0.05 & abs(deg3$logFC) > 2,]





genes <- c(rownames(deg1_filter), rownames(deg2_filter)[1:100], rownames(deg3_filter)[1:100])
genes <- unique(genes)

deg.plot(as.data.frame(t(master)),
df,
"TUMOR.PRIMARIO",
"SUBTIPO",
color.fill = c("tomato1", "cornflowerblue", "gold") ,
outdir = outdir,
flag = "primarios",
genes = genes,
filter.table = F
)


deg.plot(as.data.frame(t(master)),
df,
"SUBTIPO",
"SUBTIPO",
color.fill = c("tomato1", "cornflowerblue", "gold") ,
outdir = outdir,
flag = "completos",
genes = genes
)


master <- read.table(
paste0(outdir, "/primarios/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
header = T
)



genes <- c(rownames(deg1_filter), rownames(deg2_filter)[1:100], rownames(deg3_filter)[1:100])
genes <- unique(genes)


ntissue <- lapply(genes,function(x) {
df <- master[master$gene ==x,]
tissue <- unique(df$TUMOR.PRIMARIO)
nt <- lapply(tissue,function(y){

    tdf <- df[df$TUMOR.PRIMARIO == y,]
    if(nrow(tdf[tdf$p.adj < 0.05,])>0){
        return(1)
    }else {
       return(0)
    }

return(sum(unlist(tdf)))

})

nt <- sum(unlist(nt))

return(c("gene" = x, "tissues" = nt, "median" = median(df$p.adj)))
})

ntissue <- do.call(rbind,ntissue)

master <- merge(master, ntissue, by = "gene")
master$tissues <- as.numeric(master$tissues)
master$median <- as.numeric(master$median)

master <- master[order(master$tissues, decreasing = T),]
master <- master[master$tissues >0,]


table(ntissue[,2])

write.table(master,
paste0(outdir, "/primarios/Primario_ordenado_tissues.txt"),
sep = "\t",
row.names = F,
col.names = T
)

master <- master[order(master$median, decreasing = F),]


write.table(master,
paste0(outdir, "/primarios/Primario_ordenado_median.txt"),
sep = "\t",
row.names = F,
col.names = T
)

master <- master[order(master$p.adj, decreasing = F),]

write.table(master,
paste0(outdir, "/primarios/Primario_ordenado_p.adj.txt"),
sep = "\t",
row.names = F,
col.names = T
)



