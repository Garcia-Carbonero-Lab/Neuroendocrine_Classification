#This script needs the environment figuras
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(limma)
#load functions

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")
source("functions/cluster.analysis/corr.met.exp.functions.R")

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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura13")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura13"))

}

outdir <- paste0(wkdir, "/Figuras_tesis/Figura13")


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)

base$PRIMARY_TUMOR <- as.factor(gsub("_", " ", base$PRIMARY_TUMOR))


base$Subtype <- as.factor(base$Subtype)
levels(base$Subtype) <- c("SN1", "SN2", "SN3")


base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)
levels(base$PRIMARY_TUMOR) <- c("COLORRECTAL",
"ESTÓMAGO", "PULMÓN", "PÁNCREAS", "INTESTINO")


df <- base %>%
    transmute(
        SUBTIPO = Subtype,
        TUMOR.PRIMARIO = PRIMARY_TUMOR,
      )

data <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)

dep1 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group1-Group2.txt"),
sep = "\t",
header = T,
row.names = 1
)

dep2 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group1-Group3.txt"),
sep = "\t",
header = T,
row.names = 1
)

dep3 <- read.table(paste0(wkdir,
"/clustering/dep/Significant_anotated_dep_Group2-Group3.txt"),
sep = "\t",
header = T,
row.names = 1
)


dep1_filter <- dep1[dep1$adj.P.Val < 0.05 & abs(dep1$logFC) > 2,]
dep2_filter <- dep2[dep2$adj.P.Val < 0.05 & abs(dep2$logFC) > 2,]
dep3_filter <- dep3[dep3$adj.P.Val < 0.05 & abs(dep3$logFC) > 2,]



dep1_filter_genes <- dep1_filter[complete.cases(dep1_filter$genes),]
dep2_filter_genes <- dep2_filter[complete.cases(dep2_filter$genes),]
dep3_filter_genes <- dep3_filter[complete.cases(dep3_filter$genes),]



write.table(dep1_filter,
paste0(outdir, "/Filter_Table_dep1.txt"),
sep= "\t",
row.names = T,
col.names = NA
)

write.table(dep2_filter,
paste0(outdir, "/Filter_Table_dep2.txt"),
sep= "\t",
row.names = T,
col.names = NA
)

write.table(dep3_filter,
paste0(outdir, "/Filter_Table_dep3.txt"),
sep= "\t",
row.names = T,
col.names = NA
)

venn.diagram(
        x = list(rownames(dep1_filter), rownames(dep2_filter), rownames(dep3_filter)),
        category.names = c("SN1.vs.SN2" , "SN1.vs.SN3" , "SN2.vs.SN3"),
        filename = paste0(outdir,"/Diagrama_Venn.png"),
        output=TRUE,
        # Output features
        imagetype= "png" ,
        height = 700 , 
        width = 700 , 
        resolution = 300,
        compression = "lzw",
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = c("green","orange","magenta"),
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(0, 0, 0),
        cat.dist = c(0.055, 0.055, -0.489),
        cat.fontfamily = "sans",
        rotation = 1
)

#################################################################################################

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



epic <- read.csv(paste0(datadir,
"/methylome/anotation/infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
                       sep = ",",
                       header = T,
                       skip = 7)

met_anot <- create.met.anot(epic)


exp <- exp[genes,]

df.corr <- corr.met.exp.unique(exp, data, met_anot)
df.corr$p.adj <- p.adjust(df.corr$p.val, method = "fdr")
df.corr <- merge(df.corr, met_anot, by.x = "CpG", by.y = "ID")
df.corr <- df.corr[order(df.corr$p.val),]

write.table(df.corr,
paste0(outdir,"/Correlation_exp_met.txt"),
sep = "\t",
row.names = T,
col.names = NA)

########################################################3
#bar plot of correlation
df.corr <- read.table(paste0(outdir, "/Correlation_exp_met.txt"),
sep = "\t",
header = T
)

# select the most correlation feature

df.corr <- df.corr[!duplicated(df.corr$gene),]

df.corr <- df.corr[abs(df.corr$corr) > 0.4,]
df.corr <- df.corr[-34,]
df.corr$p.adj <- -log10(df.corr$p.adj)


df.corr <- df.corr[order(df.corr$gene, decreasing = T),]
df.corr$gene <- factor(df.corr$gene, levels = df.corr$gene)


df.corr$p.adj

pdf(paste0(outdir, "/Figura13B.pdf"), width = 20, height = 20)
ggplot(df.corr) +
geom_col(aes(x = gene, y = corr, fill= p.adj), size = 1) +
geom_line(aes(x = gene, y = r2), size = 1.5,
color="darkslategray3", group = 1) +
scale_y_continuous(
    name = "Coeficiente de Correlación",
    sec.axis = sec_axis( trans=~.,
    name="Coeficiente de Determinación")
    ) +
scale_fill_gradient(low = "white", high = "red",
limits = c(3,max(df.corr$p.adj)),
                guide = guide_colourbar(
                barwidth = 2, barheight = 10,
                title.position = "top",
                title.hjust = 0.5),
                name = "-log10(q-value)") +
theme_classic() +
theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.ticks.length.y = unit(10,"pt"),
        axis.ticks.length.x = unit(10,"pt")) +

        xlab("")+
coord_flip() 


dev.off()

#corr per primary tumor

tissues <- unique(base$PRIMARY_TUMOR)

corr.tissue <- lapply(tissues, function(x) {
    df <- base[base$PRIMARY_TUMOR == x,]


    exp.df = exp[,rownames(df)]
    met.df = data[,rownames(df)]

    corr <- corr.met.exp.unique(exp.df, met.df, met_anot)
    corr$p.adj <- p.adjust(corr$p.val, method = "fdr")
    corr$tissue <- rep(x,nrow(corr))
    corr <- merge(corr, met_anot, by.x = "CpG", by.y = "ID")
    corr <- corr[order(corr$p.val),]
    return(corr)

})


corr.tissue.df <- as.data.frame(do.call(rbind, corr.tissue))

write.table(corr.tissue.df,
paste0(outdir,"/Correlation_exp_met_primary.txt"),
sep = "\t", row.names = F, col.names = T)


df.corr <- df.corr[abs(df.corr$corr) > 0.4,]


corr.tissue.df <- corr.tissue.df[corr.tissue.df$CpG %in% df.corr$CpG,]


write.table(corr.tissue.df,
paste0(outdir,"/Correlation_exp_met_primary_filtrado.txt"),
sep = "\t", row.names = F, col.names = T)



corr.tissue.df <- corr.tissue.df[corr.tissue.df$genes %in% c("CHGA","CHGB", "PTPRN"),]

write.table(corr.tissue.df,
paste0(outdir,"/Correlation_exp_met_primary_genes.txt"),
sep = "\t", row.names = F, col.names = T)


#corr per primary tumor

subtypes <- unique(base$Subtype)

corr.subtypes <- lapply(subtypes, function(x) {
    df <- base[base$Subtype == x,]


    exp.df = exp[,rownames(df)]
    met.df = data[,rownames(df)]

    corr <- corr.met.exp.unique(exp.df, met.df, met_anot)
    corr$p.adj <- p.adjust(corr$p.val, method = "fdr")
    corr$tissue <- rep(x,nrow(corr))
    corr <- merge(corr, met_anot, by.x = "CpG", by.y = "ID")
    corr <- corr[order(corr$p.val),]
    return(corr)

})


corr.subtypes.df <- as.data.frame(do.call(rbind, corr.subtypes))

write.table(corr.subtypes.df,
paste0(outdir,"/Correlation_exp_met_subtypes.txt"),
sep = "\t", row.names = F, col.names = T)


df.corr <- df.corr[abs(df.corr$corr) > 0.4,]


corr.subtypes.df <- corr.subtypes.df[corr.subtypes.df$CpG %in% df.corr$CpG,]


write.table(corr.subtypes.df,
paste0(outdir,"/Correlation_exp_met_subtypes_filtrado.txt"),
sep = "\t", row.names = F, col.names = T)



corr.subtypes.df <- corr.subtypes.df[corr.subtypes.df$genes %in% c("CHGA","CHGB", "PTPRN"),]

write.table(corr.subtypes.df,
paste0(outdir,"/Correlation_exp_met_subtypes_genes.txt"),
sep = "\t", row.names = F, col.names = T)











#We select the highest correlated cpg per gene
df.corr <- df.corr[!duplicated(df.corr$gene),]


#significant correlation
df.corr <- df.corr[df.corr$p.adj < 0.05,]
df.corr.rep <- df.corr



l_colors<- list(
    "SUBTIPO" = c("SN1" = "gold",
    "SN2" = "cornflowerblue",
    "SN3" = "tomato1"),
    "TUMOR.PRIMARIO" = c(
        "COLORRECTAL" = "orange",
        "PULMÓN" = "blue",
        "INTESTINO" = "yellow",
        "PÁNCREAS" = "green",
        "ESTÓMAGO" = "pink"
    )
)



anot_param_legend <- list(
    "SUBTIPO" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    
    ),
    "TUMOR.PRIMARIO" = list(
        title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 2
    )
)


symbol_genes <- alias2SymbolTable(rownames(exp))
rownames(exp)[is.na(symbol_genes) == F] <- symbol_genes[is.na(symbol_genes) == F]






met.corr <- data[df.corr$CpG,]
rownames(met.corr) <- df.corr$gene



heatmap.corr(exp = exp[df.corr$gene,],
met = met.corr,
df,
outdir = outdir,
col = l_colors,
size = list("width" = 50, "height" = 75),
anot_param_legend = anot_param_legend,
flag = "Met_Exp",
max.exp = NULL,
min.exp = NULL,
max.met = NULL,
min.met = NULL,
names.size = 3,
column = "SUBTIPO",
anot_h = c(40,20),
annotation_name_gp = 60,
show_row_names = T,
cluster = "exp",
order = c("SN1", "SN2", "SN3")
)
#############################################################################
#Correlations

correlation.plot(exp,
data[df.corr$CpG,],
met_anot,
df,
genes = c("CHGA", "CHGB", "PTPRN"),
outdir = outdir,
group= "SUBTIPO",
color.group = c("gold", "cornflowerblue","tomato1"),
flag = "SN"
)

correlation.plot(exp,
data[df.corr$CpG,],
met_anot,
df,
genes = c("CHGA", "CHGB", "PTPRN"),
outdir = outdir,
group= "TUMOR.PRIMARIO",
color.group = c("orange", "blue","yellow", "green", "pink"),
flag = "PRIMARIO"
)




df.corr$genes[df.corr$corr > 0.04]


probes <- met_anot$ID[met_anot$genes == "CHGA"]


