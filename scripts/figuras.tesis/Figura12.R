library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(reshape2)
library(ggradar)
library(scales)
library(fmsb)
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura12B")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura12B"))

}

source("functions/main.functions.R")
source("functions/classification/biomarkers.functions.R")
source("functions/cluster.analysis/additional.functions.R")

outdir <- paste0(wkdir, "/Figuras_tesis/Figura12B")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


base$Subtype <- as.factor(base$Subtype)
levels(base$Subtype) <- c("SN1", "SN2", "SN3")


base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)
levels(base$PRIMARY_TUMOR) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")


##############################################################################
#RADAR


exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)



master <- read.table(paste0(wkdir,
"/master_regulators/viper/viper.results.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)


genes <- c("CHGA","CHGB",
"PTPRN","BEX1", "NOVA1",
"GSTA1", "FBP1", "FABP1","PIGR", "GJD2",
"GRB10", "HEPACAM2")

masteres <- c(
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

scale <- function(x){
df <- (x-min(x))/(max(x)-min(x))
}


data.radar <- function(data, base, genes, group){

data <- t(apply(data,1,scale))

data <- data[genes,]

#scale 0-1



df <- base[,group, drop = F]

df <- merge(df,t(data),by = "row.names")
df <- df[,2:ncol(df)]

df <- lapply(genes,function(x){
 
    res <- aggregate(df[,x],by= list(df[,group]), FUN= median)
    colnames(res) <- c(group,x)
    return(res)

})

df <- do.call(cbind,df)

df <- df[,!duplicated(colnames(df))]

return(df)
}



genes.radar <- data.radar(exp, base, genes,"Subtype")
master.radar <- data.radar(master, base, genes = masteres,"Subtype")

lcols <- c("gold", "cornflowerblue", "tomato1")

ggradar(genes.radar,
        background.circle.colour = "white",
        gridline.min.linetype = 1,
        gridline.mid.linetype = 1,
        gridline.max.linetype = 1,
        group.colours = lcols)


ggsave(
paste0(outdir, "/Genes_radar.pdf"),
device = "pdf",
dpi = 300)



ggradar(master.radar,
        background.circle.colour = "white",
        gridline.min.linetype = 1,
        gridline.mid.linetype = 1,
        gridline.max.linetype = 1,
        group.colours = lcols)


ggsave(
paste0(outdir, "/Master_radar.pdf"),
device = "pdf",
dpi = 300)





t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}

lt.tomato1 <- t_col("tomato1", perc = 80, name = "lt.tomato1")
lt.cornflowerblue <- t_col("cornflowerblue", perc = 80, name = "lt.cornflowerblue")
lt.gold <- t_col("gold", perc = 80, name = "lt.gold")


# Color vector
colors_border=c("gold", "cornflowerblue", "tomato1")
colors_in=c(lt.gold, lt.cornflowerblue, lt.tomato1)




# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!

rownames(genes.radar) <- genes.radar$Subtype
genes.radar <- genes.radar[,-1]

genes.radar <- rbind(rep(0.8,ncol(genes.radar)) , rep(0.2,ncol(genes.radar)) , genes.radar)




pdf(paste0(outdir, "/Genes.radar2.pdf"))

# plot with default options:
radarchart( genes.radar  , axistype=1 , 
    #custom polygon
    pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0.2,0.8,0.1), cglwd=1,
    #custom labels
    vlcex=0.8 
    )
dev.off()



# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!

rownames(master.radar) <- master.radar$Subtype
master.radar <- master.radar[,-1]

master.radar <- rbind(rep(0.7,ncol(master.radar)) , rep(0.2,ncol(master.radar)) , master.radar)




pdf(paste0(outdir, "/Master.radar2.pdf"))

# plot with default options:
radarchart( master.radar  , axistype=1 , 
    #custom polygon
    pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0.2,0.7,0.1), cglwd=0.8,
    #custom labels
    vlcex=0.8 
    )
dev.off()

################################################################################



gsea1 <- gsea[gsea$group1 == "G1" & gsea$group2 == "G2",]
gsea2 <- gsea[gsea$group1 == "G1" & gsea$group2 == "G3",]
gsea3 <- gsea[gsea$group1 == "G2" & gsea$group2 == "G3",]


################################################################################
#Con P valores
#GSEA1
gsea1 <- gsea1[,c("pathway", "NES", "group1")]


gsea1 <- gsea1 %>%
pivot_wider(names_from = group1, values_from = NES) %>%
replace(is.na(.), 0) %>%
mutate(
    SN2 = G1,
    SN1 =  G1,
    SN2 = replace(SN2, SN2>0 , 0),
    SN1 = replace(SN1, SN1<0, 0),
    SN2 = abs(SN2)
) %>%
transmute(
    pathway,
    SN1,
    SN2
)

pathways_gsea1 <- pathways[!pathways %in% gsea1$pathway]

gsea1 <- rbind(gsea1,data.frame( "pathway"= pathways_gsea1, 
"SN1" = rep(0,length(pathways_gsea1)),
"SN2" = rep(0,length(pathways_gsea1))))

gsea1 <- as.data.frame(gsea1)
rownames(gsea1) <- gsea1$pathway
gsea1 <- gsea1[,2:3]

#GSEA2
gsea2 <- gsea2[,c("pathway", "NES", "group1")]


gsea2 <- gsea2 %>%
pivot_wider(names_from = group1, values_from = NES) %>%
replace(is.na(.), 0) %>%
mutate(
    SN3 = G1,
    SN1 =  G1,
    SN3 = replace(SN3, SN3>0 , 0),
    SN1 = replace(SN1, SN1<0, 0),
    SN3 = abs(SN3)
) %>%
transmute(
    pathway,
    SN1,
    SN3
)

pathways_gsea2 <- pathways[!pathways %in% gsea2$pathway]

gsea2 <- rbind(gsea2,data.frame( "pathway"= pathways_gsea2, 
"SN1" = rep(0,length(pathways_gsea2)),
"SN3" = rep(0,length(pathways_gsea2))))

gsea2 <- as.data.frame(gsea2)
rownames(gsea2) <- gsea2$pathway
gsea2 <- gsea2[,2:3]
#GSEA3
gsea3 <- gsea3[,c("pathway", "NES", "group1")]


gsea3 <- gsea3 %>%
pivot_wider(names_from = group1, values_from = NES) %>%
replace(is.na(.), 0) %>%
mutate(
    SN3 = G2,
    SN2 =  G2,
    SN3 = replace(SN3, SN3>0 , 0),
    SN2 = replace(SN2, SN2<0, 0),
    SN3 = abs(SN3)
) %>%
transmute(
    pathway,
    SN2,
    SN3
)

pathways_gsea3 <- pathways[!pathways %in% gsea3$pathway]

gsea3 <- rbind(gsea3,data.frame( "pathway"= pathways_gsea3, 
"SN2" = rep(0,length(pathways_gsea3)),
"SN3" = rep(0,length(pathways_gsea3))))

gsea3 <- as.data.frame(gsea3)
rownames(gsea3) <- gsea3$pathway
gsea3 <- gsea3[,2:3]
##################################################33
gsea1 <- gsea1[pathways,]
gsea2 <- gsea2[pathways,]
gsea3 <- gsea3[pathways,]


vias <- c("Biocarta IL2RB",
"Gobp Metabolismo Hormonal",
"Hallmark Respuesta Inflamatoria",
"Hallmark Interferón Gamma",
"Hallmark Metabolismo de Xenobióticos",
"Hallmark Fosforilación Oxidativa",
"Hallmark TNFA NFKB",
"Hallmark Targets MYC",
"Kegg Ribosoma",
"Kegg Metabolismo de Ácidos Grasos",
"Kegg Apoptosis",
"Kegg Receptor Células T",
"Kegg Receptor Células B",
"Kegg Señalización MAPK",
"Kegg Ciclo Celular",
"Pid Señalización RB1",
"Pid Señalización FOXM1",
"Reactome Oxidaciones Biológicas",
"Reactome Señalización por Interleucinas",
"Reactome Sistema Neuronal",
"Reactome Replicación de ADN"
)

rownames(gsea1) <- vias
rownames(gsea2) <- vias
rownames(gsea3) <- vias

#########################################################################
#Heatmap

ha1 <- HeatmapAnnotation("Subtipo" = c("SN1","SN2"),
show_legend = T,
col = list("Subtipo"= c("SN1"="gold", "SN2"= "cornflowerblue")),
annotation_name_gp = gpar(fontsize = 60),
            border = unit(4, "cm"),
            height = unit(4, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            simple_anno_size = unit(4, "cm"),
                        annotation_legend_param = list(
        title_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 3
    )


)


ha2 <- HeatmapAnnotation("Subtipo" = c("SN1","SN3"),
show_legend = T,
col = list("Subtipo"= c("SN1"="gold", "SN3"= "tomato1")),
annotation_name_gp = gpar(fontsize = 60),
            border = unit(4, "cm"),
            height = unit(4, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            simple_anno_size = unit(4, "cm"),
            annotation_legend_param = list(
        title_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 3
    )


)

ha3 <- HeatmapAnnotation("Subtipo" = c("SN2","SN3"),
show_legend = T,
col = list("Subtipo"= c("SN2"="cornflowerblue", "SN3"= "tomato1")),
annotation_name_gp = gpar(fontsize = 60),
            border = unit(4, "cm"),
            height = unit(4, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            simple_anno_size = unit(4, "cm"),
                        annotation_legend_param = list(
        title_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        labels_gp = gpar(fontsize = 50),
        grid_height = unit(4, "cm"),
        grid_width = unit(4, "cm"),
        gap = unit(100, "mm"),
        ncol = 3
    )


)




color <- colorRamp2(
                c(1, max(abs(gsea$NES))),
                c("white", "red")
)



hm1 <- Heatmap(as.matrix(gsea1),
            top_annotation = ha1,
            show_column_names = F,
            cluster_columns = F,
            col = color,
            gap = unit(5, "mm"),
            heatmap_legend_param = list(
                title = "-log p.valor",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            height = unit(40, "cm"),
            rect_gp = gpar(col = "black", lwd = 2),
            show_heatmap_legend = F,
            show_row_dend = F,
        )



hm2 <- Heatmap(as.matrix(gsea2),
            top_annotation = ha2,
            show_column_names = F,
            gap = unit(5, "mm"),
            col = color,
            cluster_columns = F,
            heatmap_legend_param = list(
                title = "NES",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            height = unit(40, "cm"),
            rect_gp = gpar(col = "black", lwd = 2),
            show_heatmap_legend = T,
            show_row_dend = F
        )



hm3 <- Heatmap(as.matrix(gsea3),
            top_annotation = ha3,
            show_column_names = F,
            col= color,
             cluster_columns = F,
            gap = unit(5, "mm"),
            heatmap_legend_param = list(
                title = "-log p.valor",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            height = unit(40, "cm"),
            row_names_gp = gpar(fontsize = 40),
            row_names_max_width = unit(100,'cm'),
            rect_gp = gpar(col = "black", lwd = 2),
            show_heatmap_legend = F,
            show_row_dend = F


        )

hlist = hm1 + hm2 + hm3


p <- draw(hlist, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        ht_gap = unit(3, "cm"))

pdf(paste0(
            outdir, "/GSEA_Figura12B.pdf"
        ),
        width = 40,
        height = 40
        )
        print(p)
        dev.off()
