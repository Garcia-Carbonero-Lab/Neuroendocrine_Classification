library(ComplexHeatmap)
library(circlize)



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



if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura9")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura9"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura9")



base <- base[,c("Subtype",
"PRIMARY_TUMOR",
"ImmuneScore",
"TumorPurity")]


colnames(base) <- c("Subtipo",
"Tumor.Primario",
"Score.Inmune",
"Pureza.Tumoral")

base$Tumor.Primario <- as.factor(base$Tumor.Primario)
levels(base$Tumor.Primario) <- c("Colorrectal",
"Estómago",
"Pulmón", 
"Páncreas",
"I.Delgado")



base$Subtipo <- as.factor(base$Subtipo)
levels(base$Subtipo) <- c("SN1", "SN2", "SN3")

mcp <- read.csv(paste0(wkdir, "/immune/MCPcounter.tsv"),
                 sep = "\t",
                 row.names = 1,
                 header = T,
                 check.names = F)


rownames(mcp) <-c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Natural Killer",
"Células.B", 
"Linaje.monocítico",
"Células.dendríticas",
"Neutrófilos",
"Células.endoteliares",
"Fibroblastos")



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
    ),
    "Score.Inmune" =    colorRamp2(c(
                min(base$Score.Inmune, na.rm = T),
                median(base$Score.Inmune, na.rm = T),
                max(base$Score.Inmune, na.rm = T)),
            c("white", "green", "darkgreen")),
    "Pureza.Tumoral" = colorRamp2(
            c(
                min(base$Pureza.Tumoral, na.rm = T),
                median(base$Pureza.Tumoral, na.rm = T),
                max(base$Pureza.Tumoral, na.rm = T)),
            c("white", "chocolate1", "chocolate4")),
    
    "pvalor" = c("No" = "white",
    "Sí" = "black")
)



lgds1 <- lapply(names(l_colors)[1:2],function(x){

lgd <- Legend(labels = names(l_colors[[x]]),
legend_gp = gpar(fill = l_colors[[x]]),
title = x,
title_gp = gpar(fontsize = 80),
labels_gp = gpar(fontsize = 80),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm"),
title_position = "topcenter"
)
return(lgd)
})


lgds2 <- lapply(names(l_colors)[3:4],function(x){

lgd <- Legend(
col_fun = l_colors[[x]],
title = x,
title_gp = gpar(fontsize = 80),
labels_gp = gpar(fontsize = 80),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm"),
title_position = "topcenter"
)
return(lgd)
})

lgds3 <- Legend(labels = names(l_colors[[5]]),
legend_gp = gpar(fill = l_colors[[5]]),
title = "0.01<p-valor<0.05",
title_gp = gpar(fontsize = 80),
labels_gp = gpar(fontsize = 80),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1,
row_gap = unit(1.2, "cm"),
border = "black",
title_position = "topcenter"
)


pvalor <- data.frame(células = c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Natural Killer",
"Células.B", 
"Linaje.monocítico",
"Células.dendríticas",
"Neutrófilos",
"Células.endoteliares",
"Fibroblastos"),
"pvalor" = c(0.0196254030722943,
0.42755381710525,
0.0404252084740059,
0.0279371825873521,
0.0279371825873521,
0.0196254030722943,
0.0196254030722943,
0.0279371825873521,
0.42755381710525,
0.563799825448569))

rownames(pvalor) <- c("Células.T",
"T.CD8",
"Linfocitos.citotóxicos",
"Natural Killer",
"Células.B", 
"Linaje.monocítico",
"Células.dendríticas",
"Neutrófilos",
"Células.endoteliares",
"Fibroblastos")

pvalor[,2][pvalor[,2] > 0.05] <- "No"
pvalor[,2][pvalor[,2] != "No"] <- "Sí"


pvalor <- pvalor[,2,drop = F]


col <- list(
    "pvalor" = c("No" = "white",
    "Sí" = "black"))


ha_row <- rowAnnotation(
            df = pvalor,
            show_legend = F,
            show_annotation_name = c(base = F),
            col = col,
            annotation_name_gp = gpar(fontsize = 60),
            border = unit(60, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 60, col = "black"),
            annotation_name_side = "top")


data <- zscore.rows2(mcp)

anot <- base[, names(l_colors)[1:length(l_colors)-1]]

data <- data[,rownames(anot)]


ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = F,
            col = l_colors,
            annotation_name_gp = gpar(fontsize = 80),
            border = unit(40, "cm"),
            annotation_height = c(1.5,1,1,1),
            height = unit(20, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE)


hm1 <- Heatmap(as.matrix(data),
            top_annotation = ha_column,
            left_annotation= ha_row,
            show_row_names = T,
            show_column_names = F,
            col = colorRamp2(
                c(max(data), 0, min(data)),
                c("red", "white", "blue")
            ),
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            column_dend_side = "top",
            column_dend_height = unit(6, "cm"),
            cluster_columns = T,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "Z score (Infiltración Inmune)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 80),
                labels_gp = gpar(fontsize = 80,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            column_split = anot[, "Subtipo"],
            cluster_column_slices = FALSE,
            column_title  = c("","","")
        )


lgds <- append(lgds1, lgds2)
lgds <- append(lgds, lgds3)
lgds <- packLegend(lgds[[1]],lgds[[2]],lgds[[3]],
direction = "horizontal"
)

hlist <- hm1

p <- draw(hlist, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        padding = unit(c(20, 300, 20, 20), "mm"),
        annotation_legend_list = lgds)

pdf(paste0(
            outdir, "/Figura9.pdf"
        ),
        width = 60, 50
        )
        print(p)
        dev.off()






