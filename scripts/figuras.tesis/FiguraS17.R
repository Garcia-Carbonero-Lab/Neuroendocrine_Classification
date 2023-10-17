library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(reshape2)
library(patchwork)
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

if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS17")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS17"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS17")

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


gsea_cell1 <- read.table(
paste0(wkdir, "/clustering/gsea/cells_fgsea_G1vsG2.txt"),
sep = "\t",
header = T,
row.names = 1
)

gsea_cell2 <- read.table(
paste0(wkdir, "/clustering/gsea/cells_fgsea_G1vsG3.txt"),
sep = "\t",
header = T,
row.names = 1
)

gsea_cell3 <- read.table(
paste0(wkdir, "/clustering/gsea/cells_fgsea_G2vsG3.txt"),
sep = "\t",
header = T,
row.names = 1
)


cells <- c("BUSSLINGER_DUODENAL_EC_CELLS",
"BUSSLINGER_DUODENAL_GOBLET_CELLS",
"BUSSLINGER_DUODENAL_I_CELLS",
"BUSSLINGER_DUODENAL_K_CELLS",
"BUSSLINGER_DUODENAL_MATURE_ENTEROCYTES",
"BUSSLINGER_DUODENAL_PANETH_CELLS",
"BUSSLINGER_DUODENAL_STEM_CELLS",
"BUSSLINGER_DUODENAL_TRANSIT_AMPLIFYING_CELLS",
"BUSSLINGER_DUODENAL_DIFFERENTIATING_STEM_CELLS",
"BUSSLINGER_GASTRIC_D_CELLS",
"BUSSLINGER_GASTRIC_G_CELLS",
"BUSSLINGER_GASTRIC_ISTHMUS_CELLS",
"BUSSLINGER_GASTRIC_LYZ_POSITIVE_CELLS",
"BUSSLINGER_GASTRIC_MATURE_PIT_CELLS",
"BUSSLINGER_GASTRIC_METALLOTHIONEIN_CELLS",
"BUSSLINGER_GASTRIC_NECK_CELLS",
"BUSSLINGER_GASTRIC_OXYNTIC_ENTEROCHROMAFFIN_LIKE_CELLS",
"BUSSLINGER_GASTRIC_PARIETAL_CELLS",
"BUSSLINGER_GASTRIC_X_CELLS",
"GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS",
"GAO_LARGE_INTESTINE_ADULT_CH_MKI67HIGH_CELLS",
"GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS",
"GAO_LARGE_INTESTINE_24W_C10_ENTEROCYTE",
"GAO_LARGE_INTESTINE_24W_C11_PANETH_LIKE_CELL",
"GAO_LARGE_INTESTINE_24W_C8_GOBLET_CELL",
"MURARO_PANCREAS_DELTA_CELL",
"MURARO_PANCREAS_ENDOTHELIAL_CELL",
"MURARO_PANCREAS_EPSILON_CELL",
"MURARO_PANCREAS_PANCREATIC_POLYPEPTIDE_CELL",
"TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_1_CELL",
"TRAVAGLINI_LUNG_ALVEOLAR_EPITHELIAL_TYPE_2_CELL",
"TRAVAGLINI_LUNG_BASAL_CELL",
"TRAVAGLINI_LUNG_CLUB_CELL",
"TRAVAGLINI_LUNG_GOBLET_CELL",
"TRAVAGLINI_LUNG_MUCOUS_CELL",
"TRAVAGLINI_LUNG_NEUROENDOCRINE_CELL"
)




vias <- c(" Enterocromafines Duodeno",
"Copa Duodeno",
"I Duodeno",
"K Duodeno",
"Enterocitos Maduros Duodeno",
"Células Paneth Duodeno",
"Células Madre Duodeno",
"Células de tránsito Duodenales",
"Células Madre en Proceso De Diferenciación Duodenales",
"Gástricas D",
"Gástricas G",
"Gástricas del Istmo",
"Gástricas LYZ Positivas",
"Gástricas PIT Maduras",
"Gástricas con Metaloteina",
"Gástricas de Cuello",
"Gástricas Similares a las Enterocormafines",
"Gástricas Parietales",
"Gástricas X",
"Intestino Delgado Enterocromafines",
"Intestino Delgado MKI67 alto",
"Intestino Delgado Mesenquimales",
"Enterocitos Intestinales",
"Intestino Paneth",
"Intestino de Copa",
"Páncreas Delta",
"Páncreas Endoteliales",
"Páncreas Epsilon",
"Páncreas Polipeptídicas",
"Epiteliales Alveorales I",
"Epiteliales Alveolares II",
"Pulmón Basales",
"Pulmón Club",
"Pulmón Copa",
"Pulmón Mucosas",
"Pulmón Neuroendocrinas"
)



plot_cells <- function(data,cells,vias,flag,outdir){
    data <- data[cells,]
    data$pathway <- vias
    data <- data[complete.cases(data$NES),]
    data$pathway <- factor(data$pathway,
    levels= data$pathway[order(data$NES)])
    data <- data[order(data$NES),]
    
    plt <- ggplot(data = data, aes(x = NES, y = pathway, size = size,
        color = -1 * log10(padj))) +
            
            labs(title = paste0(flag)) +

            geom_point() +
            scale_size_area(max_size = 10) +
            scale_colour_gradient(low = "green", high = "red", limits = c(1.3, 48)) +
            theme(text = element_text(size = 26),
            legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(margin = margin(10, 10, 10, 20),
        colour = "black"),
        axis.ticks.length.y = unit(10,"pt"),
        axis.ticks.length.x = unit(10,"pt"),
            plot.title = element_text(hjust = 0.5,
            vjust= 1, size = 40)) 

    plt$labels$y   <- ""
    plt$labels$size   <- "Tamaño"


        ggsave(plt,
            filename = paste0(outdir,"/Dotplot_cells",flag,".pdf"),
            device = "pdf", dpi = 500, width = 22,
        height = 12
        )
        return(plt)
}

gsea1 <- plot_cells(gsea_cell1,
cells = cells,
vias = vias,
flag = "SN1vsSN2",
outdir = outdir
)


gsea2 <- plot_cells(gsea_cell2,
cells = cells,
vias = vias,
flag = "SN1vsSN3",
outdir = outdir
)


gsea3 <- plot_cells(gsea_cell3,
cells = cells,
vias = vias,
flag = "SN2vsSN3",
outdir = outdir
)


pdf(paste0(outdir, "/FiguraS18.pdf"),
width = 20,
height = 30)

gsea1 /gsea2 / gsea3 + plot_layout(guides = "collect")

dev.off()


