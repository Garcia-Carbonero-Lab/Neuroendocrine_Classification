#this script need the environment figures

library(ggplot2)
library(tidyverse)


#load functions
source("functions/MOFA/MOFA.functions.R")
source("functions/MOFA/adictional.functions.R")
source("functions/main.functions.R")

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

gsea_cell <- read.table(paste0(wkdir,
"/MOFA/analysis_results/features/Factor2/cells_GSEA_expression.txt")
, sep = "\t",
header = T
)



if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura4")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura4"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura4")

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

rownames(gsea_cell) <- gsea_cell$pathway


gsea_cell <- gsea_cell[cells,]



gsea_cell$pathway <- c("Células Enterocromafines Duodeno",
"Células Copa Duodeno",
"Células I Duodeno",
"Células K Duodeno",
"Enterocitos Maduros Duodeno",
"Células Paneth Duodeno",
"Células Madre Duodeno",
"Células de tránsito Duodenales",
"Células Madre en Proceso De Diferenciación Duodenales",
"Células Gástricas D",
"Células Gástricas G",
"Células Gástricas del Istmo",
"Células Gástricas LYZ Positivas",
"Células Gástricas PIT Maduras",
"Células Gástricas con Metaloteina",
"Células Gástricas de Cuello",
"Células Gástricas Similares a las Enterocormafines",
"Células Gástricas Parietales",
"Células Gástricas X",
"Células Intestino Delgado Enterocromafines",
"Células Intestino Delgado MKI67 alto",
"Células Intestino Delgado Mesenquimales",
"Enterocitos Intestinales",
"Células Intestino Paneth",
"Células Intestino de Copa",
"Células Páncreas Delta",
"Células Páncreas Endoteliales",
"Células Páncreas Epsilon",
"Células Páncreas Polipeptídicas",
"Células Epiteliales Alveorales I",
"Células Epiteliales Alveolares II",
"Células Pulmón Basales",
"Células Pulmón Club",
"Células Pulmón Copa",
"Células Pulmón Mucosas",
"Células Pulmón Neuroendocrinas"
)


gsea_cell$pathway <- factor(gsea_cell$pathway,
levels= gsea_cell$pathway[order(gsea_cell$NES)])



plt <- ggplot(data = gsea_cell, aes(x = NES, y = pathway, size = size,
        color = -1 * log10(padj))) +
            geom_point() +
            scale_size_area(max_size = 10) +
            scale_colour_gradient(low = "green", high = "red") +
            theme(text = element_text(size = 26),
                   legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(margin = margin(10, 10, 10, 20),
        colour = "black"),
        axis.ticks.length.y = unit(10,"pt"),
        axis.ticks.length.x = unit(10,"pt"))


plt$labels$y   <- ""
plt$labels$size   <- "Tamaño"





        ggsave(plt,
            filename = paste0(outdir,"/Dotplot_cells.pdf"),
            device = "pdf", dpi = 500, width = 22,
        height = 12
        )

write.table(gsea_cell,
paste0(outdir,"/GSEA_cells.txt"),
sep = "\t",
row.names= T,
col.names = NA )
