#This script needs the environment figuras

library(VennDiagram)
library(EnhancedVolcano)


source("functions/main.functions.R")
source("functions/cluster.analysis/corr.met.exp.functions.R")

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


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS18")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS18"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS18")


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

deg3_filter <- rownames(deg3)[deg3$adj.P.Val < 0.1 & abs(deg3$logFC) > 1]
deg2_filter <- rownames(deg2)[deg2$adj.P.Val < 0.05 & abs(deg2$logFC) > 2]
deg1_filter <- rownames(deg1)[deg1$adj.P.Val < 0.05 & abs(deg1$logFC) > 2]


deg1_filter_100 <- deg1_filter[1:100]
deg2_filter_100 <- deg2_filter[1:100]


deg2_SN1 <- rownames(deg2)[deg2$adj.P.Val < 0.05 & deg2$logFC > 2]
deg1_SN1 <- rownames(deg1)[deg1$adj.P.Val < 0.05 & deg1$logFC > 2]

deg2_SN3 <- rownames(deg2)[deg2$adj.P.Val < 0.05 & deg2$logFC < -2]
deg1_SN2 <- rownames(deg1)[deg1$adj.P.Val < 0.05 & deg3$logFC > 2]

coincidentes <- intersect(deg2_filter_100, deg1_filter_100)

intersect(coincidentes,deg2_SN1)



intersect(coincidentes, deg2_SN3)
intersect(coincidentes, deg1_SN2)

exclusivos_deg2 <- setdiff(deg2_filter_100, deg1_filter_100)
exclusivos_deg3 <- setdiff(deg2_filter_100, deg1_filter_100)

intersect(deg2_SN3, exclusivos)

intersect(deg2_SN1, exclusivos)

intersect(deg3_SN2, exclusivos_deg3)

# Chart

venn.diagram(
        x = list(deg1_filter, deg2_filter, deg3_filter),
        category.names = c("SN1.vs.SN2" , "SN1.vs.SN3" , "SN2.vs.SN3"),
        filename = paste0(outdir,"/Diagrama_Venn.png"),
        output=TRUE,
        # Output features
        imagetype= "png" ,
        height = 600 , 
        width = 600 , 
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
        cat.dist = c(0.04, 0.12, 0.09),
        cat.fontfamily = "sans",
        rotation = 1,
        reverse = T
)

######################################################

v3 <- EnhancedVolcano(deg3,
lab = rownames(deg3),
x = "logFC",
y = "adj.P.Val",
title = "SN2 vs SN3",
FCcutoff = 1,
pCutoff = 0.1,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
legendLabSize = 25,
legendLabels=c('NS','CL2','q-valor',
      'q-valor y CL2'),
caption = "",
xlim = c(min(deg3[["logFC"]], na.rm = TRUE) - 0.5, max(deg3[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg3[["adj.P.Val"]]), na.rm = TRUE) + 1))


v3 <- v3 + xlab("Cambio Logarítmico 2") + 
ylab("-Log(q-valor)") +

theme(plot.title = element_text(hjust = 0.5,
vjust= 1, size = 40),
axis.title.x = element_text(size = 25),
axis.title.y = element_text(size = 25),
axis.text.x = element_text(size = 15),
axis.text.y = element_text(size = 15),)

ggsave(plot = v3,
filename = paste0(outdir, "/Volcano3.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)


v2 <- EnhancedVolcano(deg2,
lab = rownames(deg2),
x = "logFC",
y = "adj.P.Val",
title = "SN1 vs SN3",
FCcutoff = 2,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
legendLabSize = 25,
legendLabels=c('NS','CL2','q-valor',
      'q-valor y CL2'),
caption = "",

xlim = c(min(deg2[["logFC"]], na.rm = TRUE) - 0.5, max(deg2[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg2[["adj.P.Val"]]), na.rm = TRUE) + 1))

v2 <- v2 + xlab("Cambio Logarítmico 2") + 
ylab("-Log(q-valor)") +

theme(plot.title = element_text(hjust = 0.5,
vjust= 1, size = 40),
axis.title.x = element_text(size = 25),
axis.title.y = element_text(size = 25),
axis.text.x = element_text(size = 15),
axis.text.y = element_text(size = 15),)


ggsave(plot = v2,
filename = paste0(outdir, "/Volcano2.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)

v1 <- EnhancedVolcano(deg1,
lab = rownames(deg1),
x = "logFC",
y = "adj.P.Val",
title = "SN1 vs SN2",
FCcutoff = 2,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0,
legendLabSize = 25,
legendLabels=c('NS','CL2','q-valor',
      'q-valor y CL2'),
caption = "",
xlim = c(min(deg1[["logFC"]], na.rm = TRUE) - 0.5, max(deg1[["logFC"]], na.rm = TRUE) +
0.5),
ylim = c(0, max(-log10(deg1[["adj.P.Val"]]), na.rm = TRUE) + 1))

v1 <- v1 + xlab("Cambio Logarítmico 2") + 
ylab("-Log(q-valor)") +

theme(plot.title = element_text(hjust = 0.5,
vjust= 1, size = 40),
axis.title.x = element_text(size = 25),
axis.title.y = element_text(size = 25),
axis.text.x = element_text(size = 15),
axis.text.y = element_text(size = 15),)


ggsave(plot = v1,
filename = paste0(outdir, "/Volcano1.pdf"),
device = "pdf",
dpi = 500, width = 15,
height = 10
)


#############################################################