
#necesita el ambiente paneles
library(patchwork)
library(figpatch)
library(magick)
library(rsvg)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
#load functions
source("functions/main.functions.R")
source("functions/figuras/figuras.functions.R")
source("functions/classification/biomarkers.functions.R")

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

outdir <- paste0(wkdir, "/Figuras_tesis/")

######################################################################3
#Figura2 # NO

f2a <- fig(paste0(outdir, "Figura2/Heatmap_anotated_Kruskal.clinica.pdf"))
f2b <- fig(paste0(outdir, "Figura2/overall_survival_Hazard_ratio_factors.pdf"))
f2c <- fig(paste0(outdir, "Figura2/mcp_immune.linear.model.pdf"))
f2d <- fig(paste0(outdir, "Figura2/Heatmap_Figura_2_D.pdf"))

patchwork <- (f2a + f2b) / (f2c + f2d)

pdf(paste0(outdir, "/Figura2/Figura2.pdf"))
f2a + f2b + f2c + f2d +
plot_annotation(tag_levels = 'A') + 
plot_layout(design = c(area(1,1),
area(1,2), area(2,1,3), area(2,2,3)))
dev.off()


pdf(paste0(outdir, "/Figura2/Figura2.pdf"))
f2a + f2b + f2c + f2d +
plot_annotation(tag_levels = 'A') + 
plot_layout(nrow = 2 , heights = c(1,2))
dev.off()

#################################################################3333333
#Figura 6


f6a <- fig(paste0(outdir, "Figura6/Figura6A_corregida.pdf"))
f6b <- fig(paste0(outdir, "Figura6/Figura6B_corregida.pdf"))


pdf(paste0(outdir, "Figura6/Figura6.pdf"))
f6a /f6b +  plot_annotation(tag_levels = 'A')
dev.off()
###############################################
#Figura 8

#Figura8B

f8b1 <- fig(paste0(outdir, "Figura8/Fig8B1_corregida.pdf"))


df <- data.frame("Pacientes" = c(58, 69, 67, 194),
" " = c(51, 57, 54, 162), " " = c(31, 35, 32, 98),
" " = c(10, 13, 15, 38),
"Pr SG 5 años" = c(0.91, 0.88, 0.80, 0.86))


df1 <- data.frame("Pacientes" = c(67, 68, 58, 193),
" " = c(53, 55, 49, 157), " " = c(27, 31, 24, 82),
" " = c(15, 8, 6, 29))

rownames(df) <- c("SN1", "SN2", "SN3", "Total")

pdf(paste0(outdir,"/Figura8/Figura.8B2.pdf"))
g <- tableGrob(df, theme = ttheme_minimal())
g <- gtable_add_grob(g, 
grobs = segmentsGrob( # line across the bottom
            x0 = unit(0,"npc"),
            y0 = unit(0,"npc"),
            x1 = unit(1,"npc"),
            y1 = unit(0,"npc"),
            gp = gpar(lwd = 2.0)),
        t = 4, b = 4, l = 1, r = 7)

#g <- gtable_add_grob(g, 
#grobs = segmentsGrob( # line across the bottom
#            x0 = unit(0,"npc"),
#            y0 = unit(0,"npc"),
#            x1 = unit(0,"npc"),
#            y1 = unit(1,"npc"),
#            gp = gpar(lwd = 2.0)),
#        t = 1, b = 4, l = 6, r = 6)


f8b1 / g

dev.off()


fig8a<- fig(paste0(outdir, "Figura8/Figura8A_corregida.pdf"))
fig8b<- fig(paste0(outdir, "Figura8/Figura8B2completa.pdf"))
fig8c<- fig(paste0(outdir, "Figura8/Figura8C_completa.pdf"))


pdf(paste0(outdir,"/Figura8/Figura.8.pdf"))

(fig8a | (fig8c / fig8b + plot_layout(heights = c(1,2)))) + 
plot_annotation(tag_levels = 'A')

dev.off()




############################################################
fig8a<- fig(paste0(outdir, "Figura8/Heatmap_anotated_Figura8v2.pdf"))
fig8b<- fig(paste0(outdir, "Figura8/Figura8B2completa.pdf"))
fig8c<- fig(paste0(outdir, "Figura8/Figura8Cv2.pdf"))


pdf(paste0(outdir,"/Figura8/Figura8v2.pdf"))

((fig8a | fig8c) / fig8b) + 
plot_annotation(tag_levels = 'A')

dev.off()

########################################################
#Figura 12


fig12a <- fig(paste0(outdir, "Figura12/Genes.radar2.pdf"))
fig12b <- fig(paste0(outdir, "Figura12/Master.radar2.pdf"))
fig12c <- fig(paste0(outdir, "Figura12/Figura12B.pdf"))


pdf(paste0(outdir,"/Figura12/Figura12.pdf"))

((fig12a | fig12b) / fig12c) + 
plot_annotation(tag_levels = 'A')

dev.off()



fig12a <- fig(paste0(outdir, "Figura12B/Genes.radar2.pdf"))
fig12b <- fig(paste0(outdir, "Figura12B/Master.radar2.pdf"))
fig12c <- fig(paste0(outdir, "Figura12/Figura12B.pdf"))


pdf(paste0(outdir,"/Figura12B/Figura12.pdf"))

((fig12a | fig12b) / fig12c) + 
plot_annotation(tag_levels = 'A')

dev.off()

####################################################3
#Figura13
fig13a <- fig(paste0(outdir, "Figura13/Diagrama_Venn.png"))
fig13b <- fig(paste0(outdir, "Figura13/Figura_13B_corregida.pdf"))


pdf(paste0(outdir,"Figura13/Figura13.pdf"), width = 10, height = 5)
(fig13a | fig13b) + 
plot_annotation(tag_levels = 'A')
dev.off()

#####################################33
#figS15

fig15a <- fig(paste0(outdir, "FiguraS15/Volcano1.pdf"))
fig15b <- fig(paste0(outdir, "FiguraS15/Volcano2.pdf"))
fig15c <- fig(paste0(outdir, "FiguraS15/Volcano3.pdf"))
fig15d <- fig(paste0(outdir, "FiguraS15/Diagrama_Venn.png"))

pdf(paste0(outdir,"FiguraS15/FiguraS15.1.pdf"), width = 20, height = 10)
(fig15a | fig15b | fig15c) / fig15d +  plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2,1))
dev.off()


pdf(paste0(outdir,"FiguraS15/FiguraS15.2.pdf"), width = 10, height = 10)
fig15a + fig15b + fig15c + fig15d +  plot_annotation(tag_levels = 'A') 
dev.off()

#####################################33
#figS18

figs18a <- fig(paste0(outdir, "FiguraS18/Volcano1.pdf"))
figs18b <- fig(paste0(outdir, "FiguraS18/Volcano2.pdf"))
figs18c <- fig(paste0(outdir, "FiguraS18/Volcano3.pdf"))
figs18d <- fig(paste0(outdir, "FiguraS18/Diagrama_Venn.png"))

pdf(paste0(outdir,"FiguraS18/FiguraS18.1.pdf"), width = 20, height = 10)
(figs18a | figs18b | figs18c) / figs18d +  plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2,1))
dev.off()


pdf(paste0(outdir,"FiguraS18/FiguraS18.2.pdf"), width = 10, height = 10)
figs18a + figs18b + figs18c + figs18d +  plot_annotation(tag_levels = 'A') 
dev.off()
