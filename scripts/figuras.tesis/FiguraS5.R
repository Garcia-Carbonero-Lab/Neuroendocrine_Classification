library(MOFA2)
library(patchwork)
library(ggplot2)

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


if (dir.exists(paste0(wkdir, "/Figuras_tesis/FiguraS5")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/FiguraS5"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/FiguraS5")

mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))



figs51 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor1",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")



figs52 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor2",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")

figs53 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor3",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")


figs54 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor4",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")


figs55 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor5",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")


figs56 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor6",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")

figs57 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor7",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")


figs58 <- plot_top_weights(mod,
                    view = "Expression",
                    factors = "Factor8",
                    nfeatures = 20, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 20)) +
                ylab("Peso")



pdf(paste0(outdir, "/figs5.pdf"), height = 30, width = 20)
figs51 + figs52 + figs53 + figs54 + 
figs55 + figs56 + figs57 + figs58 + plot_layout(ncol = 2)
dev.off()



s