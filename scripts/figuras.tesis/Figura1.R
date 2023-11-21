#environment figuras

library(MOFA2)
library(ggplot2)
#load functions
source("functions/MOFA/MOFA.functions.R")
source("functions/MOFA/adictional.functions.R")
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



mod <-  load_model(paste0(wkdir, "/MOFA/model/Mofa.model.hdf5"))




views_names(mod) <- c("Expresión", "Metilación")
groups_names(mod) <- c("Colorrectal",
"Estómago",
"Pulmón",
"Páncreas",
"I.Delgado")

mod <- subset_groups(mod, c("Colorrectal", "Estómago", "I.Delgado", "Páncreas", "Pulmón"))


factors_names(mod) <- gsub("Factor", "Factor ", factors_names(mod))


if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura1")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura1"))

}

pve <- plot_variance_explained(mod, x = 'group', y = 'factor')
pve <- theme_general(pve)

pve$guides$fill$title <- "Varianza %)"


ggsave(filename = paste0(wkdir, "/Figuras_tesis/Figura1/Figura1.pdf"),
device = "pdf"
,dpi =  300, width = 29, height = 15)



