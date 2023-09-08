#For this script we need  geneset.yml

library(limma)
library(GSEABase)
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

source("functions/genesets/genesets.functions.R")

#We create output folder
if (dir.exists(paste0(wkdir, "/genesets")) == F) {

    dir.create(paste0(wkdir, "/genesets"))

}

#open expression quality adjusted matrix

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/RMA_Expression_data_genes.tsv"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

hallmarks <- getGmt(paste0(wkdir, "/genesets/hallmarks.gmt"))
public <- getGmt(paste0(wkdir, "/genesets/Genesets_Subtypes.gmt"))
cells <- getGmt(paste0(wkdir, "/genesets/cells.gmt"))
gobp <- getGmt(paste0(wkdir, "/genesets/gobp.gmt"))
biocarta <- getGmt(paste0(wkdir, "/genesets/biocarta.gmt"))
kegg <- getGmt(paste0(wkdir, "/genesets/kegg.gmt"))
pid <- getGmt(paste0(wkdir, "/genesets/pid.gmt"))
reactome <- getGmt(paste0(wkdir, "/genesets/reactome.gmt"))


# Use the same gene symbol in genset and expression matrix
hallmarks <- anotate.gmt(hallmarks, expression)
public <- anotate.gmt(public, expression)
cells <- anotate.gmt(cells, expression)
gobp <- anotate.gmt(gobp, expression)
biocarta <- anotate.gmt(biocarta, expression)
kegg <- anotate.gmt(kegg, expression)
pid <- anotate.gmt(pid, expression)
reactome <- anotate.gmt(reactome, expression)




write.sampleset(public,paste0(wkdir, "/genesets/Genesets_Subtypes.gmt"))
write.sampleset(hallmarks,paste0(wkdir, "/genesets/hallmarks.gmt"))
write.sampleset(cells,paste0(wkdir, "/genesets/cells.gmt"))
write.sampleset(gobp,paste0(wkdir, "/genesets/gobp.gmt"))
write.sampleset(biocarta,paste0(wkdir, "/genesets/biocarta.gmt"))
write.sampleset(kegg,paste0(wkdir, "/genesets/kegg.gmt"))
write.sampleset(pid,paste0(wkdir, "/genesets/pid.gmt"))
write.sampleset(reactome,paste0(wkdir, "/genesets/reactome.gmt"))

