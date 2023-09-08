#This script need environment immune

library(MCPcounter)
library(xCell)

config <- read.csv("config/config.tsv",
    sep = "\t",
    header = T
)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]


clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    header = T
)

#open expression quality adjusted matrix
expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F)


if (dir.exists(paste0(wkdir, "/immune")) == F) {

    dir.create(paste0(wkdir, "/immune"))

}

#Caclutate MCPcounter scores
mcp <- MCPcounter.estimate(expression, featuresType = "HUGO_symbols")

rownames(mcp) <- gsub(" ",".",rownames(mcp))

write.table(mcp,
paste0(wkdir, "/immune/MCPcounter.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

