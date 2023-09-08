#This script need environment immunes

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


expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F)


if (dir.exists(paste0(wkdir, "/immune")) == F) {

    dir.create(paste0(wkdir, "/immune"))

}

#MCPcounter
mcp <- MCPcounter.estimate(expression, featuresType = "HUGO_symbols")
rownames(mcp) <- paste0(rownames(mcp), ".mcp")
rownames(mcp) <- gsub(" ",".",rownames(mcp))

write.table(mcp,
paste0(wkdir, "/immune/MCPcounter.tsv"),
sep = "\t",
row.names = T,
col.names = NA)

#xCell

xcell <- xCellAnalysis(expression)
rownames(xcell) <- paste0(rownames(xcell), ".xcell")
rownames(xcell) <- gsub(" ",".",rownames(xcell))
rownames(xcell) <- gsub("-",".",rownames(xcell))
rownames(xcell) <- gsub("\\+","",rownames(xcell))

write.table(xcell,
paste0(wkdir, "/immune/xCell.tsv"),
sep = "\t",
row.names = T,
col.names = NA)
