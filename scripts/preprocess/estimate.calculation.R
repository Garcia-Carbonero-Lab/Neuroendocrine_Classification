#For this script we need the enviroment estimate
library(estimate)

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
check.names = F
)


write.table(expression,
paste0(wkdir,
"/preprocess/transcriptome/expression.estimate.txt"),
row.names = T,
col.names = T
, quote = F
, sep = "\t")


#Select genes of estimate
filterCommonGenes(input.f= paste0(wkdir,
"/preprocess/transcriptome/expression.estimate.txt"),
output.f= paste0(wkdir,
"/preprocess/transcriptome/estimate_genes.txt"),
id="GeneSymbol")

#Calculate estimate score
estimateScore(paste0(wkdir,
"/preprocess/transcriptome/estimate_genes.txt"), 
paste0(wkdir,
"/preprocess/transcriptome/estimate_scores.txt"),
platform="affymetrix")


estimate <- read.table(paste0(wkdir,
"/preprocess/transcriptome/estimate_scores.txt"),
sep = "\t",
 row.names = 1, 
 header =  T,
 skip = 2)

# Rename samples code
 colnames(estimate) <- gsub("X", "", colnames(estimate))
 colnames(estimate) <- gsub("\\.", " ", colnames(estimate))

estimate <- t(estimate)
estimate <- estimate[-1,]

#add estimate results into clinical data
clinical.data <- merge(clinical.data,
estimate,
by.x = "CLINICAL_TRIAL_CODE",
by.y = "row.names")

write.table(clinical.data, paste0(datadir,"/clinical_data.txt"),
    sep = "\t", row.names = F, col.names = T
)
