#This script needs CNV environment

library(ChAMP)
# open functions

source("functions/preprocess/methylation.process.R")
source("functions/main.functions.R")

# Open config file
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

#open methylation data with quality adjustment
methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)

base <- base[colnames(methylation),]


# generate folders
if (dir.exists(paste0(wkdir, "/CNV")) == F) {

    dir.create(paste0(wkdir, "/CNV"))

}

#Segmnents calculation
cnal <- champ.CNA(intensity = methylation
                  ,pheno = base[,"PRIMARY_TUMOR"]
                  ,control = F
                  ,sampleCNA = T
                  ,groupFreqPlots = T
                  ,Rplot = T
                  ,PDFplot = T
                  ,freqThreshold = 0.3
                  ,resultsDir = paste0(wkdir, "/CNV/CNV_segments")
                  ,arraytype = 'EPIC')


# Put original samples code
cnadf <- do.call(rbind,cnal[[1]])
cnadf$ID <- gsub('X','',cnadf$ID)
cnadf$ID <- gsub('.qn','',cnadf$ID)


write.table(cnadf
            ,paste0(wkdir, "/CNV/CNV_segments/CNA_methylation.segments")
            , sep = '\t'
            ,row.names = F
            ,col.names = T
            ,quote = F)
