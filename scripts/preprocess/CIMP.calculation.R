#this script needs champ environment
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

# We create output folder
if (dir.exists(paste0(wkdir, "/MOFA")) == F) {
    dir.create(paste0(wkdir, "/MOFA"))
}




# open functions from MOFA file
source("functions/Preprocess/methylation.process.R")



clinical.data <- read.table(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    header = T
)

methylation <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F

)


methylation <- apply(methylation, c(1,2), as.numeric)

#Calculation of CIMP

CIMP <- addCIMPstatus(methylation)


#Add CIMP to clinical data
clinical.data <- merge(clinical.data, CIMP,
    by.x = "CLINICAL_TRIAL_CODE",
    by.y = "row.names"
)


write.table(clinical.data, paste0(datadir,"/clinical_data.txt"),
    sep = "\t", row.names = F, col.names = T
)
