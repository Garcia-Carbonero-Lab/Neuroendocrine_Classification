# IMPORTANT NOTES
# For this script we ned MOFA.yml enviroment



# We open MOFA library
library(MOFA2)

# we read config file
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
source("functions/MOFA/MOFA.functions.R")

#open selected expression features
exp <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression_fs.txt"),
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = F
)
#open selected methylation features

met <- read.table(paste0(wkdir,
"/preprocess/methylome/methylation_fs.txt"),
    sep = "\t",
    row.names = 1,
    header = T,
    check.names = F
)

base <- read.csv(paste0(datadir, "/clinical_data.txt"),
    sep = "\t",
    row.names = 1,
    header = T
)

#order clinical data in the same way tht expression
base <- base[colnames(exp), ]


data.list <- list(
    "Expression" = as.matrix(exp),
    "Methylation" = as.matrix(met[, colnames(exp)])
)


if (dir.exists(paste0(wkdir, "/MOFA/model")) == F) {
    dir.create(paste0(wkdir, "/MOFA/model"))
}

#Crete the MOFa model with 8 factors
Mofa.model(
    input = data.list, file = paste0(
        wkdir,
        "/MOFA/model/Mofa.model.hdf5"
    ),
    groups = base$PRIMARY_TUMOR,
    num_factors = 8,
    convergence_mode = "slow"
)

