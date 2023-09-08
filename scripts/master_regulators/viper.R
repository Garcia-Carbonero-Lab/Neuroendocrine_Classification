#To use this sicript you need the environment viper
library(viper)

source("functions/master_regulators/viper.functions.R")

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

#open expression quality adjusted matrix

expression <- read.table(paste0(datadir, "/transcriptome/expression.qual.txt"),
sep = '\t',
header = T,
row.names = 1,
check.names = F
)

paste0(datadir, "/transcriptome/expression.qual.txt")



if (dir.exists(paste0(wkdir, "/master_regulators/viper")) == F) {
  dir.create(paste0(wkdir, "/master_regulators/viper"))
}


RegProcess(paste0(wkdir, "/master_regulators/aracne/concatenated_networks.tsv"),
as.matrix(expression),
out.dir = paste0(wkdir, "/master_regulators/viper/"),
out.name = "neuroendocrine")


net <- readRDS(paste0(wkdir, "/master_regulators/viper/neuroendocrinepruned.rds"))
viper.res <- viper(as.matrix(expression), net, method = "rank")

write.table(viper.res,
paste0(wkdir, "/master_regulators/viper/viper.results.txt"),
sep = "\t",
row.names = T,
col.names = NA,
quote = F
)
