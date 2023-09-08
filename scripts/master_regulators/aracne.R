
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

# Create folder for results

if (dir.exists(paste0(wkdir, "/master_regulators")) == F) {
  dir.create(paste0(wkdir, "/master_regulators"))
  dir.create(paste0(wkdir, "/master_regulators/aracne"))
}

# INSTALL ARACNE
setwd(paste0(datadir, "/master_regulators"))

system("git clone https://github.com/califano-lab/ARACNe-AP.git")

setwd("ARACNe-AP")

system("ant main")

#open expression quality adjusted matrix

expression <- read.table(paste0(wkdir,
"/preprocess/transcriptome/expression.qual.txt"),
  sep = "\t",
  header = T,
  check.names = F
)


#prepare matrix for aracne
colnames(expression) <- c("gene", colnames(expression[2:ncol(expression)]))


write.table(expression,
  paste0(datadir, "/master_regulators/expression_aracne.tsv"),
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

reg <- list.files((paste0(datadir, "/master_regulators")))
reg <- reg[grep("hugo", reg)]

for (r in reg) {
  name <- gsub("\\.txt", "", r)
  # we calculate thresholds
  system(paste0(
    "java -Xmx5G ",
    "-jar dist/aracne.jar ",
    "-e ", datadir, "/master_regulators/expression_aracne.tsv ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --tfs ", datadir, "/master_regulators/", r,
    " --pvalue 1E-8 --seed 1 --calculateThreshold"
  ))

  # We apply boostraping
  system(paste0(
    "for i in $(seq 1 100) ; do ",
    "java -Xmx5G ",
    "-jar dist/aracne.jar ",
    "-e ", datadir, "/master_regulators/expression_aracne.tsv ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --tfs ", datadir, "/master_regulators/", r,
    " --pvalue 1E-8 --seed $i; done"
  ))
  
# We consolidate the 100 bostraping in one network

  system(paste0(
    "java -Xmx10G -jar dist/aracne.jar ",
    "-o ", wkdir, "/master_regulators/aracne/", name,
    " --consolidate"
  ))

#write result as tsv format


  aracne <- read.table(paste0(wkdir, "/master_regulators/aracne/",
  name, "/network.txt")
  , sep = "\t"
  , header = T)

  aracne <- aracne[complete.cases(aracne[, 3]),]

  write.table(aracne, paste0(wkdir, "/master_regulators/aracne/",
  name, "/network.tsv"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F)

}

#Concatenate all files
system(paste0("cat ", wkdir, "/master_regulators/aracne/*/network.tsv > ",
wkdir, "/master_regulators/aracne/concatenated_networks.tsv"))
