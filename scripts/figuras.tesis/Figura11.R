#This script needs the environment cnv.plot

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(maftools)

source("functions/main.functions.R")
source("functions/cluster.analysis/corr.met.exp.functions.R")

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


base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 header = T,
                 row.names = 1)






if (dir.exists(paste0(wkdir, "/Figuras_tesis/Figura11")) == F) {

    dir.create(paste0(wkdir, "/Figuras_tesis/Figura11"))

}


outdir <- paste0(wkdir, "/Figuras_tesis/Figura11")

cnvdis <- read.table(paste0(wkdir, "/CNV/preprocess/CNV.peaks.dis.txt"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)



fisher <- read.table(paste0(wkdir, "/clustering/cnv/CNV_discrete_table_fisher.txt"),
sep = "\t",
row.names = 1,
header = T,
check.names = F)

fisher <- fisher[which(fisher$P.value < 0.05),]
significant <- gsub("\\.Altered","",rownames(fisher))


########################################################
all.lesions <- paste0(wkdir, "/CNV/GISTIC/all_lesions.conf_99.txt")

amp.genes <- paste0(wkdir, "/CNV/GISTIC/amp_genes.conf_99.txt")


del.genes <- paste0(wkdir, "/CNV/GISTIC/del_genes.conf_99.txt")


scores.gis <- paste0(wkdir, "/CNV/GISTIC/scores.gistic")


laml.gistic = readGistic(gisticAllLesionsFile = all.lesions,
gisticAmpGenesFile = amp.genes,
gisticDelGenesFile = del.genes,
gisticScoresFile = scores.gis,
isTCGA = FALSE)


base$Tumor_Sample_Barcode = rownames(base)

base$Subtype <- as.factor(base$Subtype)
levels(base$Subtype) <- c("SN1", "SN2", "SN3")


base$PRIMARY_TUMOR <- as.factor(base$PRIMARY_TUMOR)
levels(base$PRIMARY_TUMOR) <- c("Colorrectal",
"Estómago", "Pulmón", "Páncreas", "I.Delgado")


df <- base %>%
    transmute(
        Subtipo = Subtype,
        Tumor.Primario = PRIMARY_TUMOR,
       Tumor_Sample_Barcode
      )

df$Tumor_Sample_Barcode <- gsub(" ",".", df$Tumor_Sample_Barcode)


sampleOrder <- c("1427", "903",
"1435", "1401", "BP36", "BP80", "206", "921", "BP54",
"125.NENs", "BP60", "103.NENs", "105.NENs", "11.NENs",
"111.NENs","112.NENs", "115.NENs", "121.NENs", "124.NENs",
"1305", "136.NENs", "137.NENs", "138.NENs",
"140.NENs", "142.NENs", "1420", "150.NENs", "154.NENs",
"155.NENs", "21.NENs", "24.NENs", "26.NENs", "2802", "35.NENs",
"36.NENs", "39.NENs", "40.NENs", "41.NENs", "48.NENs",
"51.NENs", "58.NENs", "65.NENs", "70.NENs", "76.NENs",
"77.NENs", "823", "828", "83.NENs", "84.NENs", "87.NENs",
"89.NENs", "908", "913", "93.NENs", "94.NENs",
"BP40", "BP44", "BP86","510","123.NENs",
"BP77", "1603", "67.NENs", "805", "835", "BP49",
"829", "220", "BP64", "BP50", "203", "1013", "1430",
"73.NENs", "BP53", "906", "838", "307", "926","BP51",
"116.NENs", "102.NENs", "1021", "104.NENs", "107.NENs",
"108.NENs", "110.NENs", "117.NENs", "119.NENs", "122.NENs",
"128.NENs", "14.NENs", "141.NENs", "1417", "1421",
"1434", "144.NENs", "146.NENs", "147.NENs", "148.NENs",
"15.NENs", "208", "215", "22.NENs", "25.NENs", "30.NENs",
"33.NENs", "3404", "37.NENs", "44.NENs", "45.NENs",
"506", "511", "53.NENs", "55.NENs", "64.NENs", "709",
"75.NENs", "78.NENs", "827", "837", "85.NENs", "88.NENs",
"92.NENs", "95.NENs", "BP42", "BP59", "GEP35","1601","1501","1409",
"1203", "130.NENs","501","62.NENs", "830", "1416", "901",
"61.NENs", "3403", "1413", "GEP26", "GEP29", "GEP27",
"120.NENs", "54.NENs", "1008", "132.NENs", "1406",
"BP38", "56.NENs", "503", "1432", "101", "165.NENs",
"911", "1317", "1015", "822", "924", "1311", "228",
"1412", "81.NENs","905", "135.NENs",
"401", "52.NENs", "101.NENs", "1016", "106.NENs",
"127.NENs", "129.NENs", "1316", "1411",
"145.NENs", "159.NENs", "216", "223", "27.NENs","31.NENs",
"32.NENs", "34.NENs", "3402", "3407", "38.NENs", "49.NENs",
"59.NENs", "72.NENs", "917", "97.NENs", "GEP25", "GEP34",
"BP56" )

df$Tumor_Sample_Barcode[!df$Tumor_Sample_Barcode %in% sampleOrder]

pdf(paste0(outdir, "/gistic.heatmap.pdf"),
height = 6, width = 6)
gisticOncoPlot(gistic = laml.gistic,
clinicalData = df,
clinicalFeatures = c("Subtipo","Tumor.Primario"),
sortByAnnotation = T,
top = nrow(cnvdis),
annotationColor = list("Subtipo" = c("SN1" = "gold",
    "SN2" = "cornflowerblue",
    "SN3" = "tomato1"), "Tumor.Primario" = c(
        "Colorrectal" = "orange",
        "Pulmón" = "blue",
        "I.Delgado" = "yellow",
        "Páncreas" = "green",
        "Estómago" = "pink"
    )),
removeNonAltered = FALSE,
sampleOrder = sampleOrder)
dev.off()







