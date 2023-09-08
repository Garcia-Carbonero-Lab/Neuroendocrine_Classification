
# For use this script we need the environment preprocess
library(biomaRt)


#load functions
source("functions/MOFA/MOFA.functions.R")
source("functions/MOFA/adictional.functions.R")
source("functions/main.functions.R")

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

if (dir.exists(paste0(wkdir, "/CNV/preprocess")) == F) {

    dir.create(paste0(wkdir, "/CNV/preprocess"))

}


outdir <- paste0(wkdir, "/CNV/preprocess")

#Create anotation table peaks

amp <- read.table(paste0(wkdir, "/CNV/GISTIC/amp_genes.conf_99.txt")
                  , sep = '\t'
                  ,header = T
                  , check.names = F)


del <- read.table(paste0(wkdir, "/CNV/GISTIC/del_genes.conf_99.txt")
                  , sep = '\t'
                  ,header = T
                  , check.names = F)

cnv_peak <- read.table(paste0(wkdir, "/CNV/GISTIC/all_lesions.conf_99.txt")
                  , sep = '\t'
                  ,header = T
                  ,row.names = 1
                  , check.names = F
                  , )

cnv_peak$Descriptor <- gsub(" *","",cnv_peak$Descriptor)

anot <- cnv_peak[,1:8]
anot <- anot[grep("CN values", invert = T, rownames(anot)),]

ensembl <- useEnsembl(biomart = 'genes', 
                       dataset = 'hsapiens_gene_ensembl',
                       version =  "GRCh37")


ids <- gsub("\\(.*", "", anot[,"Region Limits"])
ids <- gsub("chr", "",ids)


genes <-  lapply(1:length(ids),function(x){

      biom <- getBM(attributes = c('chromosome_name', 'start_position', "end_position", "band", "hgnc_symbol"),
      filters = c('chromosome_name', 'start', "end"),
      values = list(ids[[x]][[1]], ids[[x]][[2]],ids[[x]][[3]]), 
      mart = ensembl)

      biom$id <- rownames(anot)[x]
      return(biom)

}
)

genes <- do.call(rbind,genes)
genes <- genes[genes$hgnc_symbol != "",]

anot2 <- lapply(1:length(ids), function(x){

    df <- genes[genes$id == rownames(anot)[x],]
    genes <-  paste0(df[,"hgnc_symbol"], collapse = ";")

    return(data.frame(name = rownames(anot)[x],
    genes = genes,
    ngenes = nrow(df)
     ))
})

anot2 <- do.call(rbind,anot2)



anot <- merge(anot, anot2, by.x = "row.names", by.y = "name" ,all.x = T)

write.table(anot,
paste0(outdir, "/CNV.peaks.anotation.txt"),
sep = "\t",
row.names = F,
col.names = T)

################################################################################
#preprocess peaks
# We separate discrete and continous
cnvcont <- cnv_peak[grep('CN values',rownames(cnv_peak)),]
cnvdis <- cnv_peak[grep('CN values',rownames(cnv_peak), invert = T),]

#Change row.names

rownames(cnvcont) <- cnvcont$Descriptor
rownames(cnvdis) <- cnvdis$Descriptor


cnvcont <- unique(cnvcont[,9:ncol(cnvcont)])
cnvcont <- cnvcont[,colnames(cnvcont) != '']
cnvcont <- cnvcont[,complete.cases(cnvcont[1,])]

cnvdis <- unique(cnvdis[,9:ncol(cnvdis)])
cnvdis <- cnvdis[,colnames(cnvdis) != '']
cnvdis <- cnvdis[,complete.cases(cnvdis[1,])]

cnvdis <- apply(cnvdis, c(1,2) , as.character)
cnvdis <- apply(cnvdis, c(1,2) , as.numeric)

colnames(cnvcont) <- gsub('\\.',' ', colnames(cnvcont))
colnames(cnvdis) <- gsub('\\.',' ', colnames(cnvdis))


cnvdis[1,] <- ifelse(cnvdis[1,] ==1, "Amplification", "Wild.type")

cnvdis[cnvdis == "1"] <- "Deletion"
cnvdis[cnvdis == "0"] <- "Wild.type"

write.table(cnvcont,
paste0(outdir, "/CNV.peaks.con.txt"),
sep = "\t",
row.names = T,
col.names = NA)


write.table(cnvdis,
paste0(outdir, "/CNV.peaks.dis.txt"),
sep = "\t",
row.names = T,
col.names = NA)
