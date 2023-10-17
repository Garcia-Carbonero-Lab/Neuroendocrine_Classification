#For this script we need  geneset.yml
library(GSEABase)
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]

source("functions/genesets/genesets.functions.R")

#We create output folder
if (dir.exists(paste0(wkdir, "/genesets")) == F) {

    dir.create(paste0(wkdir, "/genesets"))

}




hallmarks <- getGmt(paste0(datadir, "/genesets/h.all.v2023.1.Hs.symbols.gmt"))
cells <- getGmt(paste0(datadir, "/genesets/c8.all.v2023.1.Hs.symbols.gmt"))


cells_filter <- c("BUSSLINGER",
"GAO",
"MURARO",
"TRAVAGLINI",
"VANGURP")

cells <- cells[grep(paste0(cells_filter, collapse = "|"), names(cells))]

gobp <- getGmt(paste0(datadir, "/genesets/c5.go.bp.v2023.1.Hs.symbols.gmt"))

gobp_filter <- c("HORMONE",
"INSULINE",
"ENDOCRINE",
"EPITHELIAL_CELL_PROLIFERATION",
"WNT",
"NOTCH")

gobp <- gobp[grep(paste0(gobp_filter, collapse = "|"), names(gobp))]


biocarta <- getGmt(paste0(datadir, "/genesets/c2.cp.biocarta.v2023.1.Hs.symbols.gmt"))

kegg <- getGmt(paste0(datadir, "/genesets/c2.cp.kegg.v2023.1.Hs.symbols.gmt"))


eliminate_kegg <- c("KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",
"KEGG_PATHOGENIC_ESCHERICHIA_COLI_INFECTION", 
"KEGG_PROXIMAL_TUBULE_BICARBONATE_RECLAMATION",
"KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS", "KEGG_PRIMARY_IMMUNODEFICIENCY",
"KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM",
"KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC", 
"KEGG_DILATED_CARDIOMYOPATHY",
"KEGG_PARKINSONS_DISEASE",
"KEGG_ALZHEIMERS_DISEASE",
"KEGG_LEISHMANIA_INFECTION",
"KEGG_ALLOGRAFT_REJECTION",
"KEGG_GRAFT_VERSUS_HOST_DISEASE",
"KEGG_VIRAL_MYOCARDITIS", "KEGG_AUTOIMMUNE_THYROID_DISEASE",
"KEGG_HUNTINGTONS_DISEASE",
"KEGG_OLFACTORY_TRANSDUCTION",
"KEGG_CARDIAC_MUSCLE_CONTRACTION")


kegg <- kegg[grep(paste0(eliminate_kegg, collapse = "|"), names(kegg), invert = T)]


pid <- getGmt(paste0(datadir, "/genesets/c2.cp.pid.v2023.1.Hs.symbols.gmt"))


reactome <- getGmt(paste0(datadir, "/genesets/c2.cp.reactome.v2023.1.Hs.symbols.gmt"))


filter_reactome <-  c("REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
 "REACTOME_TRANSLATION",
"REACTOME_RRNA_PROCESSING", "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
"REACTOME_SWITCHING_OF_ORIGINS_TO_A_POST_REPLICATIVE_STATE",
"REACTOME_ORC1_REMOVAL_FROM_CHROMATIN", "REACTOME_MITOCHONDRIAL_TRANSLATION",
"REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING",
"REACTOME_TP53_REGULATES_METABOLIC_GENES", "REACTOME_DNA_REPLICATION",
"REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION",
"REACTOME_REGULATION_OF_PTEN_STABILITY_AND_ACTIVITY",
"REACTOME_BIOLOGICAL_OXIDATIONS", "REACTOME_NEUTROPHIL_DEGRANULATION",
"REACTOME_SIGNALING_BY_INTERLEUKINS", "REACTOME_NEURONAL_SYSTEM")

reactome <- reactome[grep(paste0(filter_reactome, collapse = "|"), names(reactome))]


toGmt(hallmarks,paste0(wkdir, "/genesets/hallmarks.gmt"))
toGmt(cells,paste0(wkdir, "/genesets/cells.gmt"))
toGmt(gobp,paste0(wkdir, "/genesets/gobp.gmt"))
toGmt(biocarta,paste0(wkdir, "/genesets/biocarta.gmt"))
toGmt(kegg,paste0(wkdir, "/genesets/kegg.gmt"))
toGmt(pid,paste0(wkdir, "/genesets/pid.gmt"))
toGmt(reactome,paste0(wkdir, "/genesets/reactome.gmt"))


