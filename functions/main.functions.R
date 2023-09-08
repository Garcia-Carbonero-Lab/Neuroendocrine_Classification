#function:zscore.rows2
#input:
#      x: matrix with features as rows and samples as columns

#output: scaled matrix

#description: matrix scaled removing mean and dive by standard deviation



zscore.rows2 <- function(x){
  return(t(apply(x, 1, function(x) (x - mean(na.omit(x)))/sd(na.omit(x)))))
}



#function:create.met.anot
#input:
#      epic:EPIC arrays anotation 

#output: data.frame with probe names as row.names and columns with
# information about chromosome, start and end position, Island,
# fgroup, DMR genes and UCSC reference group
#description: takes Epic anotation and transform it to easy anotation matrix for
#epic probes


create.met.anot <- function(epic) {

genes <-   strsplit(as.character(epic$UCSC_RefGene_Name),';')
genes[lengths(genes) == 0] <- NA
group <- strsplit(as.character(epic$UCSC_RefGene_Group),';')
group[lengths(group) == 0] <- NA

met_anot <- data.frame('ID' = rep(epic$IlmnID, sapply(genes, length))
    , "genes" = unlist(genes)
    , "CpG_chrm" = rep(epic$CHR_hg38, sapply(genes, length))
    , "CpG_beg" = rep(epic$Start_hg38, sapply(genes, length))
    , "CpG_end" = rep(epic$End_hg38, sapply(genes, length))
    , "UCSC_RefGene_Group" = unlist(group)
    , "ISLAND" = rep(epic$Relation_to_UCSC_CpG_Island, sapply(genes, length))
    , "FGROUP" = rep(epic$Regulatory_Feature_Group, sapply(genes, length))
    , "DMR" = rep(epic$DMR, sapply(genes, length)))

met_anot <- met_anot[!duplicated(paste(met_anot$ID,met_anot$genes)),]

    met_anot <- met_anot %>% group_by(ID,CpG_chrm,CpG_beg,
    CpG_end,ISLAND,FGROUP, DMR) %>% 
        dplyr::summarise(genes = toString(genes),
        UCSC_RefGene_Group = toString(UCSC_RefGene_Group))
    
    met_anot$genes <- gsub(', ','_',met_anot$genes)
        
    
    rownames(met_anot) <- met_anot$ID
    return(met_anot)
}


theme_general <- function(plot) {
    themed <- plot +
    theme(text=element_text(size=40),
          axis.text.x = element_text(angle = 0 ,vjust = 1, hjust = 0.5),
          legend.key.height = unit(2,'cm'),
          aspect.ratio=1/1
          )
}
