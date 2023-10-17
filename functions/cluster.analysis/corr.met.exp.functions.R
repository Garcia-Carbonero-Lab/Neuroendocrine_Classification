#function: corr.function
#input:
#      exp: expression matrix with genes as rows and samples as columns
#      met: methylation matrix with probes as rows and samples as columns
#      id: id of methylation
#      gene : gene symbol

#output: matrix witjh the correlation coeficent and r2 and p value of the 
#correlation between gene and the probe 

#description: function to correlate gene expression and methylation probes


corr.function <- function(exp, met, id, gene){

df <- cbind(
    as.numeric(exp[gene,]),
    as.numeric(met[id,]))

df <- as.data.frame(df)

colnames(df) <- c("exp", "met")


cor <- cor.test(as.numeric(df[,"exp"]),as.numeric(df[,"met"]))
lmx <- summary(lm(exp ~ met, df))

r2 <- lmx$adj.r.squared[[1]]
p.val <- cor$p.value[[1]]
est <- cor$estimate[[1]]
return(list("gene" = gene,
"CpG" = id,
"corr" = est,
"r2" = r2,
"p.val" = p.val))
}

#function: corr.met.exp.unique
#input:
#      exp: expression matrix with genes as rows and samples as columns
#      met: methylation matrix with probes as rows and samples as columns
#      met_anot anotation matrix of methlatio probes

#output: matrix with the correlation coeficent and r2 and p value of the 
#correlation between gene and the probe 

#description: function to correlate gene expression and methylation probes

corr.met.exp.unique <- function(exp, met, met_anot){

met_anot <- met_anot %>%
separate_rows(genes,sep = "_")

#filter data by overlaping genes

#transform genes in symbol names
symbol_genes <- alias2SymbolTable(rownames(exp))

#Eliminate duplicated genes
anot <- data.frame(

    "symbol" = symbol_genes,
    "alias" =  rownames(exp),
    "variance" = apply(exp,1,var)

)

exp <- exp[order(anot$variance),]
anot <- anot[order(anot$variance),]

exp <- exp[!duplicated(anot$symbol),]
anot <- anot[!duplicated(anot$symbol),]

rownames(exp)[is.na(anot$symbol) == F] <- anot$symbol[is.na(anot$symbol) == F]

symbol_met <- alias2SymbolTable(met_anot$genes)
met_anot$genes[is.na(symbol_met) == F] <- symbol_met[is.na(symbol_met) == F]

met_anot <- met_anot[met_anot$genes %in% rownames(exp),]
met_anot <- met_anot[met_anot$ID %in% rownames(met),]

met <- met[met_anot$ID,]
##########################################################
exp <- exp[,colnames(met)]


corr.df <- lapply(met_anot$ID,function(x){
gene <- met_anot$genes[met_anot$ID == x]
id <- x

result <- lapply(gene,function(y){
    return(corr.function(exp, met, id, y))
})

result <- as.data.frame(do.call(rbind,result))

#eliminate list


return(result)
})

corr.df <- do.call(rbind, corr.df)

for (col in 1:ncol(corr.df)){
    corr.df[,col] <- unlist(corr.df[,col])
}
return(corr.df)
}


#function: correlation.plot
#input:
#      exp: expression matrix with genes as rows and samples as columns
#      met: methylation matrix with probes as rows and samples as columns
#      anot: anotation matrix of methlation probes and gene asociation and correlation
#      base: data.frame with samples as rows and features as columns
#      genes: genes to plot
#      outdir: path to write results
#output: matrix with the correlation coeficent and r2 and p value of the 
#correlation between gene and the probe 

#description: function to correlate gene expression and methylation probes

correlation.plot <- function(exp,
met,
anot,
base,
genes,
outdir,
group,
color.group,
flag) {


#Eliminate - in gene names

genes <- gsub("-","_", genes)
rownames(exp) <- gsub("-","_", rownames(exp))
anot$genes <- gsub("-","_", anot$genes )


exp <- exp[,rownames(base)]
met <- met[,rownames(base)]



if (dir.exists(paste0(outdir, "/corr.plots")) == F) {

    dir.create(paste0(outdir, "/corr.plots"))

}

for (gene in genes){
print(gene)
probes <- anot$ID[anot$genes == gene]
probes <- probes[probes %in% rownames(met)]
for (probe in probes){



df <- data.frame(
    group = base[,group],
    "exp" = as.numeric(exp[gene,]),
    "met" = as.numeric(met[probe,]),
    row.names = rownames(base)
)





gene_met <- paste(gene, probe, sep = "_")

colnames(df) <- c(group,
gene,
gene_met)


df[,group] <- as.factor(df[,group])


plt1 <- ggplot(df, aes_string(x = gene,
y = gene_met,
colour = group,
group = 1)) +
    theme_classic() +
    geom_point(size = 3) + 
scale_color_manual(values= color.group) +

   
    stat_cor(size = 10) +

    stat_smooth(method = "lm", color = "black") +
    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
    )

    ggsave(
        plot = plt1, filename = paste0(
            outdir,"/corr.plots/",gene, "_", probe,"_", flag,".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 17,
        height = 12
    )


    plt2 <- ggplot(df, aes_string(x = gene,
y = gene_met)) +
    theme_classic() + stat_cor(size = 10) +
stat_smooth(method = "lm")  +
geom_point(size = 3) +
    facet_wrap(as.formula(paste0("~ ", group)),
            drop = T, scales = "free_x", nrow =1
        ) +
    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
    ) 

ggsave(
        plot = plt2, filename = paste0(
            outdir,"/corr.plots/", gene, "_", probe,"_", flag, ".wrapped.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 25,
        height = 12
    )
}
}
}