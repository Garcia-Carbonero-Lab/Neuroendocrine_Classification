



#function:write.sampleset
#input:
#      sets: list of the samples sets
#      out.file Name of output file
#output: gmt file
#description: generate gmt from list of geneset 


write.sampleset <- function (sets,out.file){
  file.handle <- file(out.file, "wt")
  
  for (n in 1:length(sets)){
    write(c(names(sets)[n],'.',sets[[n]]),
          file = file.handle, append = T, sep = "\t", ncolumns = 5000)
  }
  close(file.handle)
}

#function: anotate.gmt
#input:
#      genesets: list of the gene sets
#      expression: expression matrix with genes as rows and samples as columns
#output: list of genesets
#description: In case of synonyms between geneset and expression matrix, 
#             remplace gene symbol in genset by gene symbol in expression matrix

anotate.gmt <- function(genesets, expression){

#convert the geneset into list
genesets <- geneIds(genesets)


genes_geneset <- unique(unlist(lapply(genesets,function(x){return(x)})))

#create a dataframe with the changes that we neeed to do
genes_expression <- rownames(expression)[!rownames(expression) %in% genes_geneset]

genes <- genes_geneset[!genes_geneset %in% rownames(expression)]
symbol_genes <- alias2SymbolTable(genes)


df_geneset <- data.frame("Geneset" = genes, "Symbol" = symbol_genes)

df_expression <- data.frame("Expression" = genes_expression,
"Symbol" = alias2SymbolTable(genes_expression))

df <- merge(df_geneset, df_expression, by="Symbol")
geneset_complete <- data.frame("Geneset" = genes_geneset[!genes_geneset %in% df$Geneset])

df <- merge(df,geneset_complete, by = "Geneset", all = T)
df <- df[!duplicated(df$Geneset),]
rownames(df) <- df[,1]

#apply changes in gmt
genesets_new <- lapply(genesets, function(x){
dft <- df[x,]
x[complete.cases(dft$Expression)] <- dft$Expression[complete.cases(dft$Expression)]
return(x)
})

names(genesets_new) <- names(genesets)
return(genesets)
}