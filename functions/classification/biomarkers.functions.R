
#function:deg.function
#input:
#      exp: matrix with features as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      cov: covariables to adjust the differentaill expression
#      group: name of the colum to generate the differential expression
#      color: vector with color of levels in group
#      control: indicate the label of the control in differentially expression
#      outdir: Path to write the plots
#      flag: string to add in names of the files


#output: result of the differentially expression analysis with limma

#description: function to generate differential expression analysis



deg.function <- function(exp,
base, group,
cov = NULL,
outdir,
control = NULL,
flag) {
base[,group] <- as.factor(as.character(base[,group]))
base <- base[colnames(exp),]


if (is.null(control)){

    control <- levels(base[,group])[1]
}


if (is.null(cov)){
  design <- data.frame(Group = base[,group])
  design$Group <- relevel(design$Group, control)
  design <- model.matrix( ~ 0 + Group, data = as.data.frame(design))
  }else{
    
    design <- data.frame(Group = base[,group])
    
    for (f in cov){design[,f] <- base[,f]}
    design$Group <- relevel(design$Group, control)
    design <- model.matrix(as.formula(paste0("~ 0 + Group + ",
    paste0(cov, collapse = " + "))), data = as.data.frame(design))
  }
  
  CN <- combn(colnames(design)[1:length(levels(base[,group]))],2,
  simplify = FALSE)
  
  CN <- unlist(lapply(CN,function(x){paste0(x, collapse = "-")}))
  
  
  rownames(design) <- rownames(base)
  fit<-lmFit(exp, design)
  
  

  contrast.matrix <- makeContrasts(contrasts= CN,
  levels= design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

res <- list()
for (comb in CN){

result <- topTable(fit2, coef=comb, adjust="BH",number = nrow(exp))
write.table(result,
paste0(outdir,"/TopTable_",flag, "_",comb ,".txt"),
sep = "\t",
row.names = T,
col.names = NA)

res[[comb]] <- result
}

res <- do.call(rbind, res)
return(res)
}


#function:deg.plot
#input:
#      data: matrix with features as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      group.x: name of the colum for x axis in boxplot
#      group.fill: Column name to be the fill in boxplot
#      color.fill: vector with color of levels in cgroup.fill
#      cov: covariates to adjust the limma model
#      outdir: Path to write the plots
#      flag: string to add in names of the files
#      genes: genes to show in boxplots
#      filter.table: filter table with genes in genes

#output: boxplots with limma p value and a table with ther esults of limma 

#description: function to obtain boxplots showing the results of limma


deg.plot <- function(data,
base,
group.x,
group.fill,
color.fill,
cov = NULL,
flag,
outdir,
genes,
filter.table = T) {


# eliminate problematic symbols in rownames
genes <- gsub("-","_",genes)
colnames(data) <- gsub("-", "_", colnames(data))
genes <- gsub("\\.","_",genes)
colnames(data) <- gsub("\\.", "_", colnames(data))


# select genes that are avaiable in data
genes <- genes[genes %in% colnames(data)]




if (dir.exists(paste0(outdir, "/", flag)) == F) {

    dir.create(paste0(outdir, "/", flag))

}



data <- data[rownames(base),,drop = F]

df<- cbind(data,base[,c(group.x,group.fill)])
df[,group.x] <- as.factor(as.character(df[,group.x]))
df[,group.fill] <- as.factor(as.character(df[,group.fill]))


if (group.x != group.fill){

base[,group.x] <- as.factor(as.character(base[,group.x]))
base[,group.fill] <- as.factor(as.character(base[,group.fill]))

degs <- lapply(levels(base[,group.x]), function(x){

dfp <- base[base[,group.x] == x,]

# differentially expression analysis
degp <- deg.function(as.data.frame(t(data))[,rownames(dfp)],
dfp,
group = group.fill,
cov = cov,
outdir = paste0(outdir, "/", flag),
flag = x)

#Prepare limma reults to use it gor ggpurb
degp$Row.names <- rownames(degp)
degp <- as_tibble(degp)
fdegp <- degp %>%
separate(Row.names, c("Group", "gene"), sep = "\\.") %>%
separate(Group, c("group1", "group2"), sep = "-") %>%
mutate(
  group1 = str_replace(group1, "Group", ""),
  group2 = str_replace(group2, "Group", "")

) %>%
rename(
  p = P.Value,
  p.adj = adj.P.Val 
) 

# filter genes in genes
if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}
# generate p.adj.signif column with stars
fdegp[,"p.adj.signif"] <- rep("ns", nrow(fdegp))

if (nrow(fdegp[fdegp$p.adj < 0.001,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.001,"p.adj.signif"] <- "***"
}

if (nrow(fdegp[fdegp$p.adj < 0.01,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.01 & fdegp[, "p.adj"] > 0.001,"p.adj.signif"] <- "**"
}

if (nrow(fdegp[fdegp$p.adj < 0.05,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.05 & fdegp[, "p.adj"] > 0.01,"p.adj.signif"] <- "*"
}

fdegp[,group.x] <- as.factor(rep(x,nrow(fdegp)))
return(fdegp)

})

degs <- degs %>%
bind_rows(degs) %>%
distinct

# write results
write.table(degs,
paste0(outdir, "/",flag, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)


}else{

base[,group.fill] <- as.factor(as.character(base[,group.fill]))

# differentially expression analysis

degp <- deg.function(as.data.frame(t(data)),
base,
group = group.fill,
cov = cov,
outdir = paste0(outdir, "/", flag),
flag = "Complete")
degp$Row.names <- rownames(degp)

degp <- as_tibble(degp)
#Prepare limma reults to use it gor ggpurb

fdegp <- degp %>%
separate(Row.names, c("Group", "gene"), sep = "\\.") %>%
separate(Group, c("group1", "group2"), sep = "-") %>%
mutate(
  group1 = str_replace(group1, "Group", ""),
  group2 = str_replace(group2, "Group", "")

) %>%
rename(
  p = P.Value,
  p.adj = adj.P.Val 
)

if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}

# generate p.adj.signif column with stars

fdegp[,"p.adj.signif"] <- rep("ns", nrow(fdegp))

if (nrow(fdegp[fdegp$p.adj < 0.001,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.001,"p.adj.signif"] <- "***"
}

if (nrow(fdegp[fdegp$p.adj < 0.01,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.01 & fdegp[, "p.adj"] > 0.001,"p.adj.signif"] <- "**"
}

if (nrow(fdegp[fdegp$p.adj < 0.05,]) > 0){
fdegp[fdegp[, "p.adj"] < 0.05 & fdegp[, "p.adj"] > 0.01,"p.adj.signif"] <- "*"
}

degs <- fdegp
write.table(degs,
paste0(outdir, "/", flag, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)

}


for (row in genes) {

#Prepare data for ggpurb
if (group.x != group.fill){

dfa <- df[,c(row,group.x,group.fill)]

stat.test <- dfa %>%
  group_by_at(vars(group.x)) %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.deg <- degs %>%
filter( gene == row)

stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)



stat.test <- stat.test %>% 
unite(col= "order", c( as.name(group.x), group1:group2), sep = "_",
remove = FALSE) 

stat.deg <- stat.deg %>% 
unite(col= "order", c(
group1:group2),
sep = "_",
remove = FALSE)%>%
distinct()


}else{

dfa <- df[,c(row,group.fill)]
 
stat.test <- dfa %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


stat.deg <- degs %>%
filter( gene == row)


stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)



stat.test <- stat.test %>% 
unite(col= "order", c(group1:group2), sep = "_",
remove = FALSE) 



stat.deg <- stat.deg %>% 
unite(col= "order", c(
group1:group2),
sep = "_",
remove = FALSE)%>%
distinct()



}

# add the pvalue of limma to rstatix results
stat.test <- stat.test %>%
mutate(
  p.adj.signif = stat.deg$p.adj.signif
)

stat.test[[group.fill]] <- NA

# create boxplots
bxp <- ggplot(
  df, aes(x = .data[[group.x]],
  y = .data[[row]], fill= .data[[group.fill]])) +
  geom_boxplot(width = 0.6, outlier.shape = NA,
  aes(fill= .data[[group.fill]])) +

  scale_fill_manual(values = as.list(c(color.fill,"blue", "green"))) +

  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7)) +
  
  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = T
) +


  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme_classic() +

theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
        ) +

        labs(x = "", y = row)


bxp$layers[[3]]$aes_params$size <- 1.2
bxp$layers[[3]]$aes_params$label.size <- 12 
bxp

ggsave(
        plot = bxp, filename = paste0(
            outdir, "/", flag, "/Boxplot.subtype_", row, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = 16,
        height = 12
    )

}


}