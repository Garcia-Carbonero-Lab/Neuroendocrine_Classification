#Functions adictional for MOFA analisis

#function:obtain_Z
#input:
#      mod: MOFA model
#       
#output: matrix with samples in rows and factors in columns

#description: Function to obtain the value of Factors in each sample



obtain_Z <- function(mod){
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]
    return(Z)
}



#Figura2
############################################################################

#function:kruskal
#input:
#      data: data.frame with information for columns
#      Z: matrix with samples in rows and factors in columns
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      columns: colnames of variables to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      
#output: files and plot results of the Kruskall wallis  analyze the associations between discrete variables and model
# factors

#description: Function to applay Kruskall Wallis test between Factors and discrete variables.
# 

kruskal <- function(data, Z, columns, flag, outdir, p.th) {
    # We calculate the kruskal wallis 
    ano_df <- lapply(columns, function(x) {
    vf <- lapply(colnames(Z), function(y) {
            k <- kruskal.test(as.formula(paste(y, " ~ ", x)),
                data = data)
            f <- k[[1]]
            p <- k[[3]]

            return(c("Chi-squared" = f, "p.value" = p))
        })

        names(vf) <- colnames(Z)
        vf <- as.data.frame(do.call(rbind, vf))
        colnames(vf) <- paste(x, colnames(vf))
        return(as.data.frame(vf))
    })

    ano_df <- as.data.frame(t(do.call(cbind, ano_df)))
    

    # we obtain only the p.value
    p.df <- ano_df[grep("p.value", rownames(ano_df)), ]

    for (col in colnames(p.df)) {

    p.df[, col] <- p.adjust(p.df[, col], method = "BH")

    }

    p.adj.df <- p.df
    rownames(p.adj.df) <- gsub("p.value","p.adj",rownames(p.adj.df))
    ano_df <- rbind(ano_df,p.adj.df)
    ano_df <- t(ano_df[order(rownames(ano_df)),])

    if (dir.exists(paste0(outdir, "/discrete.variables/")) == F) {

    dir.create(paste0(outdir, "/discrete.variables/"))

    }

    if (dir.exists(paste0(outdir, "/discrete.variables/", flag)) == F) {

    dir.create(paste0(outdir, "/discrete.variables/", flag))

    }

    write.table(ano_df, paste0(
        outdir, "/", flag, "_Kruskal.txt" 
    ),
    sep = "\t",
    row.names = T,
    col.names = NA
    )

    # We will show only color result in those that are significant
    # We only show variables with any significance result

    p.df[p.df > p.th] <- 1
    p.df <- p.df[apply(p.df, 1, sum) != nrow(p.df), ]
    rownames(p.df) <- gsub(" p.value", "", rownames(p.df))


    return(p.df)
}



#function: heatmap.kruskal
#input:
#      mod: MOFA model
#      columns: colnames of variables to plot.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      
#output:heatmap indicating the p-value of the features associated with factors

#description: Function to generate the heatmap with the results of Kruskal wallis
# between features adn factors of the model
# 

heatmap.kruskal <- function(mod,columns,
flag, outdir,
p.th = 0.05){
Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    df <- mod@samples_metadata

    p.df <- kruskal(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )
    
    col <- colorRamp2(
                c(max(-log10(t(p.df))), 0),
                c("red", "grey"))

    
    #Change names order
    hm1 <- Heatmap(as.matrix(-log10(t(p.df))),
            show_row_names = T,
            show_column_names = T,
            col = col,
            cluster_rows = F,
            cluster_columns = T,
            show_row_dend = F,
            show_column_dend = F,
            heatmap_legend_param = list(
                title = "-log10(q-valor)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm")
            ),
            row_names_side = "right",
            row_names_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            column_names_rot = 45,
            rect_gp = gpar(colo = "black", lwd =2)
        )




p <- draw(hm1, heatmap_legend_side = "bottom",
        padding = unit(c(2, 10, 2, 20), "mm"))


   pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = 15, height= 20
        )
        print(p)
        dev.off()


return(p)

} 

#function:surv.model
#input:
#      mod: MOFA model
#      time: name of the column with time information
#      event: name of the column with event information.
#      flag: name to add to the name of output files.
#      outdir: Path to write outputs
#output: Cox model using factor as variables and forest plot

#description: Function to create a Cox model using factors as variables to check
#the association between factors and survival

surv.model <- function(mod, time, event, flag) {
  
    SurvObject <- Surv(
        time = mod@samples_metadata[, time],
        event = mod@samples_metadata[, event]
    )

    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    fit <- coxph(SurvObject ~ Z)

    s <- summary(fit)
    coef <- s[["coefficients"]]

    sur.df <- data.frame(
        factor = factor(rownames(coef), levels = rev(rownames(coef))),
        p      = coef[, "Pr(>|z|)"],
        coef   = coef[, "exp(coef)"],
        lower  = s[["conf.int"]][, "lower .95"],
        higher = s[["conf.int"]][, "upper .95"]
    )

    write.table(sur.df,
        sep = "\t",
        paste0(outdir, "/", flag, "_Hazard_ratio_factors.txt"),
        row.names = F,
        col.names = T
    )

    sur.df$factor <- gsub("Z", "", sur.df$factor) 
    
    sur.df <- sur.df[order(sur.df$factor, decreasing = T ),]


    pt <- ggplot(sur.df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
        geom_pointrange(col = "black", size = 1.8, linewidth = 1.8) +
        coord_flip() +
        scale_x_discrete(limits= rev) +
    labs(y = "Hazard Ratio", x = "") +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        theme_bw() + 
        theme(text = element_text(size = 35, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(colour = "black"))

    ggsave(pt,
        filename = paste0(outdir, "/", flag,
        "_Hazard_ratio_factors.pdf"),
        device = "pdf",
        dpi = 900, width = 10, height = 10
    )

    return(pt)
}





#function:immune_lm
#input:
#      data: matrix with samples in rows and immune cells as columns
#      mod: MOFA model
#      outdir: Path to write outputs
#      flag: Name to add in files

#output: Results of lineal model that associates factors with immune infiltration

#description: Function to associate factors with immune cell infiltration by lineal models


immune_lm <- function(data, mod, flag, outdir){
Z <- obtain_Z(mod)
df <- merge(Z, data, by = "row.names")
p.val <- list()
r2 <- list()
estimate <- list()
for (f in colnames(Z)){
    lmx <- lm(as.formula(paste0(f, " ~ ",
    paste0(colnames(data), collapse = " + "))),
    data = df)
    sum <- summary(lmx)
    
    p.val[[f]] <- sum$coefficients[,"Pr(>|t|)"]
    r2[[f]] <- sum$r.squared
    estimate[[f]] <- sum$coefficients[,"Estimate"]
}

p.val <- do.call(rbind, p.val)

r2 <- do.call(rbind, r2)

estimate <- do.call(rbind, estimate)

colnames(r2) <- "R2"

write.table(p.val,
paste0(outdir,  "/",
            flag, "_Linear_model_p.value.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(r2,
paste0(outdir,  "/",
            flag,"_Linear_model_r2.txt"),
sep = "\t",
row.names = T,
col.names = NA)
write.table(estimate,
paste0(outdir,  "/",
            flag,"_Linear_model_estimate.txt"),
sep = "\t",
row.names = T,
col.names = NA)

estimate[p.val > 0.05] <- 0
estimate <- estimate[,-1]



col <- list("R2" = colorRamp2(c(1,0),
            c("orange", "white")
        ))


anot_legend_param <-  list(
        "R2" = list(title_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        labels_gp = gpar(fontsize = 60),
        grid_height = unit(2, "cm"),
        legend_width = unit(20, "cm"),
        direction = "horizontal",
        title_position = "topcenter",
        border = "black"
        ))


        ha_column <- HeatmapAnnotation(
            df = as.data.frame(r2),
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = 60),
            border = unit(40, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 60),
            annotation_legend_param = anot_legend_param)


hm <- Heatmap(as.matrix(t(estimate)),
            top_annotation = ha_column,
            show_row_names = T,
            show_column_names = T,
            col = colorRamp2(
                c(max(estimate), 0, min(estimate)),
                c("red", "white", "blue")
            ),
            name = "Heatmap",
            gap = unit(5, "mm"),
            cluster_rows = F,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "spearman",
            heatmap_legend_param = list(
                title = "Coeficiente",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            rect_gp = gpar(col = "black", lwd = 2),
            column_names_rot = 45
        )

        p <- draw(hm, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        padding = unit(c(2, 140, 2, 2), "mm"))

        pdf(paste0(
            outdir,  "/",
            flag, "_immune.linear.model.pdf"
        ),
        width = 20, height = 20)
        print(p)
        dev.off()

return(p)
}

#function:gsea.heatmp
#input:
#      expression_gsea: gsea results
#      outdir: Path to write outputs

#output: heatmap indicating the NES of the pathways associates with each factor

#description: Function to obtain gesa heatmap

gsea.heatmp <- function(expression_gsea, outdir){
   
pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
"HALLMARK_MYC_TARGETS_V1",
"KEGG_PROTEASOME", "KEGG_RIBOSOME",
"BIOCARTA_IL2RB_PATHWAY",
"GOBP_HORMONE_METABOLIC_PROCESS",
"HALLMARK_INFLAMMATORY_RESPONSE",
"HALLMARK_INTERFERON_GAMMA_RESPONSE",
"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
"KEGG_APOPTOSIS",
"KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
"KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
"GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
"GOBP_ENDOCRINE_SYSTEM_DEVELOPMENT",
"KEGG_FOCAL_ADHESION",
"HALLMARK_BILE_ACID_METABOLISM",
"KEGG_MAPK_SIGNALING_PATHWAY",
"PID_RB_1PATHWAY",
"PID_FOXM1_PATHWAY",
"REACTOME_NEURONAL_SYSTEM",
"REACTOME_SIGNALING_BY_INTERLEUKINS",
"REACTOME_BIOLOGICAL_OXIDATIONS",
"REACTOME_DNA_REPLICATION")

expression_gsea <- expression_gsea %>%
filter(
    padj < 0.05,
    pathway %in% pathways
) 


write.table(expression_gsea,
sep = "\t",
file= paste0(outdir, "/GSEA.txt"))



expression_gsea <- expression_gsea %>%
select(
pathway, Factor, NES
 ) %>%

pivot_wider(names_from = pathway, values_from = NES) %>%
replace(is.na(.), 0)


expression_gsea <- as.data.frame(expression_gsea)
rownames(expression_gsea) <- expression_gsea$Factor
expression_gsea <- expression_gsea[,-1]

colnames(expression_gsea) <- c("Fosforilación Oxidativa",
"Dianas de MYC",
"Ribosoma",
"Proteosoma",
"Replicación ADN",
"IL2RB",
"Metabolismo Hormonal",
"Respuesta Inflamatoria",
"Interferón Gamma",
"TNFA NFKB",
"Metabolismo de Ácidos Biliares",
"Receptor Células B",
"Apoptosis",
"Receptor Células T",
"Adhesión Focal",
"Señalización FOXM1",
"Oxidaciones Biológicas",
"Sistema Neuronal",
"Señalización por Interleucinas",
"Desarrollo Endocrino",
"WNT Canonica",
"Señalización MAPK"

)

nombres <- c("Hallmarks",
"Hallmarks",
"Kegg",
"Kegg",
"Reactome",
"Biocarta",
"Gobp",
"Hallmarks",
"Hallmarks",
"Hallmarks",
"Hallmarks",
"Kegg",
"Kegg",
"Kegg",
"Kegg",
"Pid",
"Reactome",
"Reactome",
"Reactome",
"Gobp",
"Gobp",
"Kegg")

nombres <- as.factor(nombres)
nombres <- as.data.frame(nombres)
colnames(nombres) <- "base"

col <- list(
    "base" = c("Biocarta" = "steelblue2",
    "Gobp" = "hotpink",
    "Hallmarks" = "gold",
    "Kegg" = "darkolivegreen1",
    "Pid" = "coral2",
    "Reactome" = "purple"))



anot_legend_param <-  list(
        "base" = list(title_gp = gpar(fontsize = 40),
        labels_gp = gpar(fontsize = 40),
        labels_gp = gpar(fontsize = 40),
        grid_height = unit(2, "cm"),
        grid_width = unit(2, "cm"),
        ncol = 1))

ha_row <- rowAnnotation(
            df = nombres,
            show_legend = F,
            show_annotation_name = c(base = F),
            col = col,
            annotation_name_gp = gpar(fontsize = 60),
            border = unit(60, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 60, col = "black"),
            annotation_legend_param = anot_legend_param,
            annotation_name_side = "top")

lgd = Legend(labels = levels(nombres[,1]),
legend_gp = gpar(fill = c("steelblue2",
"hotpink",
"gold",
"darkolivegreen1",
"coral2",
"purple")),
title = "Bases",
title_gp = gpar(fontsize = 60),
labels_gp = gpar(fontsize = 60),
grid_height = unit(2, "cm"),
grid_width = unit(2, "cm"),
ncol = 1
)

hm <- Heatmap(as.matrix(t(expression_gsea)),
            show_row_names = T,
            show_column_names = T,
            left_annotation= ha_row,
            col = colorRamp2(
                c(max(expression_gsea), 0, min(expression_gsea)),
                c("red", "white", "blue")
            ),
            name = "Heatmap",
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "spearman",
            heatmap_legend_param = list(
                title = "NES",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 70),
                labels_gp = gpar(fontsize = 70,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(3.3, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 70),
            column_names_gp = gpar(fontsize = 70),
            row_names_max_width = unit(20,'cm'),
            rect_gp = gpar(col = "black", lwd = 2),
            column_names_rot = 45

        )

        

        p <- draw(hm, heatmap_legend_side = "bottom",
        annotation_legend_side = "right",
        padding = unit(c(2, 120, 2, 2), "mm"),
        annotation_legend_list = lgd)

        pdf(paste0(
            outdir,  "/Heatmap_Figura_2_D.pdf"
        ),
        width = 30, height = 30)
        print(p)
        dev.off()
 
return(hm)
}

#function:complete.boxplot
#input:
#      base: data.frame with samples as rows and features as columns
#      factor: Name of factor to analyse.
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      column: colname of the column to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#               of the labels
#     remove.y: True remove name in y label
#     points: True show each data as a point
#     add.limits: true Add limit in y axis
#     title: true add variable name as title of the plot
#     remove.f: true remove x variable name
#     legend: true show the legend

#output: # Plots of boxplots of all samples that compares
# the factor between the labels of the column feature and the result of the
#Wilcoxon test

#description: Function to create boxplots that compares the levels of the factor
# between the labels of the column feature.


complete.boxplot <- function(base,
factor,
column,
colors,
outdir,
flag,
remove.y = F,
points = F,
add.limits = T,
title = T,
remove.x = F,
legend = F) {

    # we cant't have spaces in levels
    base[, column] <- gsub(" ", ".", base[, column])
    names(colors) <- gsub(" ", ".", names(colors))

    #Regenerate factors
    base[, column] <- factor(base[, column],
            levels = names(colors[[column]])
    )
    # We eliminate missing values

    base[,column] <- as.factor(base[,column])
    base <- base[complete.cases(base[, column]), ]
    base[,column] <- droplevels(base[,column])
    
   
stat.test <- base %>%
  wilcox_test(as.formula(paste0(factor, "~", column))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = column, dodge = 0.8)
    
stat.test[[column]] <- NA
    
    
    plt <- ggplot( base , aes(y = .data[[factor]],
    x = .data[[column]], fill = .data[[column]])) 
        
if (points == T){

        plt <- plt + geom_boxplot(
            width = 0.6, outlier.shape = NA
        ) + 
        geom_point(position = position_jitterdodge()) 
         

    }else{

    plt <-  plt + geom_boxplot(
            width = 0.6
        ) 
    }    

    plt <- plt + scale_fill_manual(values = colors[[column]]) +

        theme_classic() +
        stat_pvalue_manual(stat.test,
        label = "p.adj.signif",  hide.ns = T,
        step.increase = 0.1)  +

        theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(colour = "black"),
        legend.key.size = unit(3,"cm"),
        axis.ticks.length.y = unit(10,"pt"),
    plot.title = element_text(hjust = 0.5, vjust= 4)
) +
        xlab("") 

    if (legend == F){
        
        plt <- plt +  theme(legend.position = "none")

    }


    if (points == T){
    
    plt$layers[[3]]$aes_params$size <- 1
    plt$layers[[3]]$aes_params$label.size <- 20
    plt$layers[[3]]$aes_params$vjust <- 0.5

    }else{

    plt$layers[[2]]$aes_params$size <- 1
    plt$layers[[2]]$aes_params$label.size <- 20
    plt$layers[[2]]$aes_params$vjust <- 0.5
    }
    
    if (add.limits == T){
        limits <- layer_scales(plt)$y$range$range


    add <- limits[[2]] - limits[[1]] /8

    plt <- plt + ylim(limits[[1]],  add)

    }
    
    if (title == T){
        plt <- plt + ggtitle(column) 
    }
    
    if (remove.x == T){
    
    plt <- plt + theme(axis.text.x = element_blank())

    }

    if (remove.y == T){
        plt <- plt + theme(axis.title.y = element_blank())
    }
    
    

    ggsave(
        plot = plt, filename = paste0(
            outdir, "/", flag, "_Complete_heatmap.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )

    return(plt)

}

#function:wrap.boxplot
#input:
#      base: data.frame with samples as rows and features as columns
#      factor: Name of factor to analyse.
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      column: colname of the column to analyze.
#      group: column name of the variable to group the samples
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#               of the labels
#     levels: vector with levels in the order that should appear
#     wrap: True: divide samples in group using wrap function
#     remove.y: True remove name in y label
#     points: True show each data as a point
#     names: True show names in x axis


#output: # Plots of boxplots of all samples that compares
# the factor between the labels of the column feature and the result of the
#Wilcoxon test grouped by a variable

#description: Function to create boxplots that compares the levels of the factor
# between the labels of the column feature grouped by a variable

wrap.boxplot <- function(base,
factor,
column,
group,
colors,
outdir,
flag,
levels,
wrap = T,
remove.y = F,
points = F,
names = F){

    # we cant't have spaces in levels
   base[, column] <- gsub(" ", ".", base[, column])
   base[, group] <- gsub(" ", ".", base[, group])
   base[, group] <- factor(base[,group],
   levels = levels)

       
    base[,column] <- as.factor(base[,column])
    base <- base[complete.cases(base[, column]), ]
    base[,column] <- droplevels(base[,column])


df <- table(base[,group])
rmv <- names(df[df<2])
if (length(rmv > 1)){




base2 <- base[!base[,group] %in% rmv,]

stat.test <- base2 %>%
  group_by_at(vars(group)) %>%
  wilcox_test(as.formula(paste0(factor, "~", column))) %>%
  adjust_pvalue(method = "BH")

new <- stat.test[stat.test[,group] == names(df[df>2])[length(rmv)],]
new[,group] <- rmv
new$p <- rep(1,length(rmv))
new$p.adj <- rep(1,length(rmv))
stat.test <- rbind(stat.test,new)
stat.test <- stat.test %>%
    arrange(as.vector(stat.test[,group])) %>%
    add_significance("p.adj") %>%
    add_xy_position(x = group, dodge = 0.8)


}else{


stat.test <- base %>%
  group_by_at(vars(group)) %>%
  wilcox_test(as.formula(paste0(factor, "~", column))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = group, dodge = 0.8)
  
}
stat.test <- stat.test %>% 
unite(col= "order", c( as.name(group), group1:group2), sep = "_",
remove = FALSE) 

stat.test$names <- paste(stat.test$.y., 
stat.test$group1,
stat.test$group2,
stat.test[, group][[1]],
sep = "_")
    
stat.test[[column]] <- NA


    plt <-
       ggplot( base , aes(y = .data[[factor]],
    x = .data[[group]], fill = .data[[column]])) +
        
        geom_boxplot(
            width = 0.6
        ) +

        scale_fill_manual(values = colors[[column]]) +

        theme_classic() +
        stat_pvalue_manual(stat.test,
        label = "p.adj.signif",  hide.ns = T,
        step.increase = 0.1) 

if (wrap == T){

       plt <- plt + theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        legend.key.size = unit(3,"cm"),
        axis.ticks.length.y = unit(10,"pt")) +
        xlab("")
}else{

      plt <- plt +  theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        strip.text = element_blank(),
        axis.text.x = element_text(colour = "black"),
        legend.key.size = unit(3,"cm"),
        axis.ticks.length.y = unit(10,"pt")) +
        xlab("")

}
   if (remove.y == T){
        plt <- plt + theme(axis.title.y = element_blank())
    }

    if (names == F){
    
    plt <- plt + theme(axis.text.x = element_blank())

    }

    plt$layers[[2]]$aes_params$size <- 1
    plt$layers[[2]]$aes_params$label.size <- 16
    plt$layers[[2]]$aes_params$vjust <- 0.5
    

    if (points == T){

        plt <- plt + geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7)) 
    }

    ggsave(
        plot = plt, filename = paste0(
            outdir, "/", flag,  "_Warp_boxplot.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )


    return(plt)



}
#function:scatterplot
#input:
#      base: data.frame with samples as rows and features as columns
#      factor: Name of factor to analyse.
#      column: colname of the column to analyze.
#      group: column name of the variable to group the samples
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#               of the labels
#     levels: vector with levels in the order that should appear
#     wrap: True: divide samples in group using wrap function
#     remove.y: True remove name in y label
#     legend: True show legend


#output: # Plots of scatterplots indicating by color the levels of column
# and separate plots by group variable

#description: Function to create scatterplots

scatterplot <- function(base,
factor,
column,
group,
colors,
outdir,
flag,
levels,
wrap = T,
remove.y = F,
legend = F){

    base[,column] <- as.factor(base[,column])
    #base <- base[complete.cases(base[, column]), ]
    
    
    base[,column] <- droplevels(base[,column])
    base[, group] <- factor(base[,group],
   levels = levels)

  plt <- ggplot(base,
  aes(x=.data[[group]],
  y=.data[[factor]],
  fill=.data[[column]])) +

  theme_classic() +
  
  facet_wrap(as.formula(paste0("~ ", group)),
            drop = T, scales = "free_x", nrow =1) +


     geom_jitter(aes(colour = .data[[column]]),size = 5,
     alpha = 1, position = position_jitter(seed = 1)) +

   geom_hline(yintercept=0, linetype="dashed", linewidth=0.7, alpha=0.5) +
                
    scale_color_manual(values = colors[[column]]) 
        
if(wrap == T){
plt <- plt + theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.ticks.length.y = unit(10,"pt"))

        }else{
        plt <- plt + theme(text = element_text(size = 40, colour = "black" ),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.text.x = element_text(colour = "black",
        angle = 0, 
        vjust = 0,
        hjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.ticks.length.y = unit(10,"pt"))

        }
        
           if (remove.y == T){
        plt <- plt + theme(axis.title.y =  element_blank())
    }
       plt <- plt  +
        xlab("") +
    guides(color = guide_legend(override.aes = list(size = 7)))


    if (legend == T){
        plt <- plt + theme(legend.position = "right")
    }
            ggsave(
                plot = plt, filename = paste0(
                    outdir, "/", flag, "_Scatterplot.pdf"
                ),
                device = "pdf",
                dpi = 500, width = 10,
                height = 5
            )

return(plt)

}
#function:correlation.test
#input:
#      data: data.frame samples as rows and features as columns
#      Z: matrix with samples as rows and factors as columns
#      columns: colname of the columns to analyze.
#      flag: name to add to the name of output files
#      outdir: path where write results
#      p.th: threshold of p-value to show plots


#output: results of perason correlation  between factors and features
#description: Calculate pearson correlation between factors and continous columns

correlation.test <- function(data,
Z,
columns,
flag,
outdir,
p.th) {

    ano_df <- lapply(columns, function(x) {

        cdata <- data[complete.cases(data[, x]),]

        vf <- lapply(colnames(Z), function(y) {

            corr <- cor.test(cdata[, x],cdata[, y], method = "spearman")
            f <- corr$estimate[[1]]
            p <- corr$p.value
            

            return(c("estimate" = f, "p.value" = p))
        })

        names(vf) <- colnames(Z)
        vf <- as.data.frame(do.call(rbind, vf))
        colnames(vf) <- paste(x, colnames(vf))
        return(as.data.frame(vf))
    })

    ano_df <- as.data.frame(t(do.call(cbind, ano_df)))
 
    
    p.df <- ano_df[grep("p.value", rownames(ano_df)), ]


     for (col in colnames(p.df)) {

    p.df[, col] <- p.adjust(p.df[, col], method = "BH")

    }

    ano_df[grep("p.value", rownames(ano_df)), ] <- p.df

     write.table(ano_df, paste0(
        outdir , "/Correlation.txt"
    ),
    sep = "\t",
    row.names = T,
    col.names = NA
    )

     # we obtain only the p.value
    c.df <- ano_df[grep("estimate", rownames(ano_df)), ]
    rownames(c.df) <- gsub("\\.cor", "", rownames(c.df))

    # We will show only color result in those that are significant
    # We only show variables with any significance result

    c.df[p.df > p.th] <- 0

    return(c.df)
}





#function:heatmap.correlacion
#input:
#      mod: MOFA model
#      columns: colname of the columns to analyze.
#      flag: name to add to the name of output files
#      outdir: path where write results
#      p.th: threshold of p-value to show plots


#output: heatmap showing perason correlation  between factors and features
#description: Cgenerate heatmap showing results of perason correlation between
# factors and continous variables

heatmap.correlacion <- function(mod,
columns,
flag, outdir,
p.th = 0.05){

Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    df <- mod@samples_metadata

    p.df <- correlation.test(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )
    

    rownames(p.df) <- gsub(" estimate", "", rownames(p.df))

    col <- colorRamp2(
                c(1, 0, -1),
                c("red", "grey","blue" ))

    
    #Change names order
    hm1 <- Heatmap(as.matrix(t(p.df)),
            show_row_names = T,
            show_column_names = T,
            col = col,
            cluster_rows = F,
            cluster_columns = T,
            show_row_dend = F,
            show_column_dend = F,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            column_names_rot = 45,
            rect_gp = gpar(colo = "black", lwd =2),
            show_heatmap_legend = F
        )


lgd <- Legend( col_fun = colorRamp2(
                c(1, 0, -1),
                c("red", "grey","blue" )),
                title = "Correlación",
                title_position = "topcenter",
                direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm"))

pd = packLegend(lgd, direction = "horizontal")



   pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = 35, height= 35
        )
        
        draw(hm1, padding = unit(c(130, 10, 2, 20), "mm"))
        draw(lgd, x = unit(31, "cm"), y = unit(0, "cm"),
        just = c("left", "bottom"))


        dev.off()

}

#function:heatmap_classes
#input:
#      Z: matrix with samples as rows and factors as columns
#      base: data.frame with samples as rows and features as columns
#      outdir: path where write results
#      col: list indicating colors of anotation
#      size: vector indicating width in first term and height in secon term
#      factors.rm: factors that shoul remove from teh MOFA model
#      anot_param_legend: list indicating parameetrs of anotation
#      anot_h: anotation height
#      column: column to split samples
#      flag: name to add to the name of output files



#output: heatmap showing the subtypes and anotation data in col list and factors clustered
#using pearson correlation distance and hialerical clustering
#description: genearte heatmap to show factors of samples and anotated them
#by subtypes and clinical data

heatmap_classes <- function(Z, base,
outdir, col, size,
factors.rm,
anot_param_legend,
anot_h, column,
flag) {
    
Z <- Z[, !colnames(Z) %in% factors.rm]

Z <- Z[base$sample,]
data <- t(Z)
data <- zscore.rows2(data)

anot <- base[, names(col)]


ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = F,
            col = col,
            annotation_name_gp = gpar(fontsize = 120),
            border = unit(40, "cm"),
            annotation_height = anot_h,
            height = unit(30, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE)


hm1 <- Heatmap(as.matrix(data),
            top_annotation = ha_column,
            show_row_names = T,
            show_column_names = F,
            col = colorRamp2(
                c(max(data), 0, min(data)),
                c("red", "white", "blue")
            ),
            gap = unit(5, "mm"),
            cluster_rows = F,
            show_row_dend = F,
            show_column_dend = F,
            column_dend_side = "top",
            column_dend_height = unit(6, "cm"),
            cluster_columns = T,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "Z score (Factor)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 120),
                labels_gp = gpar(fontsize = 120,
                family = "Times", face = "bold"),
                legend_width = unit(50, "cm"),
                grid_height = unit(4.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 120),
            column_split = anot[, column],
            cluster_column_slices = FALSE,
            column_title  = c("","","")
        )

lgds <- packLegend(unlist(anot_param_legend),
direction = "horizontal")

hlist <- hm1

p <- draw(hlist, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        padding = unit(c(20, 80, 20, 20), "mm"),
        annotation_legend_list = lgds)

pdf(paste0(
            outdir, "/Heatmap_anotated_", flag, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()

}

#function:cut_survival
#input:
#      base: data.frame with samples as rows and clinical information in columns
#       cut: months where generate the final time to use
#      time: overall survival time from dfiagnosis to last interview
#      event : Overall survival status in 0 alive 1 dead
#      covariates: vector of colum names to adjust the model
#     
#output: base with time and event change cuting the information at months in cut 
# argument

#description:function to cut time and event to a certain number of months

cut_survival <- function(base, cut,
time = "OS.time",
event = "EXITUS" ){

base[,event][base[,time] > cut] <- 0
base[,time][base[,time] > cut] <- cut
return(base)
}

#function:Kapplan_Meyer
#input:
#      base: data.frame with samples as rows and clinical information in columns
#       class: column name of feature of interest
#      OS: overall survival time from dfiagnosis to last interview
#      OS_event : Overall survival status in 0 alive 1 dead
#      covariates: vector of colum names to adjust the model
#      outdir: Path to write results
#      flag: Name to add in the name of files
#      custcolor: vector with colors to use for levels in class variable
#      width: width of pdf
#      height: height of pdf
#      risk.table: True show risk table
#      group.by: facet by a column
#      font: size of font
#      fontsize: size of fontisize

#output: plot of Kapplan Meyer with the p-value  of log rank test and the
#survival probablity at 5 and 10 years

#description:function to create Kapplan Meyer plot

Kapplan_Meyer <- function(base,
class,
OS,
OS_event,
covariates = NULL,
outdir,
flag,
custcolor,
width= 8,
height=8,
risk.table = F,
group.by = NULL,
font = 23,
fontsize = 7
){
  colnames(base) <- gsub(OS, 'OS', colnames(base))
  colnames(base) <- gsub(OS_event, 'OS_event', colnames(base))
  base <- base[complete.cases(base$OS),]
  base <- base[complete.cases(base$OS_event),]
    # We asses that the data is numeric  
    base$OS <- as.numeric(base$OS)
    base$OS_event <- as.numeric(base$OS_event)
    base$Surv_OS <- with(base, Surv(base$OS, base$OS_event))

    #We generate survival formula
    form <- as.formula (paste( "Surv_OS ~ ", class))
  

    km.by.SG <- surv_fit (form,data = base , conf.type = "log-log")

        diff_SG <- survdiff(as.formula (paste( 'base$Surv_OS ~',
        class)),
        data = base, rho = 0)


        svc <- survfit(as.formula ( "Surv_OS ~ 1"),data = base , conf.type = "log-log")
        dfc <- summary(svc, c(60,120))
        dfc <-  data.frame("time" = dfc[[2]],
        "supervivencia" = dfc[[6]],
        "clase" = "Population")

        #calculate survival probability at 60 an 120 months
        sv <- survfit(form, data= base)
        df <- summary(sv, c(60,120))
        df <- data.frame("time" = df[[2]],
        "supervivencia" = df[[6]],
        "clase" = df[[10]])

        df <- rbind(df,dfc)
        write.table(df,
        sep = "\t",
        paste0(outdir, "/", flag, "_Supervivencia_tabla.txt"),
        row.names = F,
        col.names = T
        )

        pval <- signif(pchisq(diff_SG$chisq,
        length(diff_SG$n)-1,
        lower.tail = FALSE),
        digits = 3)
        
        pdf(paste0(outdir, "/", flag, "_",
        "Kapplan_Meyer.pdf"), width = width, height = height)

        p <- ggsurvplot(km.by.SG,
        data = base,
        pval = pval,
        risk.table = risk.table ,
        xlab = "Tiempo",
        ylab = "Probabilidad de Supervivencia Global",
        palette = custcolor,
        font.main = c(font, "bold"),
        font.x = c(font, "bold"),
        font.y = c(font, "bold"),
        font.caption = c(font, "bold"), 
        font.legend = c(font, "bold"), 
        font.tickslab = c(font, "bold"),
        pval.size = 8,
        fontsize = fontsize,
        facet.by = group.by) 

        print(p)
        dev.off()
        return(p)
}



#function:heatmap.molecular
#input:
#      exp: matrix with genes as rows and samples as columns
#      met: matrix with methylation probes as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      outdir: Path to write results
#      col: list indicating colors of anotation
#      size: vector indicating width in first term and height in secon term
#      anot_param_legend: list indicating parameetrs of anotation
#      flag: string to add in name of the file
#      max.exp: top number for color NULL select the maximum in the expression matrix
#      min.exp: minimum number for the color NULL select the minimum number in the expression matrix
#      max.met: top number for color NULL select the maximum in the matrix
#      min.met: minimum number for the color NULL select the minimum number in the matrix
#      names.size: size of rownames
#      column: column to split samples
#      anot.h: anotation height
#      annotation_name_gp: size of anotation names
#      show_row_names: True rownames is shown False row names are not shown
#      scale_met: True if scale methylation data
#      scale_exp: True if scale methylation data
#      height_anotation: height of anotation data
#      cluster: exp or met to select which omic use to applay clustering
#      order: order of the dendogram
#      nlegend: n legend in legend_param

#output: heatmap

#description: function to create a heatmap with expression and methylation data
heatmap.molecular <- function(exp,
met,
base,
outdir, col, size,
anot_param_legend,
flag,
max.exp = NULL,
min.exp = NULL,
max.met = NULL,
min.met = NULL,
names.size = 3,
column,
anot_h,
annotation_name_gp = 60,
show_row_names = T,
scale_met = T,
scale_exp = T,
height_anotation = 25,
cluster = "exp",
order,
nlegend = NULL) {

exp <- exp[, rownames(base)]


base[,column] <- as.factor(as.character(base[,column]))


if (scale_exp == T){
exp <- zscore.rows2(exp)
}



met <- met[, rownames(base)]

if (scale_met == T){
met <- zscore.rows2(met)
}

if (cluster == "exp"){
cluster_cols <- cluster_within_group(exp,base[,column])

}

if (cluster == "met"){
cluster_cols <- cluster_within_group(met,base[,column])
}

#order data as dendogram
base <- base[labels(cluster_cols),]
base <- base[c(grep(order[[1]],base[,column]),
grep(order[[2]],base[,column]),
grep(order[[3]],base[,column])),]
exp <- exp[, rownames(base)]
met <- met[, rownames(base)]
anot <- base[, names(col)]

ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = F,
            col = col,
            annotation_name_gp = gpar(fontsize = annotation_name_gp),
            border = unit(40, "cm"),
            annotation_height = anot_h,
            height = unit(height_anotation, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            annotation_legend_param = anot_param_legend,
            simple_anno_size = unit(4, "cm"))

if (is.null(max.exp)){

    max.exp <- max(exp)
}
if (is.null(max.met)){

    max.met <- max(met)
}

if (is.null(min.exp)){

    min.exp <- min(exp)
}

if (is.null(min.met)){

    min.met <- min(met)
}

col.exp <- colorRamp2(
                c(max.exp, 0, min.exp),
                c("red", "white", "blue")
)



col.met <- colorRamp2(
                c(max.met, 0, min.met),
                c("green", "white", "blue")
)

hm1 <- Heatmap(as.matrix(exp),
            top_annotation = ha_column,
            show_row_names = show_row_names,
            show_column_names = T,
            col = col.exp,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "Z score (expresión)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 80),
                labels_gp = gpar(fontsize = 80,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
#            column_split = anot[, column],
            height = unit(40, "cm")

        )

hm2 <- Heatmap(as.matrix(met),
            show_row_names = show_row_names,
            show_column_names = F,
            col = col.met,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "Z score (metilación)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 80),
                labels_gp = gpar(fontsize = 80,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
            height = unit(40, "cm")
        )


if (nlegend == 4){

   lgds <-  packLegend(anot_param_legend[[1]],
    anot_param_legend[[2]],
    anot_param_legend[[3]],
    anot_param_legend[[4]],
    direction = "horizontal",
    column_gap = unit(9, "cm")

    )
}else{

lgds <- packLegend(unlist(anot_param_legend))
}
hlist <- hm1 %v% hm2

p <- draw(hlist, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        ht_gap = unit(3, "cm"),
        padding = unit(c(20, 20, 20, 20), "mm"),
         annotation_legend_list = lgds)

pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()

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



deg.boxplot <- function(data,
base,
group.x,
group.fill,
color.fill,
cov = NULL,
flag,
outdir,
genes,
filter.table = T) {



genes <- gsub("-","_",genes)
colnames(data) <- gsub("-", "_", colnames(data))
genes <- gsub("\\.","_",genes)
colnames(data) <- gsub("\\.", "_", colnames(data))



genes <- genes[genes %in% colnames(data)]


data <- data[rownames(base),,drop = F]

df<- cbind(data,base[,c(group.x,group.fill)])
df[,group.x] <- as.factor(as.character(df[,group.x]))
df[,group.fill] <- as.factor(as.character(df[,group.fill]))


if (group.x != group.fill){

base[,group.x] <- as.factor(as.character(base[,group.x]))
base[,group.fill] <- as.factor(as.character(base[,group.fill]))

degs <- lapply(levels(base[,group.x]), function(x){

dfp <- base[base[,group.x] == x,]

degp <- deg.function(as.data.frame(t(data))[,rownames(dfp)],
dfp,
group = group.fill,
cov = cov,
outdir = outdir,
flag = x)

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

if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}
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

write.table(degs,
paste0(outdir, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)


}else{

base[,group.fill] <- as.factor(as.character(base[,group.fill]))


degp <- deg.function(as.data.frame(t(data)),
base,
group = group.fill,
cov = cov,
outdir = outdir,
flag = "Complete")
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

if (filter.table == T) {
fdegp <- fdegp%>%
filter(
  gene %in% genes
)
}
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
paste0(outdir, "/Limma_genes_primary_tumor_Subtype.txt"),
sep = "\t",
row.names = F,
col.names = T)

}

plots <- list()

for (row in genes) {


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

stat.test <- stat.test %>%
mutate(
  p.adj.signif = stat.deg$p.adj.signif
)

stat.test[[group.fill]] <- NA


bxp <- ggplot(
  df, aes(x = .data[[group.x]],
  y = .data[[row]], fill= .data[[group.fill]])) +
  geom_boxplot(width = 0.6, outlier.shape = NA,
  aes(fill= .data[[group.fill]])) +

  scale_fill_manual(values = as.list(c(color.fill,"blue", "green"))) +

  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7)) +
  
  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", hide.ns = T,  step.increase = 0.1
) +


  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme_classic() +

theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.ticks.length.y = unit(10,"pt")
        ) +

        labs(x = "", y = row)


bxp$layers[[3]]$aes_params$size <- 1.2
bxp$layers[[3]]$aes_params$label.size <- 12 

ggsave(
        plot = bxp, filename = paste0(
            outdir, "/Boxplot.subtype_", row, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = 16,
        height = 12
    )


 plots[[row]] <- bxp
   

}





return(plots)
}


#function:correlation.plot
#input:
#      exp: transcriptomic matrix with genes s rows and samples as columns
#      met: methylation matrix with probes as rows and samples as columns
#      anot: anot matrix for methylation data
#      base: data.frame with samples as rows and features as columns
#      genes: genes to apply the corelation
#      probes: methylation probes to applay de correlation
#      outdir: Path to write results
#      group: Column name to facet the correlation plots
#      flag: string o add in file names
#      stat.size: size of the pearson correlation results

#output: boxplots with limma p value and a table with ther esults of limma 

#description: function to obtain boxplots showing the results of limma


correlation.plot <- function(exp,
met,
anot,
base,
genes,
probes,
outdir,
group,
flag,
stat.size = 6) {


#Eliminate - in gene names

genes <- gsub("-","_", genes)
rownames(exp) <- gsub("-","_", rownames(exp))
anot$genes <- gsub("-","_", anot$genes )

exp <- exp[genes,]
met <- met[probes,]

exp <- exp[,rownames(base)]
met <- met[,rownames(base)]



crp <- list()

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


plt2 <- ggplot(df, aes_string(x = gene,
    y = gene_met)) +
    theme_classic() + 
    geom_smooth(method = "lm")  +
    geom_point(size = 3) +
    facet_wrap(as.formula(paste0("~ ", group)),
            drop = T, scales = "free_x", nrow =1
        ) +
    stat_cor(label.y= max(df[,gene_met]) + max(df[,gene_met]) / 5, label.x = Inf,
    vjust = 0.5, hjust = 1.1,
    size = stat.size) +

    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(margin = margin(20, 10, 20, 20),
        colour = "black"),
        axis.ticks.length.y = unit(10,"pt")
    ) 


ggsave(
        plot = plt2, filename = paste0(
            outdir,"/", gene, "_", probe,"_", group, ".wrapped.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 25,
        height = 12
    )

crp[[gene]] <- plt2

}
}
return(crp)
}

#function:generate_medians 
#input:
#      data: matrix with features to add the model as rows and samples as columns 
#      features: features to generate medians

#output: data.frame with samples in rows an features in columns indicating if
# samples are up or down the median of the feature
#description:function to obtain which sampels are up or down the median of 
#some genes
generate_medians <- function(data,features){

    medians <- lapply(features,function(x){

        med <- median(as.numeric(data[x,]))
        
        result <- ifelse(as.numeric(data[x,])> med, "Alto", "Bajo")
        return(result)
    })


    names(medians) <- features
    medians <- do.call(cbind, medians)
    rownames(medians) <- colnames(data)
    return(medians)

}


#function:Hazard_Cox
#input:
#      base: data.frame with samples as rows and clinical information in columns
#      class: column name of feature of interest
#      OS: overall survival time from dfiagnosis to last interview
#      OS_event : Overall survival status in 0 alive 1 dead
#      covariates: vector of colum names to adjust the model
#      outdir: Path to write results
#      flag: Name to add in the name of files

#output: table with hazard ratio the interval and p-value of each level in each variable

#description:function to create Cox models

Hazard_Cox <- function(base,
class,
OS,
OS_event,
covariates = NULL,
outdir,
flag
) {

  colnames(base) <- gsub(OS, 'OS', colnames(base))
  colnames(base) <- gsub(OS_event, 'OS_event', colnames(base))

    # We asses that the data is numeric  
    base$OS <- as.numeric(base$OS)
    base$OS_event <- as.numeric(base$OS_event)
    #We generate the rurvival objectbase
    base$Surv_OS <- with(base, Surv(base$OS, base$OS_event == 1))
    if (is.null(covariates)){
    fit <- coxph( as.formula(paste('Surv_OS ~ ', class)), data = base)

    }else{
      fit <- coxph( as.formula(paste('Surv_OS ~ ', class,' + ',paste(covariates, collapse = ' + '))) , data = base)
    }
    sum <- summary(fit)
    coef <- sum[["coefficients"]]

    
    coef_model <- coef(fit)
    cov <- vcov(fit)


   global_p <- lapply(c(class,covariates), function(x) {
        
        indices <- grep(paste0("^",x),names(coef_model))

        wald <- wald.test(b = coef_model[indices],
        Sigma = as.matrix(cov[indices,indices]),
        L = diag(length(indices)))
        
        return(rep(wald$result$chi2[["P"]],length(indices)))
    })


    sur.df <- data.frame(
        cluster = factor(rownames(coef), levels = rev(rownames(coef))),
        p      = coef[, "Pr(>|z|)"],
        coef   = coef[, "exp(coef)"],
        lower  = sum[["conf.int"]][, "lower .95"],
        higher = sum[["conf.int"]][, "upper .95"],
        global_p = unlist(global_p)
    )

    write.table(sur.df,
        sep = "\t",
        paste0(outdir, "/", flag, "_Hazard_ratio_factors.txt"),
        row.names = F,
        col.names = T
    )
    
    ggplot(sur.df, aes(x = cluster, y = coef, ymin = lower, ymax = higher)) +
        geom_pointrange(col = "#619CFF", size = 1.3) +
        coord_flip() +
        scale_x_discrete() +
    labs(y = "Hazard Ratio", x = "") +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        theme_bw() + 
        theme(text = element_text(size = 26))

    ggsave(
        filename = paste0(outdir, "/", flag, "Hazard_ratio_factors.pdf"),
        device = "pdf",
        dpi = 900, width = 10, height = 10
    )
    return <- sur.df
    }