#Functions adictional for MOFA analisis

#function:obtain_Z
#input:
#      mod: MOFA model
#       
#output: matrix with samples in rows and factors in columns

#description: Function to obtain the value of Factors in each sample
# 


obtain_Z <- function(mod){
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]
    return(Z)
}

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
        outdir,
        "/discrete.variables/", flag, "/Kruskal_Features_Factors_", flag, ".txt"
    ),
    sep = "\t",
    row.names = T,
    col.names = NA
    )

    # We will show only color result in those that are significant
    # We only show variables with any significance result

    p.df[p.df > p.th] <- 1
    p.df <- p.df[apply(p.df, 1, sum) != nrow(p.df), ]

    col <- colorRampPalette(c("lightgrey", "red"))(n = 100)

    rownames(p.df) <- gsub(" p.value", "", rownames(p.df))

    pheatmap::pheatmap(-log(t(p.df)),
        cluster_rows = F, color = col,
        filename = paste0(
            outdir,
            "/discrete.variables/", flag, "/Heatmap_Kruskal_", flag, ".pdf"
        )
    )

    return(p.df)
}

#function:wrap.boxplot
#input:
#      mod: MOFA model or clinical data with Factor information.
#      factor: Name of factor to analyse.
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      column: colname of the column to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#               of the labels
#      is.model: True mod is MOFA model or False if it is clinical data.frame
#      group: column name to stratify samples
#      levels: Label order for de group column.  
#output: # Plots of boxplots separating samples by the group variable that compare
# the factor between the labels of the column feature and the result of the
#Wilcoxon test

#description: Function to create boxplots that compares the levels of the factor
# between the labels of the column feature in each of the labels of the group column
#separately

wrap.boxplot <- function(mod, column,
factor,
flag,
colors,
outdir,
is.model = TRUE,
group ,
levels )
 {
    
    # Obtain clinical data.frame
    if (is.model == T){
    cli.df <- mod@samples_metadata
    }else{
        cli.df <- mod
    }


    # we cant't have spaces in levels
    cli.df[, column] <- gsub(" ", "_", cli.df[, column])
    cli.df[, group] <- gsub(" ", "_", cli.df[, group])
    levels <- gsub(" ", "_", levels)
    cli.df[, group] <- factor(cli.df[,group],
    levels = levels)

    # We eliminate missing values
        
        df <- table(base[,group], base[,column])
        l <- rownames(df)[df[,1] == 0 | df[,2] == 0]
        if(length(l) > 0){
        cli.df <- cli.df[cli.df[,group] != l, ]
        }

    cli.df[,column] <- as.factor(cli.df[,column])
    cli.df <- cli.df[complete.cases(cli.df[, column]), ]
    cli.df[,column] <- droplevels(cli.df[,column])

    # Obtain results of Wilcoxon test
    stat_df <- compare_means(as.formula(paste0(factor, " ~ ", column)),
        group.by = group, data = cli.df, p.adjust.method = "fdr"
    )

    # Calculate effect size
    we <- cli.df  %>% 
    group_by_at(vars(group)) %>%
    wilcox_effsize(as.formula(paste0(factor, " ~ ", column)))
    
    we$group <- gsub(" ", "_", we$group)
    
    we$names <- paste(we$.y., 
    we$group1,
    we$group2,
    we$group,
    sep = "_")

    stat_df$names <- paste(stat_df$.y., 
    stat_df$group1,
    stat_df$group2,
    stat_df[, group][[1]],
    sep = "_")
    
    # Add effect size to the results of Wilcoxon test
    stat_df <- merge(stat_df, we[,c(4:9)], by = "names")

# Generate boxplot
    plt <-
        ggboxplot(cli.df,
            y = factor,
            x = group, fill = column,
            width = 0.6
        ) +

        scale_fill_manual(values = colors[[column]]) +

        theme_classic() +
        
        # add p values
        geom_pwc(
            aes(group = .data[[column]]),
            method = "wilcox.test", label = "p.adj.signif",
            hide.ns = T,
            method.args = list(p.adjust.method = "BH", exact = FALSE),
        ) +

        # separate by group column
        facet_wrap(as.formula(paste0("~ ", group)),
            drop = T, scales = "free_x", nrow =1) +

        # change sizes
        theme(
            text = element_text(size = 40),
            axis.text.x = element_blank(),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 40),
            strip.background = element_rect(
            color="black", size=3, linetype="solid"
     )
        )

    plt$layers[[2]]$aes_params$size <- 1
    plt$layers[[2]]$aes_params$label.size <- 26


    ggsave(
        plot = plt, filename = paste0(
            outdir, "/discrete.variables/",
            flag, "/", column, "/Factors_Boxplot_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )


    return(stat_df)
}

#function:complete.boxplot
#input:
#      mod: MOFA model or clinical data with Factor information.
#      factor: Name of factor to analyse.
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      column: colname of the column to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#               of the labels
#      is.model: True mod is MOFA model or False if it is clinical data.frame

#output: # Plots of boxplots of all samples that compares
# the factor between the labels of the column feature and the result of the
#Wilcoxon test

#description: Function to create boxplots that compares the levels of the factor
# between the labels of the column feature.

complete.boxplot <- function(mod, column,
factor,
flag,
colors,
outdir,
is.model = T) {
    
    if (is.model == T){
    cli.df <- mod@samples_metadata
    }else{
        cli.df = mod
    }
    # we cant't have spaces in levels
    cli.df[, column] <- gsub(" ", ".", cli.df[, column])
    names(colors) <- gsub(" ", ".", names(colors))

    #Regenerate factors
    cli.df[, column] <- factor(cli.df[, column],
            levels = names(colors[[column]])
    )
    # We eliminate missing values

    cli.df[,column] <- as.factor(cli.df[,column])
    cli.df <- cli.df[complete.cases(cli.df[, column]), ]
    cli.df[,column] <- droplevels(cli.df[,column])
    
    stat_df <- compare_means(as.formula(paste0(factor, " ~ ", column)),
    data = cli.df, p.adjust.method = "BH"
    )
    
    we <- wilcox_effsize(data = cli.df,
    as.formula(paste0(factor, " ~ ", column)))
    
    we$names <- paste(we$.y., 
    we$group1,
    we$group2,
    sep = "_")

    stat_df$names <- paste(stat_df$.y., 
    stat_df$group1,
    stat_df$group2,
    sep = "_")
    
    stat_df <- merge(stat_df, we[,c(4:8)], by = "names")

    plt <-
        ggboxplot(cli.df,
            y = factor,
            x = column, fill = column,
            width = 0.6
        ) +

        scale_fill_manual(values = colors[[column]]) +

        theme_classic() +
        geom_pwc(
            aes(group = .data[[column]]),
            method = "wilcox.test", label = "p.adj.signif",
            hide.ns = T,
            method.args = list(p.adjust.method = "BH", exact = FALSE),
        ) +

        theme(
            text = element_text(size = 40),
            axis.text.x = element_blank(),
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 40)
            
        )

    plt$layers[[2]]$aes_params$size <- 1
    plt$layers[[2]]$aes_params$label.size <- 26


    ggsave(
        plot = plt, filename = paste0(
            outdir, "/discrete.variables/",
            flag, "/", column, "/Factors_Boxplot_Complete_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )

    return(stat_df)
}

#function:correlation.test
#input:
#      data:  clinical data.frame.
#      Z: matrix with samples as rows and factors as columns
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      columns: colnames of the columns to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      

#output: # Result of pearson correlation between factor and continous factor

#description: Function to create heatmap of spearman correlations between factors and continous features
# and the data.frame with this information

correlation.test <- function(data, Z, columns, flag, outdir, p.th) {

    ano_df <- lapply(columns, function(x) {

        cdata <- data[complete.cases(data[, x]),]

        # correlation between each factor and each variable
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
 
    # Obtain pnly p values
    p.df <- ano_df[grep("p.value", rownames(ano_df)), ]


     for (col in colnames(p.df)) {
    
    # adjust p-values
    p.df[, col] <- p.adjust(p.df[, col], method = "BH")

    }

    ano_df[grep("p.value", rownames(ano_df)), ] <- p.df

     write.table(ano_df, paste0(
        outdir,
        "/continous.variables/", flag,
        "/Correlation_Features_Factors_", flag, ".txt"
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

breaksList <- seq(from = -1 , to =1, by = 0.01)

col <- colorRampPalette(c("blue", "lightgrey", "red"))(length(breaksList))

    rownames(c.df) <- gsub(" estimate", "", rownames(c.df))

    # generate heatmap
    pheatmap::pheatmap(t(c.df),
        cluster_rows = F, color = col, breaks = breaksList,
        
        filename = paste0(
            outdir,
            "/continous.variables/", flag, "/Heatmap_Correlation_", flag, ".pdf"
        )
    )

    return(c.df)
}

#function:correlation.plot
#input:
#      mod: MOFA model or clinical data with Factor information.
#      factor: Name of factor to analyse.
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      column: colname of the column to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      group: column name to stratify samples
#      levels: Label order for de group column.  

#output: # Correlation plots between factors and variables using all sampples and 
# separate samples by group column.

#description: Function to generate correlation plots between factors and variables using all sampples and 
# separate samples by group column.


correlation.plot <- function(mod, column, factor,
flag, outdir, group, levels) {


 cli.df <- mod@samples_metadata

    cli.df <- cli.df[complete.cases(cli.df[, column]), ]
    cli.df[,group] <- as.factor(cli.df[,group])

plt1 <- ggplot(cli.df, aes_string(x = column, y = factor, group = 1)) +
    theme_classic() +
    geom_point(size = 3) +

    stat_cor(size = 10) +

    stat_smooth(method = "lm", color = "black") +
    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
    )

    ggsave(
        plot = plt1, filename = paste0(
            outdir, "/continous.variables/",
            flag, "/", column, "/Factors_Corr_plot_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 17,
        height = 12
    )


    plt2 <- ggplot(cli.df, aes_string(x = column, y = factor)) +
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
            outdir, "/continous.variables/",
            flag, "/", column, "/Factors_Corr_plot_wrapper_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 22,
        height = 12
    )

}





#function:Anot_path
#input:
#      df: matrix with genes as rownamaes or a column with genes to anotate in pathways
#      gmt: gmt file 
#     

#output: # df matrix with pathway anotations

#description:Function to anotate the genes of df matrix

Anot_path <- function(df, gmt) {
    
    genes <- unlist(gmt)
    genes <- unique(genes)

    genes_l <- list()

    for (x in genes) {
        for (y in names(gmt)) {
            if (x %in% gmt[[y]]) {
                genes_l[[x]] <- append(genes_l[[x]], y)
            }
        }
    }

    genes_l <- lapply(genes_l, toString)


    genes <- do.call(rbind, genes_l)

    colnames(genes) <- "Pathways"
    genes <- as.data.frame(genes)

    if ("genes" %in% colnames(genes) == F) {
        genes$genes <- rownames(genes)
    }

    if (!"genes" %in% colnames(df)) {
        df$genes <- rownames(df)
    }
    merge(df, genes, by.x = "genes", by.y = "genes", all.x = T)
}

#function:Anot_path
#input:
#      genes_l: list of genes to be anotated
#      gmt: gmt file 
#     

#output: # df matrix with pathway anotations

#description:Function to anotate the genes of genes_l

MOFA_list2DF <- function(genes_l, gmt) {

    df_all <- as.data.frame(unlist(genes_l))
    features <- lapply(rownames(df_all), function(x) {
        strsplit(x, "\\.")
    })
    factor <- unlist(lapply(features, function(x) {
        return(x[[1]][1])
    }))
    omic <- unlist(lapply(features, function(x) {
        return(x[[1]][2])
    }))
    omic <- gsub("[0-9]*", "", omic)

    df_all$Factor <- factor
    df_all$Omic <- omic
    colnames(df_all) <- c("genes", colnames(df_all)[2:ncol(df_all)])
    rownames(df_all) <- NULL
    df_all <- transform(df_all, lab = ave(Omic, genes, Factor, FUN = toString))

    df_all <- df_all[, c(1, 2, 4)]
    df_all <- df_all[!duplicated(df_all), ]

    df_all <- dcast(data = df_all, formula = genes ~ Factor, value.var = "lab")

    df_all$Total <- apply(df_all, 1, function(x) {
        return(length(x[complete.cases(x)]))
    })


    df_all <- Anot_path(df_all, gmt)

    df_all <- df_all[order(df_all$Total, decreasing = T), ]
    df_all$Total <- df_all$Total - 1
    rownames(df_all) <- df_all$genes

    df_all <- df_all[, -1]

    return(df_all)
}


#function:Corr_factor
#input:
#      omic: matrix indicating in rows features and in columns samples
#      Z: matrix with rows samples and in columns factors
#     

#output: # Correlation between genes and factors
#description: Function to correlate genes and factors

Corr_factor <- function(omic, Z, f) {
    crexp <- lapply(as.data.frame(t(omic)), function(x) {
        crt <- cor.test(x, Z[, f])
        return(list("cor" = crt$estimate[[1]], "p.value" = crt$p.value[[1]]))
    })

    dfcorr <- as.data.frame(do.call(rbind, crexp))
    dfcorr$cor <- unlist(dfcorr$cor)
    dfcorr$p.value <- unlist(dfcorr$p.value)
    return(dfcorr)
}



#function: gsea
#input:
#      rnk: ranking of genes
#      gmt: gmt file with pathway information
#      title: name of the output file
#     

#output: # gsea of the ranking in rnk
#description: Function to applay gsea over a gene ranking base 
# on correlatioon of factors

gsea <- function(rnk, gmt, title) {
    
    gmt <- gmtPathways(gmt)
    gsea <- fgsea(
        pathways = gmt,
        stats = rnk,
        minSize = 15,
        maxSize = 500
    )

    gsea$leadingEdge <- unlist(lapply(gsea$leadingEdge, function(x) {
        paste(x, collapse = "/")
    }))
    

    if (nrow(gsea[gsea$padj < 0.05, ]) > 0) {
        n <- nrow(gsea[gsea$padj < 0.05, ])

        if (n > 10) {
            n <- 10
        }

        gsea <- gsea[which(gsea$padj < 0.05), ]
        gsea <- gsea[order(abs(gsea$padj)), ]
        gsea$pathway <- factor(gsea$pathway, levels = gsea$pathway)

        ggplot(data = gsea[1:n, ], aes(x = NES, y = pathway, size = size,
        color = -1 * log10(padj))) +
            geom_point() +
            scale_size_area(max_size = 10) +
            scale_colour_gradient(low = "green", high = "red") +
            theme(text = element_text(size = 26))


        ggsave(
            filename = title,
            device = "pdf", dpi = 500, width = 22,
        height = 12
        )

    }

    return(gsea)
}