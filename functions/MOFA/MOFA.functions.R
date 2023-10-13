# Functions to create and analisis MOFA models
#function:Mofa.model
#input:
#      input: named list with the input matrix where samples 
#             are columns and featurs are rows
#       file: file where save the model generated
#      convergence_mode: convergence mode slow or fast
#      groups : To use group mode introduce a factor indicating the levels
#               of samples in the same order in which they are in the input
#               matrix

#output: MOFA model

#description: function to create the MOFA model


Mofa.model <- function(input,
                       file,
                       num_factors,
                       convergence_mode,
                       groups = NULL) {
    
    
    #Create MOFA model with or without groups
    if (is.null(groups) == T) {
        mofa1 <- create_mofa(input)
    } else {
        mofa1 <- create_mofa(input, groups = groups)
    }
    
    
    # Change default model options
    ModelOptions <- get_default_model_options(mofa1)
    ModelOptions$num_factors <- num_factors
    data_Options <- get_default_data_options(mofa1)

    TrainOptions <- get_default_training_options(mofa1)
    TrainOptions$seed <- 57
    TrainOptions$convergence_mode <- convergence_mode

    mofa1 <- prepare_mofa(
        mofa1,
        data_options = data_Options,
        model_options = ModelOptions,
        training_options = TrainOptions
    )

    # Create model
    mofa1 <- run_mofa(mofa1, outfile = file, use_basilisk = TRUE)
}


# Functions to create and analisis MOFA models
#function:anotate.metilation
#input:
#      model: MOFA model
#       anotation: anotation matrix with information 
#                  of genes in the methylated CpG

#output: MOFA model with genes anotated in CpG names

#description: We add genes anotated in the probes names of CpG

anotate.metilation <- function(model, anotation) {
    
    updated_features_names <- features_names(model)
    mt <- anotation[updated_features_names$Methylation, ]
    mt <- as.data.frame(mt)

    # Chage feature names adding gene symbol
    updated_features_names[["Methylation"]] <- as.vector(paste(mt$ID,mt$genes,
        sep = "_"
    ))

  features_names(model) <- updated_features_names

return(model)
}

#function:calculate.mean.variance
#input:
#      model: MOFA model
#       outdir: output path to write results
#output: matrix with factors as rows and views as colum indicating the mean
#        variance explained by each factor in each view of all tissues.

#description: To obtain the mean variance explained of each view by each factor
#             by the model

calculate.mean.variance <- function(model, outdir) {

    ve <- calculate_variance_explained(model)

    ve <- Reduce(rbind, ve$r2_per_factor)
    factors <- colnames(model@expectations$Z[[1]])

    fve <- lapply(factors, function(x) {
        apply(ve[rownames(ve) == x, ], 2, mean)
    })


    fve <- do.call(rbind, fve)
    rownames(fve) <- factors


    write.table(fve,
        paste0(outdir, "/Variance_Mean_views.txt"),
        sep = "\t",
        row.names = T,
        col.names = NA
    )
}

#function:discrete.features
#input:
#      mod: MOFA model
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      columns: colnames of variables to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#      of the variables
#      select: True or False. True to select variables which
#      kruskal wallis p-value less than threshold
#      group: column name to stratify samples
#      levels: Label order for de group column.  
#output: files and plot results of the Kruskall wallis and Wilcoxon test
# to analyze the associations between discrete variables and model
# factors

#description: Function to study the associations between the discrete variables
# and model factors.

discrete.features <- function(mod,
p.th = 0.05,
columns,
outdir,
flag,
colors,
select = T,
group = "TUMOR_PRIMARIO",
levels = c("INTESTINO", "PULMÓN", "COLORRECTAL", "ESTÓMAGO", "PÁNCREAS")) {

if (dir.exists(paste0(outdir, "/discrete.variables/", flag)) == F) {

    dir.create(paste0(outdir, "/discrete.variables/", flag))

}

    # Obtain factors
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    df <- mod@samples_metadata

    #Kruskal wallis results
    p.df <- kruskal(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

    # PLOT SIGNIFICANT ASSOCIATED FACTOR

    num <- apply(p.df,1,sum)
    if (select == T){
    columns <- rownames(p.df)[num != ncol(p.df)]
    }else{
    columns <- rownames(p.df)
    }
    for (column in columns) {

        if (dir.exists(paste0(
            outdir, "/discrete.variables/",
            flag, "/", column
        )) == F) {
            dir.create(paste0(
                outdir, "/discrete.variables/",
                flag, "/", column
            ))
        }
        mod@samples_metadata[, column] <- factor(mod@samples_metadata[, column],
            levels = names(colors[[column]])
        )
        
        if (select == T) {
        factors <- colnames(p.df)[p.df[column, ] < 1]
        }else{
        factors <- colnames(p.df)
        }

        # Plot factors all toghether
        plt <- plot_factor(mod,
            factors = as.numeric(gsub("Factor", "", factors)),
            color_by = column
        ) +

            scale_fill_manual(values = colors[[column]]) +

            theme(
                text = element_text(size = 18),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            ) +
            guides(color = guide_legend(override.aes = list(size = 15)))

        ggsave(
            plot = plt, filename = paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factors_Scatterplot_",
                flag, "_", column, ".pdf"
            ),
            device = "pdf",
            dpi = 500, width = 15,
            height = 5
        )



        # Plot factors each one separately

        stat_f <- list()
        stat_c <- list()
        for (factor in factors) {
            plt <- plot_factor(mod,
                factors = factor,
                color_by = column
            ) +

                scale_fill_manual(values = colors[[column]]) +

                theme(
                    text = element_text(size = 18),
                    axis.text.x = element_text(angle = 90,
                    vjust = 0.5, hjust = 1)
                ) +
                guides(color = guide_legend(override.aes = list(size = 15)))

            ggsave(
                plot = plt, filename = paste0(
                    outdir, "/discrete.variables/",
                    flag, "/", column, "/Factors_Scatterplot_",
                    flag, "_", column, "_", factor, ".pdf"
                ),
                device = "pdf",
                dpi = 500, width = 10,
                height = 5
            )

            # Plot boxplot stratified by group and obtai results of Wilcoxon
            stat_df <- wrap.boxplot(
                mod = mod, column = column, factor = factor,
                flag = flag, colors = colors, outdir = outdir,
                group = group, levels = levels
            )

            stat_f[[factor]] <- as.data.frame(stat_df)

            # Plot boxplot using all samples togheter and obtain results of Wilcoxon

            stat_complete <- complete.boxplot(
                mod = mod, column = column, factor = factor,
                flag = flag, colors = colors, outdir = outdir
            )

            stat_c[[factor]] <- as.data.frame(stat_complete)


        }

        
        # Create data.frame fom list with the results of all selected variables
        stat_f <- do.call(rbind, stat_f)
        stat_c <- do.call(rbind, stat_c)

        write.table(stat_f,
            paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factor_Wilcox_", flag, "_", column, ".txt"
            ),
            sep = "\t",
            row.names = F,
            col.names = T
        )

        write.table(stat_c,
            paste0(
                outdir, "/discrete.variables/",
                flag, "/", column, "/Factor_Wilcox_Complete_", flag,
                "_", column, ".txt"
            ),
            sep = "\t",
            row.names = F,
            col.names = T
        )
    }
}


#function:continous.features
#input:
#      mod: MOFA model
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      columns: colnames of variables to analyze.
#      outdir: Path to write outputs
#      flag: name to add to the name of output files
#      colors: list of named vector indicating the color
#      of the variables
#      select: True or False. True to select variables which
#      person correlation p-value less than threshold.
#      group: column name to stratify samples
#      levels: Label order for de group column.  
#output: files and plot results of the Kruskall wallis and Wilcoxon test
# to analyze the associations between discrete variables and model
# factors

#description: Function to study the associations between the discrete variables
# 
continous.features <- function(mod,
p.th = 0.05,
columns,
outdir,
flag,
group = "TUMOR_PRIMARIO",
levels = c("COLORRECTAL", "ESTÓMAGO", "PULMÓN", "PÁNCREAS", "INTESTINO")) {

if (dir.exists(paste0(outdir, "/continous.variables/", flag)) == F) {

    dir.create(paste0(outdir, "/continous.variables/", flag))

}
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    df <- mod@samples_metadata

    c.df <- correlation.test(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

    rownames(c.df) <- gsub("\\.rho", "", rownames(c.df))

    factors <- colnames(Z)

    # PLOT SIGNIFICANT ASSOCIATED FACTOR


    
    num <- apply(c.df,1,sum)
    columns <- rownames(c.df)[num != 0]

    for (column in columns) {

        if (dir.exists(paste0(
            outdir, "/continous.variables/",
            flag, "/", column
        )) == F) {
            dir.create(paste0(
                outdir, "/continous.variables/",
                flag, "/", column
            ))
        }

        # Plot all toghether

        plt <- plot_factor(mod,
            factors = as.numeric(gsub("Factor", "", factors)),
            color_by = column
        ) +
            theme(
                text = element_text(size = 18),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            )

        ggsave(
            plot = plt, filename = paste0(
                outdir, "/continous.variables/",
                flag, "/", column, "/Factors_Scatterplot_",
                flag, "_", column, ".pdf"
            ),
            device = "pdf",
            dpi = 500, width = 15,
            height = 5
        )



        # Plot each one separately


        for (factor in factors) {

        if (abs(c.df[column, factor]) > 0){

            plt <- plot_factor(mod,
                factors = factor,
                color_by = column
            ) +
                theme(
                    text = element_text(size = 18),
                    axis.text.x = element_text(angle = 90,
                    vjust = 0.5, hjust = 1)
                ) 

            ggsave(
                plot = plt, filename = paste0(
                    outdir, "/continous.variables/",
                    flag, "/", column, "/Factors_Scatterplot_",
                    flag, "_", column, "_", factor, ".pdf"
                ),
                device = "pdf",
                dpi = 500, width = 10,
                height = 5
            )


        #Plot significant variables
            correlation.plot(
                mod = mod, column = column, factor = factor,
                flag = flag, outdir = outdir, group = group, levels = levels
            )

        }

        }

    }
}

#function:heatmap_continous
#input:
#      mod: MOFA model
#      p.th: p-value threshold to select variables in Kruskal Wallis.
#      columns: colnames of variables to show in .
#      outdir: Path to write outputsheatmap
#      flag: name to add to the name of output files
#      colorv: list of named vector indicating the color
#      of the variables
#      size: size of the pdf 
#output: heatmap showing samples as columns and variables from columns input as rows.
# ordered by factors and

#description: Function to create the correlation heatmap between factors and variables.

# 

heatmap_continous <- function(mod, columns,
                              flag, outdir, colorv, size, anot_legend_param) {
    

    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    dbase <- mod@samples_metadata

    factors <- colnames(dbase)[grep("Factor", colnames(dbase))]

    # create anotation

    for (factor in factors) {


        # We order all by factor


        Z <- Z[order(Z[, factor]), ]


        df <- dbase[match(rownames(Z), dbase$sample), ]

        data <- df[, columns]
        rownames(data) <- df$sample

        data <- t(data)
        
        data <- zscore.rows2(data)


        col <- colorv


        col[[factor]] <- colorRamp2(
            c(
                min(df[, factor], na.rm = T),
                0,
                max(df[, factor], na.rm = T)
            ),
            c("blue", "white", "red")
        )



    names(anot_legend_param) <- c(names(anot_legend_param)[1:length(anot_legend_param)-1],
            as.character(factor))

        anot <- df[, names(col)]
        rownames(anot) <- rownames(df)


        ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = 40),
            border = unit(40, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 40),
            annotation_legend_param = anot_legend_param)




        hm <- Heatmap(as.matrix(data),
            top_annotation = ha_column,
            show_row_names = T,
            show_column_names = F,
            col = colorRamp2(
                c(max(data), 0, min(data)),
                c("red", "white", "blue")
            ),
            name = "Heatmap",
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "spearman",
            heatmap_legend_param = list(
                title = "Absolute infiltration",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 36),
        )

        p <- draw(hm, heatmap_legend_side = "bottom",
        padding = unit(c(2, 120, 2, 2), "mm"))

        pdf(paste0(
            outdir, "/continous.variables/",
            flag, "/Heatmap_anotated_", flag,
            "_", factor, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()
    }
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


surv.model <- function(mod, time, event, flag, outdir) {
  
    if (dir.exists(paste0(outdir, "/survival/", flag)) == F) {
        dir.create(paste0(outdir, "/survival/", flag))
    }


    SurvObject <- Surv(
        time = mod@samples_metadata[, time],
        event = mod@samples_metadata[, event]
    )

    #obtain factors data
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    #Cox model using factors
    fit <- coxph(SurvObject ~ Z)

    s <- summary(fit)
    coef <- s[["coefficients"]]

    #r¡create results table
    sur.df <- data.frame(
        factor = factor(rownames(coef), levels = rev(rownames(coef))),
        p      = coef[, "Pr(>|z|)"],
        coef   = coef[, "exp(coef)"],
        lower  = s[["conf.int"]][, "lower .95"],
        higher = s[["conf.int"]][, "upper .95"]
    )

    write.table(sur.df,
        sep = "\t",
        paste0(outdir, "/survival/", flag, "/Hazard_ratio_factors.txt"),
        row.names = F,
        col.names = T
    )

    #Forest plot
    ggplot(sur.df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
        geom_pointrange(col = "#619CFF", size = 1.3) +
        coord_flip() +
        scale_x_discrete() +
    labs(y = "Hazard Ratio", x = "") +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        theme_bw() + 
        theme(text = element_text(size = 26))

    ggsave(
        filename = paste0(outdir, "/survival/", flag,
        "/Hazard_ratio_factors.pdf"),
        device = "pdf",
        dpi = 900, width = 10, height = 10
    )

}


#function:important.features
#input:
#      mod: MOFA model
#      gmt: gmt.file with pathways
#      nfeatures:number of top features to show in plots
#      th: threshold in importance
#      met_anot: anotation matrix of the CpGs information
#output: Cox model using factor as variables and forest plot

#description: Function to obtain the cox models using factors as features



important.features <- function(mod, outdir, gmt,
nfeatures = 20, th = 0.6, met_anot) {

    genes_l <- list()


    for (v in names(mod@dimensions$D)) {
        # Obtain weights and scale it within each omic and Fator

        w <- get_weights(mod)[[v]]
        w <- as.data.frame(apply(w, 2, function(x) {
            x / max(abs(x))
        }))
        wi <- as.data.frame(w)

        if (v == "Expression") {
            wi <- Anot_path(wi,gmt)
        }

        if (v == "Methylation") {
            wi$id <- gsub("_.*", "", rownames(wi))
            wi <- merge(wi,met_anot, by.x = "id", by.y = "ID")


            head(met_anot)
            
            wi <- wi %>% separate_rows(genes, sep = "_")

            wi <- Anot_path(wi, gmt)
        }

            write.table(wi,
                paste0(outdir, "/features/", v, "_Importance_features.txt"),
                sep = "\t",
                col.names = T,
                row.names = F
            )

            for (f in colnames(w)) {

                nfactor <- as.numeric(gsub("Factor", "", f))
                plot_top_weights(mod,
                    view = v,
                    factors = nfactor,
                    nfeatures = nfeatures, # Top number of features to highlight
                    scale = T # Scale weights from -1 to 1
                ) +
                theme(text = element_text(size = 18))

            if (dir.exists(paste0(outdir, "/features/", f)) == F) {

                dir.create(paste0(outdir, "/features/", f))

                }

                ggsave(
                    filename = paste0(outdir, "/features/", f, "/", v,
                    "_top_features.pdf"),
                    device = "pdf", dpi = 600, width = 10, height = 6
                )

                if (v == "Expression") {
                    
                    #select features by threshold
                    genes_l[[f]][["Expression Positive"]] <- 
                    as.character(unique(wi$genes[wi[, f] >= th]))
                    genes_l[[f]][["Expression Negative"]] <- 
                    as.character(unique(wi$genes[wi[, f] <= -th]))
                    
                    write.table(wi[order(abs(wi[, f]), decreasing = T),
                    c("genes", f, "Pathways")],
                        paste0(outdir, "/features/", f,"/", v,
                        "_top_genes.txt"),
                        sep = "\t",
                        row.names = F,
                        col.names = T
                    )
                }

                if (v == "Methylation") {
                    
                    #select features by threshold

                    genes_l[[f]][["Methylation Positive"]] <-
                    as.character(unique(wi$genes[wi[, f] >= th]))
                    genes_l[[f]][["Methylation Negative"]] <-
                    as.character(unique(wi$genes[wi[, f] <= -th]))

                    write.table(wi[order(abs(wi[, f]), decreasing = T),
                    c("id", f, "Pathways",
                    colnames(met_anot)[2:ncol(met_anot)])],
                        paste0(outdir, "/features/", f,"/", v,
                        "_top_genes.txt"),
                        sep = "\t",
                        row.names = F,
                        col.names = T
                    )
                }
            }
        }

        # Resume table
        tabla_final <- MOFA_list2DF(genes_l, gmt)
        write.table(tabla_final,
            sep = "\t",
            paste0(outdir, "/features/Resume_table_features"),
            row.names = T, col.names = NA
        )
    }


#function:gsea_factors
#input:
#      omic: type of omic expression or methylation
#      mod: MOFA model
#      gmts: gmt.file with pathways
#      nfeatures:number of top features to show in plots
#      th: threshold in importance
#      outdir: Path to write outputs
#      type: select expression or methylation
#      met_anot: If type is mehthylation put anbotyation matrix of CpGs
#      flag: Name to add in files
#output: data.frame with the results of gsea

#description: Function to obtain the results of the gsea

gsea_factors <- function(omic,
mod,
gmts,
outdir,
type,
met_anot = NULL,
flag) {


    
    if (!type %in% c("expression", "methylation")){
    stop("type should be expression or methylation")
    }

    if (type == "methylation"){

        if (is.null(met_anot)){
    stop("provide methylation anotation")
    }
    }


    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample,]
    Z <- as.data.frame(Z)
    Z <- Z[colnames(omic),]

    for (f in colnames(Z)){
        

        df_corr <- Corr_factor(omic,Z,f)
        
        
        if (type == "methylation"){

            df_corr <- merge(df_corr, met_anot, by.x = "row.names",
            by.y = "ID")

            df_corr <- df_corr %>% separate_rows(genes, sep = "_")
            df_corr <- df_corr[!duplicated(df_corr$genes),]
            df_corr <- as.data.frame(df_corr)
            rownames(df_corr) <- df_corr$genes 

        }
        
        df_corr$p.adj <- p.adjust(df_corr$p.value, method = "fdr")
        df_corr <- df_corr[order(df_corr[,'p.value'], decreasing = F),]
        
        
        
        write.table(df_corr
                    ,sep = '\t'
                    , paste0(outdir, "/features/", f, "/",
                    "Correlation_",flag,".txt")
                    , row.names = T
                    ,col.names = NA)

        rnk <- df_corr[order(df_corr[,'cor'], decreasing = T),'cor']
        names(rnk) <- rownames(df_corr[order(df_corr[,'cor'], decreasing = T),])
        rm(df_corr)
        for (gmt in gmts) {


        name <- gsub(".*/|.gmt", "", gmt)

        gseahk <- gsea(rnk = rnk,
                       gmt = gmt
                       , title = paste0(outdir, "/features/", f, "/", name,
                    "_dotplot_",flag,".pdf"))


        gseahk$Factor <- f

        write.table(gseahk
                    , paste0(outdir, "/features/", f, "/", name,
                    "_GSEA_",flag,".txt")
                    ,sep = '\t'
                    ,row.names = F
                    ,col.names = T)



        }

}

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
paste0(outdir,  "/continous.variables/",
            flag, "/Linear_model_p.value.txt"),
sep = "\t",
row.names = T,
col.names = NA)

write.table(r2,
paste0(outdir,  "/continous.variables/",
            flag,"/Linear_model_r2.txt"),
sep = "\t",
row.names = T,
col.names = NA)
write.table(estimate,
paste0(outdir,  "/continous.variables/",
            flag,"/Linear_model_estimate.txt"),
sep = "\t",
row.names = T,
col.names = NA)

p.val[p.val > 0.05] <- 1
p.val <- -log(p.val)
p.val <- p.val[,-1]

col <- list("R2" = colorRamp2(c(1,0),
            c("red", "white")
        ))


anot_legend_param <-  list(
        "R2" = list(title_gp = gpar(fontsize = 40),
        labels_gp = gpar(fontsize = 40),
        labels_gp = gpar(fontsize = 40),
        legend_height = unit(8, "cm"),
        grid_width = unit(2, "cm"),
        title_position= "leftcenter-rot")
        )


        ha_column <- HeatmapAnnotation(
            df = as.data.frame(r2),
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = 40),
            border = unit(40, "cm"),
            simple_anno_size = unit(2, "cm"),
            gp = gpar(fontsize = 40),
            annotation_legend_param = anot_legend_param)


hm <- Heatmap(as.matrix(t(p.val)),
            top_annotation = ha_column,
            show_row_names = T,
            show_column_names = T,
            col = colorRamp2(
                c(max(p.val), 0),
                c("red", "white")
            ),
            name = "Heatmap",
            gap = unit(5, "mm"),
            cluster_rows = F,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "spearman",
            heatmap_legend_param = list(
                title = "-log(p.value)",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 40),
                labels_gp = gpar(fontsize = 40,
                family = "Times", face = "bold"),
                legend_width = unit(20, "cm"),
                grid_height = unit(2, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 36),
            column_names_gp = gpar(fontsize = 36),
            rect_gp = gpar(col = "black", lwd = 2)
        )

        p <- draw(hm, heatmap_legend_side = "bottom",
        padding = unit(c(2, 120, 2, 2), "mm"))

        pdf(paste0(
            outdir,  "/continous.variables/",
            flag, "/immune.linear.model.pdf"
        ),
        width = 20, height = 20)
        print(p)
        dev.off()


}
