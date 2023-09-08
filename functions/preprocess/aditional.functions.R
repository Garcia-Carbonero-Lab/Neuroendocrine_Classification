

theme_general <- function(plot) {
    themed <- plot +
        theme(
            text = element_text(size = 15),
            legend.key.height = unit(10, "mm")
        )
}



pca_corr <- function(data, base, columns, outpath, scale = T, flag,
                     anotate1 = c(80, 150),
                     anotate2 = c(80, 140)) {
    if (scale == T) {
        data <- zscore.rows2(data)
    }


    # We calculate PCA and variance
    pcal <- prcomp(t(data), scale. = F)
    pca <- as.data.frame(pcal$x[, 1:10])
    prop_varianza <- pcal$sdev^2 / sum(pcal$sdev^2)


    base_num <- base[, columns]
    base_num <- merge(base_num, pca, by = "row.names")

    corrl <- list()
    est <- list()

    for (c in columns) {
        for (p in colnames(pca)) {

            # We calculate correlation between PCAs and variables
            bs <- base_num[complete.cases(base_num[c]), ]
            corr <- cor.test(bs[, p], bs[, c])
            corrl[[c]][[p]] <- paste0(
                "corr=", round(corr$estimate[[1]], 2),
                "/pval=", round(corr$p.value, 6)
            )
            est[[c]][[p]] <- round(corr$estimate, 2)
        }
        # We plot the PCA1 and PCA2 indicating as color the variables
        plt <- ggplot(aes_string("PC1", "PC2", color = c), data = bs) +
            geom_point() +
            theme_classic() +
            xlab(paste("PCA1", round(prop_varianza[[1]], 2))) +
            ylab(paste("PCA2", round(prop_varianza[[2]], 2))) +
            scale_color_gradient2(
                midpoint = median(base_num[, c], na.rm = T),
                low = "blue", mid = "white", high = "red",
                name = c
            ) +
            annotate("text",
                x = anotate1[1], y = anotate1[2], size = 4,
                label = paste0(
                    "correlación PC1 = ",
                    est[[c]][["PC1"]]
                )
            ) +

            annotate("text",
                x = anotate2[1], y = anotate2[2], size = 4,
                label = paste0(
                    "correlación PC2 = ",
                    est[[c]][["PC2"]]
                )
            ) +

    theme(text=element_text(size=20),
          legend.key.height = unit(1,'cm'),
          aspect.ratio=1/1
          )

        ggsave(plot = plt, filename = paste0(
            outpath, "/", c,
            "_", flag, "_PCA.pdf"
        ), height = 8, width = 8, device = "pdf"
           ,dpi =  300)
    }
    df_corr <- do.call(rbind, corrl)
    return(plt)
}

umap_corr <- function(data, base, columns, outpath, scale = T, flag,
                      anotate1 = c(2, 4),
                      anotate2 = c(2, 3)) {
    if (scale == T) {
        data <- zscore.rows2(data)
    }
    base <- base[colnames(data), ]
    set.seed(1232)
    data.umap <- umap(t(data), preserve.seed = T)
    data.umapl <- data.umap$layout
    data.umapl <- data.frame("UMAP1" = data.umapl[, 1], "UMAP2" = data.umapl[, 2], row.names = rownames(data.umapl))
    data.umapl <- merge(data.umapl, base[, columns], by = "row.names")

    corrl <- list()
    est <- list()
    for (c in columns) {
        for (p in c("UMAP1", "UMAP2")) {
            bs <- data.umapl[complete.cases(data.umapl[c]), ]
            corr <- cor.test(bs[, p], bs[, c])
            corrl[[c]][[p]] <- paste0("corr=", round(corr$estimate[[1]], 2), "/pval=", round(corr$p.value, 6))
            est[[c]][[p]] <- round(corr$estimate, 2)
        }
        plt <- ggplot(aes_string("UMAP1", "UMAP2", color = c), data = bs) +
            geom_point() +
            theme_classic() +
            xlab("UMAP1") +
            ylab("UMAP2") +
            scale_color_gradient2(
                midpoint = median(bs[, c], na.rm = T),
                low = "blue", mid = "white", high = "red", name = "value"
            ) +
            annotate("text",
                x = anotate1[1], y = anotate1[2], size = 4,
                label = paste0(
                    "corr.UMAP1=",
                    est[[c]][["UMAP1"]]
                )
            ) +
            annotate("text",
                x = anotate2[1], y = anotate2[2], size = 4,
                label = paste0(
                    "corr.UMAP2=",
                    est[[c]][["UMAP2"]]
                )
            ) +
            ggtitle(c)


        plt <- theme_general(plt)
        ggsave(plot = plt, filename = paste0(
            outpath, "/", c,
            "_", flag, "_UMAP.pdf"
        ), height = 8, width = 8)
    }
    df_corr <- do.call(rbind, corrl)
    return(plt)
}


generate_pca <- function(data) {
    pca_data <- prcomp(t(data), scale. = F)
    pca_data <- as.data.frame(t(pca_data$x))
    pca_data <- t(pca_data)
    pca_data <- pca_data[, 1, drop = F]
    return(pca_data)
}


pca_corr_discrete <- function(data, base, columns, outpath, scale = T, flag,
                              anotate1 = c(100, 150),
                              anotate2 = c(100, 140)) {
    if (scale == T) {
        data <- zscore.rows2(data)
    }
    # We generate the PCAs 1:10 and we add this to clinical data

    base <- base[colnames(data), ]
    pcal <- prcomp(t(data), scale. = F)
    pca <- as.data.frame(pcal$x[, 1:10])
    prop_varianza <- pcal$sdev^2 / sum(pcal$sdev^2)
    base <- merge(base, pca, by = "row.names")
    corrl <- list()
    est <- list()
    # For ech column
    for (c in columns) {
        df <- base[complete.cases(base[, c]), ]


        # We calculate a linear model with ecah variable
        # and correlate the fittted values with the PCAs

        for (p in colnames(pca)) {
            df <- base[complete.cases(base[, c]), ]
            model <- aov(as.formula(paste(p, " ~ ", c)), data = df)
            corr <- cor.test(df[, p], model$fitted.values)
            corrl[[c]][[p]] <- paste0(
                "corr=", round(corr$estimate[[1]], 2),
                "/pval=", round(corr$p.value, 6)
            )
            est[[c]][[p]] <- round(corr$estimate, 2)
        }


        # If the levels are > 12 que use the palette Paired
        if (length(unique(df[, c])) <= 12) {
            plt <- ggplot(aes_string("PC1", "PC2", color = c), data = df) +
                geom_point() +
                scale_color_brewer(palette = "Paired", name = "label") +
                theme_classic() +
                xlab(paste("PCA1", round(prop_varianza[[1]], 2))) +
                ylab(paste("PCA2", round(prop_varianza[[2]], 2))) +
                annotate("text",
                    x = anotate1[1], y = anotate1[2], size = 4,
                    label = paste0(
                        "corr.PC1=",
                        est[[c]][["PC1"]]
                    )
                ) +
                annotate("text",
                    x = anotate2[1], y = anotate2[2], size = 4,
                    label = paste0(
                        "corr.PC2=",
                        est[[c]][["PC2"]]
                    )
                ) +
                ggtitle(c)


            plt <- theme_general(plt)
        }

        if (length(unique(df[, c])) > 12) {
            qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
            col_vector <- unlist(mapply(
                brewer.pal, qual_col_pals$maxcolors,
                rownames(qual_col_pals)
            ))

            plt <- ggplot(aes_string("PC1", "PC2", color = c), data = df) +
                geom_point() +
                scale_color_manual(values = sample(
                    col_vector,
                    length(unique(df[, c]))
                ), name = "label") +
                theme_classic() +
                xlab(paste("PCA1", round(prop_varianza[[1]], 2))) +
                ylab(paste("PCA2", round(prop_varianza[[2]], 2))) +
                annotate("text",
                    x = anotate1[1], y = anotate1[2], size = 4,
                    label = paste0(
                        "corr.PC1=",
                        est[[c]][["PC1"]]
                    )
                ) +
                annotate("text",
                    x = anotate2[1], y = anotate2[2], size = 4,
                    label = paste0(
                        "corr.PC2=",
                        est[[c]][["PC2"]]
                    )
                ) +
                ggtitle(c)


            plt <- theme_general(plt)
        }



        ggsave(plot = plt, filename = paste0(
            outpath, "/", c,
            "_", flag, "_PCA.pdf"
        ), height = 8, width = 8)

        df <- do.call(rbind, corrl)
        colnames(df) <- colnames(pca)
    }
    return(df)
}


umap_corr_discrete <- function(data, base, columns, outpath, scale = T, flag,
                               anotate1 = c(2, 4),
                               anotate2 = c(2, 3)) {
    library(RColorBrewer)
    if (scale == T) {
        data <- zscore.rows2(data)
    }

    # We generate the UMAp1 and UMAP2 and we add this to clinical data
    base <- base[colnames(data), ]
    set.seed(1232)
    data.umap <- umap(t(data), preserve.seed = T)
    data.umapl <- data.umap$layout
    data.umapl <- data.frame(
        "UMAP1" = data.umapl[, 1],
        "UMAP2" = data.umapl[, 2], row.names = rownames(data.umapl)
    )

    data.umapl <- merge(data.umapl, base[, columns], by = "row.names")

    corrl <- list()
    est <- list()

    for (c in columns) {
        df <- data.umapl[complete.cases(data.umapl[, c]), ]

        for (p in c("UMAP1", "UMAP2")) {
            model <- aov(as.formula(paste(p, " ~ ", c)), data = df)
            corr <- cor.test(df[, p], model$fitted.values)
            corrl[[c]][[p]] <- paste0(
                "corr=", round(corr$estimate[[1]], 2),
                "/pval=", round(corr$p.value, 6)
            )
            est[[c]][[p]] <- round(corr$estimate, 2)
        }


        if (length(unique(df[, c])) <= 12) {
            plt <- ggplot(aes_string("UMAP1", "UMAP2", color = c),
                data = df
            ) +
                geom_point() +
                scale_color_brewer(palette = "Paired", name = "label") +
                theme_classic() +
                annotate("text",
                    x = anotate1[1], y = anotate1[2], size = 4,
                    label = paste0(
                        "corr.UMAP1=",
                        est[[c]][["UMAP1"]]
                    )
                ) +
                annotate("text",
                    x = anotate2[1], y = anotate2[2], size = 4,
                    label = paste0(
                        "corr.UMAP2=",
                        est[[c]][["UMAP2"]]
                    )
                ) +
                ggtitle(c)


            plt <- theme_general(plt)
        }

        if (length(unique(df[, c])) > 12) {
            qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
            col_vector <- unlist(mapply(
                brewer.pal, qual_col_pals$maxcolors,
                rownames(qual_col_pals)
            ))

            plt <- ggplot(aes_string("UMAP1", "UMAP2", color = c), data = df) +
                geom_point() +
                scale_color_manual(values = sample(
                    col_vector,
                    length(unique(df[, c]))
                ), name = "label") +
                theme_classic() +
                annotate("text",
                    x = anotate1[1], y = anotate1[2], size = 4,
                    label = paste0(
                        "corr.UMAP1=",
                        est[[c]][["UMAP1"]]
                    )
                ) +
                annotate("text",
                    x = anotate2[1], y = anotate2[2], size = 4,
                    label = paste0(
                        "corr.UMAP2=",
                        est[[c]][["UMAP2"]]
                    )
                ) +
                ggtitle(c)


            plt <- theme_general(plt)
        }



        ggsave(plot = plt, filename = paste0(
            outpath, "/", c,
            "_", flag, "_UMAP.pdf"
        ), height = 8, width = 8)
    }
    df <- do.call(rbind, corrl)
    colnames(df) <- c("UMAP1", "UMAP2")

    return(df)
}

filter_tissue <- function(data, base, tissues, n = 5000, flag) {
    
    #we applied the selection in each primary tumor separately
    mad_tissue <- lapply(tissues, function(x) {
    #
    patients <- rownames(base)[base$BIOPSY_LOCATION == x &
    
    #we use only primary tumors
    base$BIOPSY_SITES == "PRIMARY"]

        df <- data[, colnames(data) %in% patients]
        v <- apply(df, 1, mad)
        df <- df[order(v, decreasing = T), ]


        #we order genes by MAD
        res <- data.frame(rownames(df), 1:nrow(df))
        colnames(res) <- c("genes", x)
        res <- res[order(res$genes), ]
        return(res)
    })

    #we obtain the genes with highest MAD in all tissues
    mad_tissue <- bind_cols(mad_tissue)
    rownames(mad_tissue) <- mad_tissue$genes...1
    mad_tissue <- mad_tissue[, grep("genes", invert = T, colnames(mad_tissue))]
    mad_tissue$GLOBAL <- apply(mad_tissue, 1, median)
    write.table(mad_tissue, paste0(flag, "_MAD_TISSUE.txt"),
        sep = "\t", row.names = T, col.names = NA
    )
    genes <- rownames(mad_tissue[order(mad_tissue$GLOBAL), ])[1:n]
    data <- data[genes, ]
    return(data)
}