correlation.nonlinear.test <- function(data,
Z,
columns,
flag,
outdir,
p.th = 0.05) {


    # We calculate the anova
    ano_df <- lapply(columns, function(x) {

        cdata <- data[complete.cases(data[, x]),]

        vf <- lapply(colnames(Z), function(y) {

            corr <- nlcor(cdata[, x],cdata[, y], plt = F)
            f <- corr$cor.estimate
            p <- corr$adjusted.p.value
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

    p.df[, col] <- p.adjust(p.df[, col], method = "fdr")

    }

    ano_df[grep("p.value", rownames(ano_df)), ] <- p.df

     write.table(ano_df, paste0(
        outdir,
        "/continous.variables/", flag,
        "/non_linear_correlation_Features_Factors_", flag, ".txt"
    ),
    sep = "\t",
    row.names = T,
    col.names = NA
    )

     # we obtain only the p.value
    c.df <- ano_df[grep("estimate", rownames(ano_df)), ]
    # We will show only color result in those that are significant
    # We only show variables with any significance result


    c.df[p.df > p.th] <- 0

   breaksList <- seq(from = 0 , to =1, by = 0.01)

   col <- colorRampPalette(c("lightgrey", "red"))(length(breaksList))

    rownames(c.df) <- gsub(" estimate", "", rownames(c.df))

    pheatmap::pheatmap(t(c.df),
        cluster_rows = F, color = col, breaks = breaksList,
        filename = paste0(
            outdir,
            "/continous.variables/", flag, "/Non_linear_heatmap_Correlation_", flag, ".pdf"
        )
    )

    return(c.df)
}


correlation.non.linear.plot <- function(mod, column, factor,
flag, outdir, estimate) {

 cli.df <- mod@samples_metadata

    cli.df <- cli.df[complete.cases(cli.df[, column]), ]


plt1 <- ggplot(cli.df, aes_string(x = column, y = factor)) +
    theme_classic() +
    geom_point(size = 3) +
    stat_smooth() +
    theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
    )

    ggsave(
        plot = plt1, filename = paste0(
            outdir, "/continous.variables/",
            flag, "/", column, "/Factors_non_linear_corr_plot_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 17,
        height = 12
    )
  

    plt2 <- ggplot(cli.df, aes_string(x = column, y = factor)) +
    theme_classic() +
stat_smooth()  +
geom_point(size = 3) +
    facet_wrap(as.formula("~ PRIMARY_TUMOR"),
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
            flag, "/", column, "/Factors_non_linear_Corr_plot_wrapper_", flag,
            "_", column, "_", factor, ".pdf"
        ),
        device = "pdf",
        dpi = 500, width = 22,
        height = 12
    )

}


continous.features.nonlinear <- function(mod,
p.th = 0.05,
columns,
outdir,
flag,
select = T) {

if (dir.exists(paste0(outdir, "/continous.variables/", flag)) == F) {

    dir.create(paste0(outdir, "/continous.variables/", flag))

}
    Z <- Reduce(rbind, mod@expectations$Z)
    Z <- Z[mod@samples_metadata$sample, ]

    df <- mod@samples_metadata

    c.df <- correlation.nonlinear.test(
        data = df, Z = Z, columns = columns,
        flag = flag, outdir = outdir, p.th = p.th
    )

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
        # Plot each one separately

        for (factor in factors) {
        
        if (select == T){
        if (c.df[column,factor] > 0){

        #Plot significant variables
            correlation.non.linear.plot(
                mod = mod, column = column, factor = factor,
                flag = flag, outdir = outdir
            )

        }} else{

         correlation.non.linear.plot(
                mod = mod, column = column, factor = factor,
                flag = flag, outdir = outdir
            )
        }   
        

        }

    }
}
