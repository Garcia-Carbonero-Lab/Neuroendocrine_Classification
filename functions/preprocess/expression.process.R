Norm_exp <- function(targetfile, inpath, outputpath, norm = "rma", flag) {
    Target <- readTargets(targetfile, row.names = "FileName")
    data <- ReadAffy(filenames = Target$FileName, cdfname = "clariomshumancdf")

    if (norm == "vsn") {
        ndata <- vsnrma(data)
    }

    if (norm == "rma") {
        ndata <- expresso(data,
            bg.correct = TRUE,
            bgcorrect.method = "rma",
            normalize = TRUE,
            normalize.method = "quantiles",
            pmcorrect.method = "pmonly",
            summary.method = "medianpolish",
            verbose = TRUE
        )
    }


    if (norm == "none") {
        ndata <- expresso(data,
            bg.correct = TRUE,
            bgcorrect.method = "rma",
            normalize = FALSE,
            pmcorrect.method = "pmonly",
            summary.method = "medianpolish",
            verbose = TRUE
        )
    }


    data_rma_tab <- as.data.frame(exprs(ndata))
    colnames(data_rma_tab) <- gsub(
        "_\\(Clariom_S_Human\\).CEL", "",
        colnames(data_rma_tab)
    )
    colnames(data_rma_tab)[grep("NENs", colnames(data_rma_tab))] <- gsub(
        "^[0-9]*\\.", "",
        colnames(data_rma_tab)[grep("NENs", colnames(data_rma_tab))]
    )
    sample_anot <- read.table(paste0(inpath, "/sample_anotation.txt"),
        sep = "\t",
        header = T
    )


    data_rma_tab <- t(data_rma_tab)



    data_rma_tab <- merge(sample_anot, data_rma_tab,
        by.x = "Muestras", by.y = "row.names"
    )
    rownames(data_rma_tab) <- data_rma_tab$Código.paciente
    data_rma_tab <- data_rma_tab[, 3:ncol(data_rma_tab)]
    data_rma_tab <- as.data.frame(t(data_rma_tab))

    # Anotamos con la info de clariomshuman
    anotation <- read.table(paste0(inpath, "/Anotacion_Clariom_Human_S.txt"),
        sep = "\t",
        header = T
    )
    data_rma_tab <- merge(anotation, data_rma_tab,
        by.x = "Probe",
        by.y = "row.names"
    )

    # Colapsamos por varianza
    data_rma_tab <- data_rma_tab[, -1]
    data_rma_tab <- data_rma_tab[complete.cases(data_rma_tab$gene), ]
    num <- apply(data_rma_tab[, 2:ncol(data_rma_tab)], c(1, 2), as.numeric)
    vari <- apply(num, 1, var)
    data_rma_tab <- data_rma_tab[order(vari), ]
    data_rma_tab <- data_rma_tab[!duplicated(data_rma_tab$gene), ]


    rownames(data_rma_tab) <- data_rma_tab$gene

    data_rma_tab <- data_rma_tab[, -1]
    write.table(data_rma_tab,
        paste0(outputpath, "/", flag, "_Expression_data_genes.tsv"),
        sep = "\t",
        quote = F,
        col.names = NA,
        row.names = T
    )
}