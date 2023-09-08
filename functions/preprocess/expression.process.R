#function: norm_exp
#input:
#      targetfile: path were the targetfile is located. This must be a tsv file
#                  with one column called FileName and in each 
#                  row the location of each .CEL file.
#          inpath: Path where there are the anotation files "sample_anotation.txt"
#                  and Anotacion_Clariom_Human_S.txt
#         outpath: Path where the output will be writted
#          norm:   Normalization method either rma , vsn or none
#          flag:   Name to add in the final files. 
#          
#output: expression matrix with genes as rows and samples in 
         #columns which is writted in the ouptut path.

#description: function that generates the normalized expression matrix
#             and write it in the outputh path


norm_exp <- function(targetfile, inpath, outputpath, norm = "rma", flag) {
    
    
    #Generate the raw matrix from .CEL files
    Target <- readTargets(targetfile, row.names = "FileName")
    data <- ReadAffy(filenames = Target$FileName, cdfname = "clariomshumancdf")


    #Normalization
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


    #Armonization of samples code
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

    # Gene Anotation
    anotation <- read.table(paste0(inpath, "/Anotacion_Clariom_Human_S.txt"),
        sep = "\t",
        header = T
    )
    data_rma_tab <- merge(anotation, data_rma_tab,
        by.x = "Probe",
        by.y = "row.names"
    )

    # Select probes with highest variance when they select the same gene
    data_rma_tab <- data_rma_tab[, -1]
    data_rma_tab <- data_rma_tab[complete.cases(data_rma_tab$gene), ]
    num <- apply(data_rma_tab[, 2:ncol(data_rma_tab)], c(1, 2), as.numeric)
    vari <- apply(num, 1, var)
    data_rma_tab <- data_rma_tab[order(vari), ]
    data_rma_tab <- data_rma_tab[!duplicated(data_rma_tab$gene), ]


    rownames(data_rma_tab) <- data_rma_tab$gene

    data_rma_tab <- data_rma_tab[, -1]

    #write data
    write.table(data_rma_tab,
        paste0(outputpath, "/", flag, "_Expression_data_genes.tsv"),
        sep = "\t",
        quote = F,
        col.names = NA,
        row.names = T
    )
}