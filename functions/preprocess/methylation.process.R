



#function: qual.met
#input:
#      import: list obtained from champ.import function

#output: data.frame with sample names as rows and the
#        the -log of the mean of the p value as the only column
#        called qual.met

#description: Calculation of the methylation quality by
#             using the -log of the mean of the p value


qual.met <- function(import) {



    pval <- import$detP
    #calculate the -log of the mean of the p value
    qual <- as.data.frame(-log(apply(pval, 2, mean)))
    colnames(qual) <- "qual.met"
    return(qual)
}







#function: addCIMPstatus
#input:
#      values: matrix with methylation values with
#              probes in rows and samples in columns

#output: data.frame with sample names as rows and the 
#        the -log of the mean of the p value as the only column
#        called qual.met

#description: Calculation of the methylation quality by
#             using the -log of the mean of the p value

addCIMPstatus <- function(values){
    
    #Calculate median of the total methylation
    median_score <- median(values)
    
    #Classify samples into high and 
    #low by comparing they median with the total median

    f <- function(x){
        m <- median(x)
        if (m >= median_score){l <- 'High' }
        if (m <= median_score){l <- 'Low' }
        return(l)
    }
    
    CIMP <- apply(values,2,f)
    CIMP <- as.data.frame(CIMP)
    colnames(CIMP) <- "CIMP.status"
    return(CIMP)
}

