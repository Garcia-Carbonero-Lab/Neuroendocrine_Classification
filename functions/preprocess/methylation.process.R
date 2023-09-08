


qual.met <- function(import) {

    # methylation quality is calculated
    # using the -log of the mean of
    # the p value

    pval <- import$detP
    qual <- as.data.frame(-log(apply(pval, 2, mean)))
    colnames(qual) <- "qual.met"
    return(qual)
}

addCIMPstatus <- function(beta){
    
    median_score <- median(beta)
    
    f <- function(x){
        m <- median(x)
        if (m >= median_score){l <- 'High' }
        if (m <= median_score){l <- 'Low' }
        return(l)
    }
    
    CIMP <- apply(beta,2,f)
    CIMP <- as.data.frame(CIMP)
    colnames(CIMP) <- "CIMP.status"
    return(CIMP)
}

scale.beta <- function(beta){

      (beta - min(beta)) / (max(beta) - min(beta))

}