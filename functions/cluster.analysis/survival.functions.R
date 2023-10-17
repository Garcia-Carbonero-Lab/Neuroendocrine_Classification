#function:Mofa.model
#input:
#      base: data.frame with samples as rows and clinical information in columns
#       class: column name of feature of interest
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
    
    #Fit cox model
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

    #generate table with results
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
    
    #Forest plot
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
height=8
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
        risk.table = T,
        xlab = "Time",
        ylab = "Overall Survival",
        palette = custcolor,
        font.main = c(18, "bold"),
        font.x = c(18, "bold"),
        font.y = c(18, "bold"),
        font.caption = c(18, "bold"), 
        font.legend = c(18, "bold"), 
        font.tickslab = c(18, "bold"),
        pval.size = 10,
        fontsize = 7) 



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


#function:cox.lasso 
#input:
#      base: data.frame with samples as rows and clinical information in columns
#      data: matrix with features to add the model as rows and samples as columns 
#      OS: overall survival time from dfiagnosis to last interview
#      OS_event : Overall survival status in 0 alive 1 dead
#      covariates: vector of colum names of base to add in cox models
#      outdir: Path to write results
#      flag: Name to add in the name of files
#      alpha: alpha of penalization 1 is Lasso
#      n: number of models to obtain the lambda
#      forced: T forced the model to introduce the covariates. F if lasso cox model
#      selects also the covariates

#output: features selected by cox lasso model and the regression coeficient
#description:function to obtain the features selected by lasso cox models

cox.lasso <- function(base,
data,
OS,
OS_event,
covariates = NULL,
outdir,
flag,
alpha = 1,
n = 1000,
forced = F) {

#complete cases
base <- base[complete.cases(base[,OS]),]
base <- base[complete.cases(base[,OS_event]),]

#zscore data

data <- zscore.rows2(data)

#base and data in the same order
data <- data[,rownames(base)]
data <- as.data.frame(t(data))



if (is.null(covariates) == F){
    base <- base[complete.cases(base[,covariates]),]
    data <- data[rownames(base),]
    genes <- colnames(data)
    data <- cbind(data, base[, covariates])

}

#create the y variable with time and status
y <- base[,c(OS, OS_event)]
colnames(y) <- c("time", "status")


#eliminate - in colnames
fits <- lapply(1:n, function(x){
if (is.null(covariates) == F){
if(forced == T){
fit <- cv.glmnet( x = as.matrix(data),
y = Surv(base[,OS], base[,OS_event] == 1),
family = "cox",
alpha = alpha,
maxit = 1000,
penalty.factor=c(rep(1, length(genes)),
rep(0, ncol(data) - length(genes))
))
}else{

fit <- cv.glmnet( x = as.matrix(data),
y = Surv(base[,OS], base[,OS_event] == 1),
family = "cox",
alpha = alpha,
maxit = 1000)

}
}else{
fit <- cv.glmnet( x = as.matrix(data),
y = Surv(base[,OS], base[,OS_event] == 1),
family = "cox",
alpha = alpha,
maxit = 1000)

}
return(list("lamda.min" = fit$lambda.min, "fit" = fit))
})

lambda.min <- unlist(lapply(fits,function(x){
return(x$lamda.min)}))

fit <- lapply(fits,function(x){
return(x$fit)})

tb <- table(lambda.min)
lambda <- names(tb)[tb == max(tb)]
fit <- fit[lambda.min == lambda][[1]]

pdf(paste0(outdir, "/CV.cox.plot_",flag,".pdf"))
plot(fit)
dev.off()
#obtenin coef != 0 using lambda min
coef.df <- coef(fit, s = "lambda.min")
coef <- coef.df[coef.df[,1] != 0,]

if(length(coef) == 1){
names(coef) <- rownames(coef.df)[coef.df[,1] != 0]
}
return(coef)
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