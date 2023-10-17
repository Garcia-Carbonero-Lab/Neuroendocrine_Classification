#This script need the environment surveval
library(survival)  
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
set.seed(123)

source("functions/cluster.analysis/surv.eval.R")

# we read config file
config <- read.csv("config/config.tsv",
sep = "\t",
header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
wkdir <- values["wkdir"]
datadir <- values["datadir"]


if (dir.exists(paste0(wkdir, "/clustering/survival.eval")) == F) {

    dir.create(paste0(wkdir, "/clustering/survival.eval"))

}


outdir <- paste0(wkdir, "/clustering/survival.eval")



base <- read.csv(paste0(datadir, "/clinical_data.txt"),
                 sep = "\t",
                 row.names = 1,
                 header = T)


base$Subtype <- as.character(base$Subtype)
base$PRIMARY_TUMOR <- factor(base$PRIMARY_TUMOR,
levels = c("GASTRIC","COLORECTAL", "LUNG", "PANCREAS", "SMALL INTESTINE"))



# We asses that the data is numeric  
base$OS.time <- as.numeric(base$OS.time)
base$OS.event <- as.numeric(base$EXITUS)
#We generate the rurvival objectbase
base$Surv_OS <- with(base, Surv(base$OS.time, base$OS.event == 1))

features <- c("GRADE", "AGE", "PRIMARY_TUMOR", "GENDER","STAGE_IV")

cindexm1 <- list()
cindexm2 <- list()

# Generate models 100 times
for (i in 1:100){

# Random samples
index_data = sample(1:nrow(base), nrow(base)*0.8)

# create Cox models with and without subtypes
m1 = coxph( as.formula(paste('Surv_OS ~ ',paste(features, collapse = ' + '))), data = base[index_data,], x = T)
m2 = coxph( as.formula(paste('Surv_OS ~ Subtype +',paste(features, collapse = ' + '))), data = base[index_data,], x = T)

#obtain the cindex
cindex <- concordance(m1,m2)

cindexm1 <- append(cindexm1,cindex[[1]][1])
cindexm2 <- append(cindexm2,cindex[[1]][2])
}


# Compare the C index between models with Subtype information and without subtype information
df <- data.frame("Cox.model" = unlist(cindexm1), "Cox.model.NS" = unlist(cindexm2))

df <- df %>% pivot_longer(c(Cox.model, Cox.model.NS),names_to = "Model", values_to = "C.index")

stat.test <- df  %>% 
  t_test(C.index ~ Model, paired = T) %>%
  add_significance() %>% add_xy_position(x = "Model")

stat.test[["Model"]] <- NA



plt <- ggplot( df , aes(y = .data[["C.index"]],
    x = .data[["Model"]], fill = .data[["Model"]])) 
        

        plt <- plt + geom_boxplot(
            width = 0.6, outlier.shape = NA
        ) + 
        geom_point(position = position_jitterdodge()) 
         
plt <- plt + scale_fill_manual(values = list("Cox.model" = "darkslategray3",
    "Cox.model.NS" = "slateblue3")) +

        theme_classic() +
        stat_pvalue_manual(stat.test,
        label = "p.signif",  hide.ns = T,
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

    plt$layers[[3]]$aes_params$size <- 1
    plt$layers[[3]]$aes_params$label.size <- 20
    plt$layers[[3]]$aes_params$vjust <- 0.5


    ggsave(
        plot = plt, filename = paste0(
            outdir, "/Cindex.pdf"
        ),
        device = "pdf",
        dpi = 500, width = 20,
        height = 15
    )



write.table(stat.test[,c(1:9)],
paste0(outdir,"/T_test_results,txt"),
sep = "\t",
row.names = F,
col.names = T)





