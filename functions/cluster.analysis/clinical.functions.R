#function: clinical.table
#input:
#      data: data.frame with samples as rows and features as columns
#      features: column names to analyze
#      group: column name of the feature to generate groups
#      na : True if na samples are taken into account. False if Na samples are 
#      removed

#output: data.frame indicating the percent of each level pf features in each label of group
# and the p value of fisher test

#description: function to evaluate the distribution of features in the group levels


clinical.table <- function(data, features, group, na = F) {

    dfl <- lapply(features, function(feature) {
        if (na == F){
        count <- table(data[, feature], data[, group])
        }else{
        
        count <- table(data[, feature], data[, group],useNA = "a")

        }
        f <- try(fisher.test(count, simulate.p.value = T))

        if (class(f) == 'try-error'){
            p <- 'Not calculated'
        }else{
            p <- round(f$p.value,5)
        }
        if (na == F){
        total <- table(data[, feature])
        }else{
        
        total <- table(data[, feature], useNA = "a")

        }
        count <- cbind(count,total)
        pc <- round(t(t(count) / apply(count, 2, sum)) * 100, 2)
        
        df <- do.call(rbind,
                      lapply(seq_len(nrow(count)),
                             function(x) {
                                 paste0(count[x, ], " (", pc[x, ], "%)")
                             }))
        
        df <- cbind(df, c(p,rep('',nrow(df) - 1)))
        
        colnames(df) <- c(paste0(colnames(count), " (", c(table(data[, group]),
        sum(nrow(data))), ")")
                          ,'P.value')
        rownames(df) <- paste(feature, rownames(count), sep = ".")
        
        return(df)
    })
    
    table <- as.data.frame(do.call(rbind, dfl))
    
    return(table)

}


#function: kruskal
#input:
#      base: data.frame with samples as rows and features as columns
#      features: column names to analyze
#      group: column name of the feature to generate groups
#      outdir: path where write the results

#output: matrix indicating the p.value and adjusted p.value of kruskal wallis

#description: function to evalluate the differences among subtypes of the 
#continous features

kruskal <- function(base,
group, 
features,
outdir,
flag){

vf <- lapply(features, function(x) {
            df <- base[complete.cases(base[,x]),]
            
            k <- kruskal.test(as.formula(paste(x, " ~ ", group)),
                data = base)
            p <- k[[3]]

            return(c("p.value" = p))
        })


df <- as.data.frame(do.call(rbind,vf))
rownames(df) <- features
df$p.adjust <- p.adjust(df$p.value, method = "BH")

write.table(df,
paste0(outdir, "/Kruskal_table_", flag, ".txt"),
row.names = T,
col.names = NA,
sep = "\t")

return(df)
}

#function: boxplots
#input:
#      base: data.frame where samples are in rows and features in columns
#      group.x:  variable x in boxplot
#      group.fill: group to separate and generate colors
#      color.fill: vector of colors for group.fill levels
#      features: names of columns to study
#      th: threshold of the p.value to use it to create the boxplot
#      outdir: Path where write the results
#      flag: name to add in files
#      width: width of pdf
#      height: height of pdf

#output: MOFA model

#description: function to create the MOFA model

boxplots <- function(base,
group.x,
group.fill,
color.fill,
features,
th = 0.05,
outdir,
flag,
width = 16,
height = 12){


if (dir.exists(paste0(outdir, "/", flag)) == F) {

    dir.create(paste0(outdir, "/", flag))

}




base <- base[,c(group.x,group.fill,features)]
base <- base[complete.cases(base[,group.x]),]
base <- base[complete.cases(base[,group.fill]),]

base[,group.x] <- as.factor(as.character(base[,group.x]))
base[,group.fill] <- as.factor(as.character(base[,group.fill]))

stat.df <- list()
for (feature in features) {


if (group.x != group.fill){

dfa <- base[,c(feature,group.x,group.fill)]
dfa <- dfa[complete.cases(dfa[,feature]),]
#Wilcox test
stat.test <- dfa %>%
  group_by_at(vars(group.x)) %>%
  wilcox_test(as.formula(paste0(feature, "~", group.fill))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = group.x, dodge = 0.8)

median <- dfa[,c(feature,group.x, group.fill)] %>% 
    group_by_at(vars(one_of(c(group.x,group.fill)))) %>%
    get_summary_stats(type = "median")

we <- dfa %>%
  group_by_at(vars(group.x)) %>%
  wilcox_effsize(as.formula(paste0(feature, "~", group.fill)))


stat.test <- stat.test %>% 
unite(col= "order", c( as.name(group.x), group1:group2), sep = "_",
remove = FALSE) 

we$names <- paste(we$.y., 
    we$group1,
    we$group2,
    we[,group.x][[1]],
    sep = "_")

stat.test$names <- paste(stat.test$.y., 
stat.test$group1,
stat.test$group2,
stat.test[, group.x][[1]],
sep = "_")
    
stat.test <- merge(stat.test, we[,c(4,8:9)], by = "names")

}else{

dfa <- base[,c(feature,group.fill)]
dfa <- dfa[complete.cases(dfa[,feature]),]

stat.test <- dfa %>%
  wilcox_test(as.formula(paste0(feature, "~", group.fill))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = group.x, dodge = 0.8)

median <- dfa[,c(feature, group.fill)] %>% 
    group_by_at(vars(group.fill)) %>%
    get_summary_stats(type = "median") 



stat.test <- stat.test %>% 
unite(col= "order", c(group1:group2), sep = "_",
remove = FALSE) 

we <- dfa %>%
  wilcox_effsize(as.formula(paste0(feature, "~", group.fill)))


we$names <- paste(we$.y., 
    we$group1,
    we$group2,
    sep = "_")

stat.test$names <- paste(stat.test$.y., 
stat.test$group1,
stat.test$group2,
sep = "_")
stat.test <- merge(stat.test, we[,c(4,7:8)], by = "names")

}

stat.df[[feature]] <- stat.test
stat.test[[group.fill]] <- NA


bxp <- ggplot(
  dfa, aes(x = .data[[group.x]],
  y = .data[[feature]], fill= .data[[group.fill]])) +
  
  geom_boxplot(width = 0.6, outlier.shape = NA,
  aes(fill= .data[[group.fill]])) +

  scale_fill_manual(values = as.list(c(color.fill,"blue", "green"))) +

  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.7)) +
  
  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = T
) +


  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme_classic() +

theme(
        text = element_text(size = 40),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 40)
        ) +

        labs(x = "", y = feature)


bxp$layers[[3]]$aes_params$size <- 1.2
bxp$layers[[3]]$aes_params$label.size <- 12 


ggsave(
        plot = bxp, filename = paste0(
            outdir, "/", flag, "/Boxplot.subtype_", feature, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = width,
        height = height
    )

}
stat.df <- stat.df %>% bind_rows()
stat.df <- stat.df[!duplicated(stat.df),]

write.table(
stat.df[,c(1:12,18:19)],
file = paste0(outdir, "/", flag, "/Wilcoxon.table.txt"),
sep = "\t",
row.names = F,
col.names = T)

write.table(
median,
file = paste0(outdir, "/", flag, "/Median.table.txt"),
sep = "\t",
row.names = F,
col.names = T)
}