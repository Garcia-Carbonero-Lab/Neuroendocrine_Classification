
#function:heatmap
#input:
#      data: matrix with features as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      outdir: Path to write results
#      col: list indicating colors of anotation
#      size: vector indicating width in first term and height in secon term
#      anot_param_legend: list indicating parameetrs of anotation
#      flag: string to add in name of the file
#      max: top number for color NULL select the maximum in the matrix
#      min: minimum number for the color NULL select the minimum number in the matrix
#      names.size: size of rownames
#      column: column to split samples
#      anot.h: anotation height
#      annotation_name_gp: size of anotation names
#      show_row_names: True rownames is shown False row names are not shown
#      color.max: color name for top numbers
#      color.min: color name of minimum numbers

#output: heatmap

#description: function to create a heatmap

heatmap <- function(data,
base,
outdir, col, size,
anot_param_legend,
flag,
max = NULL,
min = NULL,
names.size = 3,
column,
anot_h,
annotation_name_gp = 60,
show_row_names = T,
color.max = "red",
color.min = "blue") {

data <- data[, rownames(base)]
data <- zscore.rows2(data)

anot <- base[, names(col)]

ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = annotation_name_gp),
            border = unit(40, "cm"),
            annotation_height = anot_h,
            height = unit(10, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            annotation_legend_param = anot_param_legend)

if (is.null(max)){

    max <- max(data)
}


if (is.null(min)){

    min <- min(data)
}
col <- colorRamp2(
                c(max, 0, min),
                c(color.max, "white", color.min)
)

hm1 <- Heatmap(as.matrix(data),
            top_annotation = ha_column,
            show_row_names = show_row_names,
            show_column_names = F,
            col = col,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            column_dend_side = "top",
            column_dend_height = unit(6, "cm"),
            cluster_columns = T,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "Z score",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            column_split = anot[, column],
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE
        )



hlist <- hm1

p <- draw(hlist, heatmap_legend_side = "bottom",
        annotation_legend_side = "top",
        padding = unit(c(20, 240, 20, 20), "mm"))

pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()

}





#function:wilcox_plots
#input:
#      data: matrix with features as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      group.x: name of the colum for x axis in boxplot
#      group.fill: Column name to be the fill in boxplot
#      color.fill: vector with color of levels in cgroup.fill
#      outdir: Path to write the plots
#      flag: string to add in names of the files


#output: boxplots with wilcoxon p value and a table with ther esults of wilcoxon 
# test

#description: function to compare continous variables using wilcoxon rank test


wilcox_plots <- function(data,
base,
group.x,
group.fill,
color.fill,
outdir,
flag) {

data <- data[rownames(base),,drop = F]

df<- cbind(data,base[,c(group.x,group.fill)])
df[,group.x] <- as.factor(as.character(df[,group.x]))
df[,group.fill] <- as.factor(as.character(df[,group.fill]))

stats <- list()

for (row in colnames(data)) {

if (group.x != group.fill){
stat.test <- df %>%
  group_by(.dots=group.x) %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

   median <- df[,c(row,group.x, group.fill)] %>% 
    group_by_at(vars(one_of(c(group.x,group.fill)))) %>%
    get_summary_stats(type = "median")
}else{
stat.test <- df %>%
  wilcox_test(as.formula(paste0(row, "~", group.fill))) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

  median <- df[,c(row, group.fill)] %>% 
    group_by_at(vars(group.fill)) %>%
    get_summary_stats(type = "median") 

}


median_l <- lapply(1:nrow(stat.test), function(x){
dfmedian <- data.frame(median_group_1 = 
median$median[median$Subtype == stat.test[[x,"group1"]]],
median_group_2 = 
median$median[median$Subtype == stat.test[[x,"group2"]]]
)
})

median <- do.call(rbind,median_l)

stats[[row]] <- as.data.frame(cbind(stat.test,median))

stat.test <- stat.test %>%
  add_xy_position(x = group.x, dodge = 0.8)

bxp <- ggboxplot(
  df, x = group.x, y = row,
  fill = group.fill,
  width = 0.6 ,
  palette = color.fill
  ) +


  stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0, hide.ns = T
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +

theme(
        text = element_text(size = 26),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 26)
        ) +

        labs(x = "", y = row)

    bxp$layers[[2]]$aes_params$size <- 1
    bxp$layers[[2]]$aes_params$label.size <- 20

ggsave(
        plot = bxp, filename = paste0(
            outdir, "/Boxplot.primary.subtype_", row, ".pdf"
        ),
        device = "pdf",
        dpi = 300, width = 16,
        height = 12
    )


}

stats <- do.call(rbind,stats)

write.table(stats,
paste0(outdir, "/Wilcox_", flag, ".txt"),
sep = "\t",
row.names = F,
col.names = T)

}


#function:heatmap.molecular
#input:
#      exp: matrix with genes as rows and samples as columns
#      met: matrix with methylation probes as rows and samples as columns
#      base: data.frame with samples as rows and features as columns
#      outdir: Path to write results
#      col: list indicating colors of anotation
#      size: vector indicating width in first term and height in secon term
#      anot_param_legend: list indicating parameetrs of anotation
#      flag: string to add in name of the file
#      max.exp: top number for color NULL select the maximum in the expression matrix
#      min.exp: minimum number for the color NULL select the minimum number in the expression matrix
#      max.met: top number for color NULL select the maximum in the matrix
#      min.met: minimum number for the color NULL select the minimum number in the matrix
#      names.size: size of rownames
#      column: column to split samples
#      anot.h: anotation height
#      annotation_name_gp: size of anotation names
#      show_row_names: True rownames is shown False row names are not shown
#      scale_met: True if scale methylation data
#      scale_exp: True if scale methylation data
#      color.max: color name for top numbers
#      color.min: color name of minimum numbers
#      height_anotation: height of anotation data
#      cluster: exp or met to select which omic use to applay clustering
#      order: order of the dendogram

#output: heatmap

#description: function to create a heatmap with expression and methylation data
heatmap.molecular <- function(exp,
met,
base,
outdir, col, size,
anot_param_legend,
flag,
max.exp = NULL,
min.exp = NULL,
max.met = NULL,
min.met = NULL,
names.size = 3,
column,
anot_h,
annotation_name_gp = 60,
show_row_names = T,
scale_met = T,
scale_exp = T,
height_anotation = 25,
cluster = "exp",
order) {

exp <- exp[, rownames(base)]


base[,column] <- as.factor(as.character(base[,column]))


if (scale_exp == T){
exp <- zscore.rows2(exp)
}



met <- met[, rownames(base)]

if (scale_met == T){
met <- zscore.rows2(met)
}

if (cluster == "exp"){
cluster_cols <- cluster_within_group(exp,base[,column])

}

if (cluster == "met"){
cluster_cols <- cluster_within_group(met,base[,column])
}

#order data as dendogram
base <- base[labels(cluster_cols),]
base <- base[c(grep(order[[1]],base[,column]),
grep(order[[2]],base[,column]),
grep(order[[3]],base[,column])),]
exp <- exp[, rownames(base)]
met <- met[, rownames(base)]
anot <- base[, names(col)]

ha_column <- HeatmapAnnotation(
            df = anot,
            show_legend = T,
            col = col,
            annotation_name_gp = gpar(fontsize = annotation_name_gp),
            border = unit(40, "cm"),
            annotation_height = anot_h,
            height = unit(height_anotation, "cm"),
            gp = gpar(fontsize = 80),
            simple_anno_size_adjust = TRUE,
            annotation_legend_param = anot_param_legend,
            simple_anno_size = unit(4, "cm"))

if (is.null(max.exp)){

    max.exp <- max(exp)
}
if (is.null(max.met)){

    max.met <- max(met)
}

if (is.null(min.exp)){

    min.exp <- min(exp)
}

if (is.null(min.met)){

    min.met <- min(met)
}

col.exp <- colorRamp2(
                c(max.exp, 0, min.exp),
                c("red", "white", "blue")
)



col.met <- colorRamp2(
                c(max.met, 0, min.met),
                c("green", "white", "blue")
)

hm1 <- Heatmap(as.matrix(exp),
            top_annotation = ha_column,
            show_row_names = show_row_names,
            show_column_names = T,
            col = col.exp,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "expression Z score",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
#            column_split = anot[, column],
            height = unit(40, "cm")

        )

hm2 <- Heatmap(as.matrix(met),
            show_row_names = show_row_names,
            show_column_names = F,
            col = col.met,
            gap = unit(5, "mm"),
            cluster_rows = T,
            show_row_dend = F,
            show_column_dend = F,
            cluster_columns = F,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            heatmap_legend_param = list(
                title = "methylation Z score",
                title_position = "topcenter",
                legend_direction = "horizontal",
                title_gp = gpar(fontsize = 60),
                labels_gp = gpar(fontsize = 60,
                family = "Times", face = "bold"),
                legend_width = unit(30, "cm"),
                grid_height = unit(2.5, "cm")
            ),
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 100),
            row_names_max_width = unit(names.size,'cm'),
            cluster_column_slices = FALSE,
            height = unit(40, "cm")
        )

hlist <- hm1 %v% hm2

p <- draw(hlist, heatmap_legend_side = "top",
        annotation_legend_side = "top",
        padding = unit(c(20, 240, 20, 20), "mm"),
        ht_gap = unit(3, "cm"))

pdf(paste0(
            outdir, "/Heatmap_anotated_",flag, ".pdf"
        ),
        width = size[["width"]], size[["height"]]
        )
        print(p)
        dev.off()

}
