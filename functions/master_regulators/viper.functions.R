
#function:RegProcess
#input:
#      a.file: aracne network file
#      exp.matrix: expression matrix with genes in rows and samples in columns
#      out.dir: Output path to write results
#      out.name: name of output file

#output: ggplot object

#description: function to change text size and legend size in ggplot object


RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  require(viper)
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}
