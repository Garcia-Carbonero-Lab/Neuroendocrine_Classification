# Neuroendocrine_Classification
This repository has all the code used to develop our lab's cross-site classification in well-differentiated neuroendocrine tumors. Unfortunately, the data is not public as these results have not yet been published.

## Preprocess

The first step is to preprocess all data used in this project. The order in which scripts are described is the order to run them

### preprocess folder

**expression.matrix.R**: script to read .cel files and generate the normalized expression matrix with genes as rows and samples as columns.\
**methylation.matrix.R**: script to read .idat files and generate the normalized methylation matrix with probes as rows and samples as columns.\
**estimate.calculation.R**: script to calculate ESTIMATE scores and add this information to clinical data.frame where samples are rows and features are columns.\
**CIMP.calculation.R**: script to calculate the CIMP using methylation data and add this data to the clinical data.frame.\
**immune.infiltration.R**: script to calculate the scores of immune infiltration using MCPcounter R packages.\
**umap.pca.discrete.data.R**: script to create umap and PCA plots showing the distribution of samples by discrete variables.\
**quality.correction.R**: script to remove variance in data related to PCA1 with high correlation with quality.\
**feature.selection.R**: script to select variables for build MOFA model.

### CNV folder

**cnv.champ.R**: script to obtain segments of cnvs using methylation data.\
**gistic.sh**: script to obtain the CNVs peaks from segments using GISTIC2 algorithm.\
**preprocess.cnv.R**: script to process gistic result to obtain a matrix with cnv in rows sample in columns and indicate if there is an amplification, deletion or is wild type.

### genesets folder

**filter.genesets**: script to select genesets to evaluate in by gsea.\ 
**anotate.genesets**: script to change gene names in genesets by the aliases annotated in expression data.

### master_regulators folder:

**aracne.R**: script to obtain the network using ARACNE algorithm.\
**viper.R**: script to obtain master regulator scores using VIPER algorithm.

## MOFA analysis

Once all data is ready to be analyzed, the second step is to build the MOFA model and interrogate each factor to see which clinical and biological features are associated to each factor.
### MOFA folder

**MOFA.creation.R**:script to create the MOFA model.\ 
**MOFA.analysis.R**: script to associate MOFA factors with clinical and biological variables.\
**MOFA.non.linear.analysis.R**: script to calculate the no linear association between factors and continuous features.

## Create a neuroendocrine molecular classification

The third step is to use an ensemble strategy to obtain the consensus classification 

### clustering folder

**clustering**: script to obtain the neuroendocrine subtypes.\
**heatmap.clustering**: script to generate the heatmap showing factors and the annotation of molecular subtypes

## Characterization of molecular subtypes
Finally, we have characterized each of the molecular subtypes from a clinical and biological perspective.

### Cluster.analysis folder

**clinical**: script to investigate associations between clinical and molecular data.\
**survival**: script to study the prognostic value of the neuroendocrine molecular classification.\
**surv.eval**: script to evaluate the prediction capacity of Cox models with and without the molecular classification.\
**immune**: script to study differences in the immune infiltration scores among clusters.\
**cnv**: script to check differences in CNVs among subtypes.\
**deg**: script to apply differential expression analysis.\
**deg**: script to apply differential methylation analysis.\
**corr_met_exp**: script to calculate the correlation between methylation probes and genes annotated in the same DNA region.\
**prognostic.expression**: script to select gene whose expression has a prognostic value using Lasso Cox models.\


## Figuras

In the folder figuras.tesis there is the code used to obtain the figures that appear in the manuscript.

