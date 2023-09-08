#Consruct build environment

conda env create -f envs/build.yml -n build

#build estimate environment
conda activate build
conda mambabuild envs/estimate/conda-recipe --output-folder envs/estimate/estimate_package
conda env create -f envs/estimate.yml -n estimate
conda activate estimate
conda install --use-local --update-deps envs/estimate/estimate_package/linux-64/estimate-1.0.13-r35_6.tar.bz2

#build a MOFA.nlcor environment

conda activate build
conda mambabuild envs/nlcor/conda-recipe --output-folder envs/nlcor/nlcor_package
conda env create -f envs/clustering.ensemble.yml -n clustering.ensemble
conda activate MOFA.nlcor
conda install --use-local --update-deps envs/nlcor/nlcor_package/linux-64/r-nlcor-2.03-r40_6.tar.bz2

#build a clustering environment
conda activate build
conda mambabuild envs/mixedClust/conda-recipe --output-folder envs/mixedClust/mixedClust_package
conda env create -f envs/clustering.ensemble.yml -n clustering.ensemble
conda activate clustering.ensemble
conda install --use-local --update-deps envs/mixedClust/mixedClust_package/linux-64/mixedclust-1.0.2-r42_6.tar.bz2


#build a figpatch environment
conda activate build
conda mambabuild envs/figpatch/conda-recipe --output-folder envs/figpatch/figpatch_package
conda env create -f envs/paneles.yml -n clustering.ensemble
conda activate panales
conda install --use-local --update-deps envs/figpatch/figpatch_package/linux-64/r-figpatch-0.2-r42_6.tar.bz2


#build the other environments
conda activate build
conda env create -f envs/arrays_m.yml -n arrays_m
conda env create -f envs/champ.yml -n champ
conda env create -f envs/MOFA.yml -n MOFA
conda env create -f envs/preprocess.yml -n preprocess
conda env create -f envs/viper.yml -n viper
conda env create -f envs/MOFA.compare.yml -n MOFA.compare
conda env create -f -f envs/immune.yml -n immune
conda env create -f -f envs/CNV.yml -n CNV
conda env create -f -f envs/geneset.yml -n geneset
conda env create -f -f envs/survival.yml -n survival
