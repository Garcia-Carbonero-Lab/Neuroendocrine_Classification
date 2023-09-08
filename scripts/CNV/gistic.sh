#!/bin/sh
## run example GISTIC analysis

## output directory
echo --- creating output directory ---
datadir=$(cat config/config.tsv | grep datadir |cut -f2)
wkdir=$(cat config/config.tsv | grep wkdir |cut -f2)

basedir=$wkdir"/CNV/GISTIC"
mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
segfile=$wkdir"/CNV/CNV_segments/CNA_methylation.segments"

refgenefile=$datadir"/methylome/anotation/hg19.UCSC.add_miR.140312.refgene.mat"

## call script that sets MCR environment and calls GISTIC executable 

export LANGUAGE
export LC_ALL="C"
gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1  -conf 0.99 -armpeel 1 -savegene 1 -gcm mean -savedata 1 --maxseg 10000
