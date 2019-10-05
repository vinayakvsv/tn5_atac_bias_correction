#!/bin/bash
#
#
#This script will take in a set of folders containing BAM files. This script will then feed into a core script that create bedgraph files out of eacn BAM, merge the bedgraphs together, and create a bigwig out of the result
#Moreover, the core script will generate "stranded" bedgraphs, in which only the positive-strand or negative strand reads in the input BAM files will be put through the script
#The bigwig file that results from this script will provide an aggregate of many samples' coverage of the 5'ends of their read over the genome
#This is suited for hg38
#
##0. Setup
#module load seq/BEDtools/2.23.0
#module load seq/UCSC-tools
#hg38sizes="/n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/chroms/hg38.23chrom.sizes"
#
#1. Take in the input, a name for your job and a list of BAM files
core_script="/n/data2/mgh/ragon/pillai/pipelines/epigenetics_pipeline/motif_footprints/atacseq_footprinting/transpos_site_cov/calculate_5primecov_sereads_tn5correct.sh"
all=("$@")
name=${all[0]}
posCorr=${all[1]}
negCorr=${all[2]}
species=${all[3]}
bamdir=${all[4]}
#
#2. Loop through the list of BAM files, converting each BAM file into a bedgraph. Save up these bedgraph file names
bedgraphs=() #this is not entirely necessary
bedgraphs_pos=()
bedgraphs_neg=()
for i in $(echo $bamdir/S*); do
    wd=$(basename $i)
	mkdir $wd #create a working directory
	bamfiles=$(echo $i/*.bam)
	cd $wd
    bsub -q mcore -n 6 -W 24:00 -o out.txt -e err.txt \
    bash $core_script \
    $name \
    $posCorr \
    $negCorr \
    $species \
    $bamfiles
    cd ..
done
