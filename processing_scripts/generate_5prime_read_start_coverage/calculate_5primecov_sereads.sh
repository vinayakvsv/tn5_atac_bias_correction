#!/bin/bash
#
#
#This script combines elements of the following scripts:
#/n/data2/mgh/ragon/pillai/analysis/bcell_subsets/gc/gc.atac-seq/Run_195/human.samples/alignments.q10filtered.sorted.hg38-autosom-sex-chr/5prime.cov.bdgs/test_5primecov_calculate/test_vsm_script/5primecovBedgraphs_merge.sh
#/n/data2/mgh/ragon/pillai/analysis/bcell_subsets/gc/gc.atac-seq/Run_195/human.samples/alignments.q10filtered.sorted.hg38-autosom-sex-chr/5prime.cov.bdgs/test_5primecov_calculate/test_vsm_script/5primecov_calculate.stranded.vsm.sh
#This script will take in a list of BAM files, create bedgraph files out of them, merge the bedgraphs together, and create a bigwig out of the result
#Moreover, the script will generate "stranded" bedgraphs, in which only the positive-strand or negative strand reads in the input BAM files will be put through the script
#The bigwig file that results from this script will provide an aggregate of many samples' coverage of the 5'ends of their read over the genome
#This is suited for hg38
#
#0. Setup
module load seq/BEDtools/2.23.0
module load seq/UCSC-tools
hg38sizes="/n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/chroms/hg38.23chrom.sizes"
#
#1. Take in the input, a name for your job and a list of BAM files
all=("$@")
name=${all[0]}
bamfiles=${all[*]:1}
#
#2. Loop through the list of BAM files, converting each BAM file into a bedgraph. Save up these bedgraph file names
bedgraphs=() #this is not entirely necessary
bedgraphs_pos=()
bedgraphs_neg=()
for i in $bamfiles; do
    bamfile=$(basename $i) #get the bam file with its extension and no path
    bamname=${bamfile%.bam} #get the name of the bam file
    bedtools genomecov -ibam $i -5 -bg | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.23chrs.bedGraph
    bedgraphs+=($bamname.5primecov.23chrs.bedGraph) #Add the name of the bigwig with positive- and negative-stranded reads to the list
    bedtools genomecov -ibam $i -5 -bg -strand + | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.23chrs.pos.bedGraph
    bedgraphs_pos+=($bamname.5primecov.23chrs.pos.bedGraph) #Add the name of the bigwig with positive-stranded reads to the list
    bedtools genomecov -ibam $i -5 -bg -strand - | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.23chrs.neg.bedGraph
    bedgraphs_neg+=($bamname.5primecov.23chrs.neg.bedGraph) #Add the name of the bigwig with negative-stranded reads to the list
    #NB 1: the grep -v "_" command is intended to remove nucleotides on chromosomes not part of the primary assembly. We could just remove all instances of "_random" (mitochondrial and unmapped regions were removed earlier in alignment quality control)
    bedGraphToBigWig $bamname.5primecov.23chrs.bedGraph $hg38sizes $bamname.5primecov.23chrs.bw
    bedGraphToBigWig $bamname.5primecov.23chrs.pos.bedGraph $hg38sizes $bamname.5primecov.23chrs.pos.bw
    bedGraphToBigWig $bamname.5primecov.23chrs.neg.bedGraph $hg38sizes $bamname.5primecov.23chrs.neg.bw
done
#
#3. Sum up the bedgraphs in each of the lists
bedtools unionbedg -i ${bedgraphs[*]} | awk -v OFS='\t' '{print $1,$2,$3,$4+$5+$6}' > $name.merged.5primecov.23chrs.bedGraph
bedtools unionbedg -i ${bedgraphs_pos[*]} | awk -v OFS='\t' '{print $1,$2,$3,$4+$5+$6}' > $name.merged.5primecov.23chrs.pos.bedGraph
bedtools unionbedg -i ${bedgraphs_neg[*]} | awk -v OFS='\t' '{print $1,$2,$3,$4+$5+$6}' > $name.merged.5primecov.23chrs.neg.bedGraph
#
#4. Convert each of the summed-up bedgraphs into a bigwig file
bedGraphToBigWig $name.merged.5primecov.23chrs.bedGraph /n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/chroms/hg38.23chrom.sizes $name.merged.5primecov.23chrs.bw
bedGraphToBigWig $name.merged.5primecov.23chrs.pos.bedGraph /n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/chroms/hg38.23chrom.sizes $name.merged.5primecov.23chrs.pos.bw
bedGraphToBigWig $name.merged.5primecov.23chrs.neg.bedGraph /n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/chroms/hg38.23chrom.sizes $name.merged.5primecov.23chrs.neg.bw
#