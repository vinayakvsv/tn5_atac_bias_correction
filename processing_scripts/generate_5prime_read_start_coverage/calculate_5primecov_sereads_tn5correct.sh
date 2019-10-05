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
#
#1. Take in the input, a name for your job and a list of BAM files
all=("$@")
name=${all[0]}
posCorr=${all[1]}
negCorr=${all[2]}
species=${all[3]}
bamfiles=${all[*]:4}
gsizes="/n/data2/mgh/ragon/pillai/referenceGenomes/"$species"_ucsc/chroms/"$species".chrom.sizes"
#
echo "5' end correction for positive-strand reads: $posCorr bases upstream" >> comments.txt
echo "5' end correction for negative-strand reads: $negCorr bases upstream" >> comments.txt
echo "Genome sizes file: $gsizes" >> comments.txt
#
#2. Loop through the list of BAM files, converting each BAM file into a bedgraph. Save up these bedgraph file names
bedgraphs=() #this is not entirely necessary
bedgraphs_pos=()
bedgraphs_neg=()
for i in $bamfiles; do
    bamfile=$(basename $i) #get the bam file with its extension and no path
    bamname=${bamfile%.bam} #get the name of the bam file
    #bedtools genomecov -ibam $i -5 -bg | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.23chrs.bedGraph
    #bedgraphs+=($bamname.5primecov.23chrs.bedGraph) #Add the name of the bigwig with positive- and negative-stranded reads to the list
    #positive strand
    bedtools genomecov -ibam $i -5 -bg -strand + | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.pos.bedGraph
    bedtools slop -l -$posCorr -r $posCorr -g $gsizes -i $bamname.5primecov.pos.bedGraph > $bamname.5primecov.pos.tn5adjust.bedGraph #shift the positions of the positive-strand Tn5 cut sites downstream by $posCorr bases. Since "slop" will expand intervals in both directions, the positive-strand has a negative $posCorr for its "start coordinate"
    bedgraphs_tn5adjust_pos+=($bamname.5primecov.pos.tn5adjust.bedGraph) #Add the name of the bigwig with positive-stranded reads to the list #This may not be necessary
    #negative strand
    bedtools genomecov -ibam $i -5 -bg -strand - | grep -v "_" | uniq | sort -k1,1 -k2,2n > $bamname.5primecov.neg.bedGraph
    bedtools slop -l $negCorr -r -$negCorr -g $gsizes -i $bamname.5primecov.neg.bedGraph > $bamname.5primecov.neg.tn5adjust.bedGraph #shift the positions of the negative-strand Tn5 cut sites upstream by 5 bases. Since "slop" will expand intervals in both directions, the negative-strand has a negative $posCorr for its "end coordinate"
    bedgraphs_tn5adjust_neg+=($bamname.5primecov.neg.tn5adjust.bedGraph) #Add the name of the bigwig with negative-stranded reads to the list
    #NB 1: the grep -v "_" command is intended to remove nucleotides on chromosomes not part of the primary assembly. We could just remove all instances of "_random" (mitochondrial and unmapped regions were removed earlier in alignment quality control)
    #Make bigwigs for each sample fed into the script
    #bedGraphToBigWig $bamname.5primecov.bedGraph $hg38sizes $bamname.5primecov.bw
    #bedGraphToBigWig $bamname.5primecov.pos.bedGraph $hg38sizes $bamname.5primecov.pos.bw
    #bedGraphToBigWig $bamname.5primecov.neg.bedGraph $hg38sizes $bamname.5primecov.23chrs.neg.bw
    rm $bamname.5primecov.pos.bedGraph
    rm $bamname.5primecov.neg.bedGraph
done
#
#3. Sum up the bedgraphs in each of the lists
##bedtools unionbedg -i ${bedgraphs[*]} | awk -v OFS='\t' '{print $1,$2,$3,$4+$5+$6}' > $name.merged.5primecov.bedGraph #We will need to modify the awk command to add all entries after the third column
##awk '{ for(i=4; i<=NF;i++) j+=$i; print j; j=0 }' <input bedgraph> to sum up the bedgraph scores in each row
##bedtools unionbedg -i ${bedgraphs_tn5adjust_pos[*]} | awk -v OFS='\t' '{print $1,$2,$3,$4+$5+$6}' > $name.merged.5primecov.pos.tn5adjust.bedGraph #We will need to modify the awk command to add all entries after the third column
#
#4. Merge the strand-adjusted bedgraphs and sort the result
#bedtools unionbedg -i $name.merged.5primecov.pos.tn5adjust.bedGraph $name.merged.5primecov.neg.tn5adjust.bedGraph | awk -v OFS='\t' '{print $1,$2,$3,$4+$5}' > $name.merged.5primecov.tn5adjust.bedGraph
bedtools unionbedg -i *.5primecov.pos.tn5adjust.bedGraph *.5primecov.neg.tn5adjust.bedGraph | awk '{ for(i=4; i<=NF;i++) j+=$i; print $1,$2,$3,j; j=0 }' > $name.5primecov.tn5adjust.merged.bedGraph
#
#5. Convert each of the summed-up bedgraphs into a bigwig file
bedGraphToBigWig $name.5primecov.tn5adjust.merged.bedGraph $gsizes $name.5primecov.tn5adjust.merged.bw
#
#include operations to delete the Uncorrected sample bedgraphs

# Only do this if there are multiple samples for a strand
if [ ${#bedgraphs_tn5adjust_pos[@]} -gt 1 ]; then
    bedtools unionbedg -i *.5primecov.pos.tn5adjust.bedGraph | awk '{ for(i=4; i<=NF;i++) j+=$i; print $1,$2,$3,j; j=0 }' > $name.5primecov.pos.tn5adjust.merged.bedGraph #We will need to modify the awk command to add all entries after the third column
    bedGraphToBigWig $name.5primecov.pos.tn5adjust.merged.bedGraph $gsizes $name.5primecov.pos.tn5adjust.merged.bw
fi
#
if [ ${#bedgraphs_tn5adjust_neg[@]} -gt 1 ]; then
    bedtools unionbedg -i *.5primecov.neg.tn5adjust.bedGraph | awk '{ for(i=4; i<=NF;i++) j+=$i; print $1,$2,$3,j; j=0 }' > $name.5primecov.neg.tn5adjust.merged.bedGraph #We will need to modify the awk command to add all entries after the third column
    bedGraphToBigWig $name.5primecov.neg.tn5adjust.merged.bedGraph $gsizes $name.5primecov.neg.tn5adjust.merged.bw
fi
#