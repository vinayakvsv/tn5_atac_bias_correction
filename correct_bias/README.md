# Introduction

This folder contains the main scripts to impelment our seqbias model


# Instructions to apply existing models

These instructions will help you apply our models as they are. This page will be updated to describe how you can train your own model from your own ATAC-seq data.

## Inputs

* **BAM file** We require that you split your ATAC-seq `.bam` file into two separate BAM files, one with positive strand-aligned reads and one with negative strand-aligned reads.
* **Motif intervals** You will need a `.bed` file with intervals containing all transcription factor binding site instances at which you are interested in observing footprints. For ATAC-seq data, we recommend you search for transcription factor binding site instances in accessible intervals -- i.e. peaks of ATAC-seq read coverage.
* **Reference FASTA**

## Motif intervals
If you wish to identify transcription factor binding site motif instances, you can use a position-weight matrix (PWM) to search the genome (or some desired set of intervals) using `FIMO` from the `MEME` suite (see http://meme-suite.org/doc/fimo.html). We have provided scripts for you to search instances across intervals you wish.

## Command

We have provided a single command to implement a model of interest

```
Rscript /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/seqbias_atac_orchestra_setmodels.R \
<reference genome FASTA> \
<positive-strand BAM> \
<negative-strand BAM> \
<BED file with motif instances> \
<length of motif> \
<number of bases to clip from the positive and negative strands of the intervals> \
<name of output> \
<seqbias model that you wish to apply>
```

This script will automatically compute sequnces

Example:
```
Rscript /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/seqbias_atac_orchestra_setmodels.R \
/n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/hg38.fa \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/bamfiles/SRR891268_buenrostro/SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.pos.sorted.bam \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/bamfiles/SRR891268_buenrostro/SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.neg.sorted.bam \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/intervals/motifs_intersect_macs2peaks/intersect_SRR891268.narrow/JASPAR_core_vertebrate_nonoredundnant/bedfiles/BED6/MA00/MA0027.2_EN1.motif.meme.fimo_out.BED6.fimo.bed \
7 \
MA0027.2_EN1.motif.meme.fimo_out.BED6.fimo \
4
```
