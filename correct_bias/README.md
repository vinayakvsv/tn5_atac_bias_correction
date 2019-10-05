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

We have provided a single command to implement a model of interest. Note that we require two models, one for the positive-stranded reads and one for the negative stranded reads. Given our observation that the Tn5 bias is identical for positive and negative strands, you can input the same model for those options

```
Rscript ./seqbias_atac_orchestra.R \
<reference genome FASTA> \
<positive-strand BAM> \
<negative-strand BAM> \
<BED file with motif instances> \
<length of motif> \
<name of output> \
<seqbias model that you wish to apply for positive strand> \
<seqbias model that you wish to apply for negative strand>
```

This script will automatically compute transposition frequencies for 100 bases on either side of your motif instances.  


