# Using Seqbias models computed from ATAC-seq data

A Seqbias model (explained in Jones et al 2012) represents the dependencies in some read count features between positions in a set of sequences. The model is appropriate for representing positional biases in read counts -- for instance, if two genomic bases in an ATAC-seq experiment tend to be transposed together whenever Tn5 is in the general region, then the two sites will tend to yield similar transposition frequencies across all cells studied.
One way this model can be used is to determine whether certain bases around a given read's start site tend to also contain a large number of read start sites -- suggesting that the bases tend to be selected for fragmentation and ultimately for generating sequencing reads in an experiment. Ideally, in a genomic sequencing experiment, we would like all bases to have an equal chance of being captured with a read start site. However, if certain positions around tend to be represented more than others, we can observe a small bias that can be problematic whenever we analyze features at the genomic distance scale as the bases around reads (~10^1 bases). Transcription factor binding sites are within this order of magnitude.

Each seqbias model is meant to capture the biases in transposition frequency within given stretches of sequences around the entire genome. Seqbias uses this model to adjust the genome-wide transposition frequencies seen within the given sequences to reflect what they would ideally look like without bias. 

## Our computed models

We have presented two types of models, both of which should capture similar bias profiles for Tn5
1. A model trained on Tn5 transposition patterns upon naked genomic DNA; the data used to train this model came from Adey et al 2010. Since this genomic DNA is stripped of all native chromatin marks, it should present the bias of Tn5 transposition without the influence of intervening chromatin (although the bias profiles for native chromatin and naked DNA should be very similar). This data does not contain any information that could identify the person from whom the DNA came from.
2. A model trained on Tn5 transposition patterns upon native chromatin in GM12878 cells; the data used to train this model came from Buentrostro et al 2013 (i.e. the original ATAC-seq paper). In theory, any ATAC-seq dataset should be able to largely reproduce the bias profile trained in this model, since (as we argue in our paper) Tn5 bias in inherent. 

For each type of model, we variants that capture bias over different distances from ATAC-seq read start sites.  We also trained separate models based on the strand (one positive-strand vs one negative-strand). The bias is the same on either strand, but since Tn5 leaves an overhang where it cuts, the insertion sites are slightly offset on each Tn5 (4 bases of 5' overhang).  

You would only need to use one of these models (Adey or Buenrostro) and the model with the desired size profile. We would recommend applying each strand-specific model to a count matrix generated from reads aligning to each strand separately and centered on your TF of interest. You can then combine the two strand-specific corrected profiles into a single unified transposition track to examine the final corrected footprint. 

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


# Instructions to build your own models

Although our ATAC-seq models should help correct bias at transcription factor binding sites for all types of ATAC-seq experiments, you may choose to build a different model. If so, we have provided instructions.

## Inputs

* **BAM file**
* **Reference FASTA**

##

