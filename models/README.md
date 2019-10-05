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

Please see `correct_bias/` (one folder up) for instructions on how to implement the models in this folder here. 

# Instructions to build your own models

Although our ATAC-seq models should help correct bias at transcription factor binding sites for all types of ATAC-seq experiments, you may choose to build a different model. If so, we have provided instructions.

## Inputs

* **BAM file**
* **Reference FASTA**

##
