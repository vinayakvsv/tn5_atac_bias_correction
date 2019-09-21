# Using Seqbias models computed from ATAC-seq data

A Seqbias model (explained in Jones et al 2012) represents the dependencies in some read count features between positions in a set of sequences. The model is appropriate for representing positional biases in read counts -- for instance, if two genomic bases in an ATAC-seq experiment tend to be transposed together whenever Tn5 is in the general region, then the two sites will tend to yield similar transposition frequencies across all cells studied.
One way this model can be used is to determine whether certain bases around a given read's start site tend to also contain a large number of read start sites -- suggesting that the bases tend to be selected for fragmentation and ultimately for generating sequencing reads in an experiment. Ideally, in a genomic sequencing experiment, we would like all bases to have an equal chance of being captured with a read start site. However, if certain positions around tend to be represented more than others, we can observe a small bias that can be problematic whenever we analyze features at the genomic distance scale as the bases around reads (~10^1 bases). Transcription factor binding sites are within this order of magnitude.

Each seqbias model is meant to capture the biases in transposition frequency within given stretches of sequences around the entire genome. Seqbias uses this model to adjust the genome-wide transposition frequencies seen within the given sequences to reflect what they would ideally look like without bias. 

## Our computed models


# Instructions
