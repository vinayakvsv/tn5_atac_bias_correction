# Introduction

We ran the following script -- we used a version of the Rscript in which we hardcoded the input of models (this script is available in our `correct_bias` folder

```
Rscript ./seqbias_atac_orchestra_setmodels.R \
/n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/hg38.fa \
SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.pos.sorted.bam \
SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.neg.sorted.bam \
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.bed \
9 \
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo
# the cliplen was set to 4. 
```

# Output
```
# files containing footprint stats before footprint
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection_neg.kmer_freq_before_correction.pdf # k-mer frequencies before correction on the negative strand
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection_pos.kmer_freq_before_correction.pdf # k-mer frequencies before correction on the negative strand
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.ppm # a table containing base transposition frequencies over each position of the aggregated intervals
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.sb.ppm.neg.before.Rlist # Rlist containing information of .ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.sb.ppm.pos.before.Robject # Robject containing information of .ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.seqbias.inputs.neg.Rlist # Inputs in an Rlist for the negative-strand correction
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.seqbias.inputs.pos.Rlist # Inputs in an Rlist for the positive-strand correction
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.before_correction.shifted_freq.txt # contains the Tn5 transposition frequencies over each position of the aggregated intervals. The "before" footprint is computed from this file
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimobefore.footprint.pdf # PDF containing the footprint computed

# results of applying the Buenrostro model 
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.after_buenrostromodel.shifted_freq.txt
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimoafter.footprint.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro_neg.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro_pos.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.sb.ppm.neg.after.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.sb.ppm.pos.after.Robject
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.seqbias.corrected.neg.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.seqbias.corrected.pos.Rlist

# results of applying the Shendure model
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.after_shendurechr22model.shifted_freq.txt
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22_neg.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22_pos.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.sb.ppm.neg.after.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.sb.ppm.pos.after.Robject
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.seqbias.corrected.neg.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.seqbias.corrected.pos.Rlist

# depicts side-by-side plots of the footprint before and after model application
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.footprint.pdf

# depicts the footprint after applying the native-ATAC-seq model, along with spline models fitted to the transposition pattern for occupancy profiling
../MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.pdf

```
