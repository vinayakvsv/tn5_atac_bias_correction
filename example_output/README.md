# Introduction

We ran the following script -- we used a version of the Rscript in which we hardcoded the input of models (this script is available in our `correct_bias` folder

```
Rscript ./seqbias_atac_orchestra_setmodels.R \
/n/data2/mgh/ragon/pillai/referenceGenomes/hg38_ucsc/hg38.fa \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/bamfiles/SRR891268_buenrostro/SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.pos.sorted.bam \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/bamfiles/SRR891268_buenrostro/SRR891268.atac-seq.bowtie2.q10filtered.sorted.hg38-autosom-sex-chr.sorted.neg.sorted.bam \
/home/vvv3/pillaifolder/analysis/tn5_cutsite/buenrostro-et-al-2013/seqbias_correct_footprint/intervals/motifs_intersect_macs2peaks/intersect_SRR891268.narrow/JASPAR_core_vertebrate_nonoredundnant/bedfiles/BED6/MA00/MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.bed \
9 \
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo
# the cliplen was set to 4. 
```

# Output
```
# files containing footprint stats before footprint
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection_neg.kmer_freq_before_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection_pos.kmer_freq_before_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.sb.ppm.neg.before.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.sb.ppm.pos.before.Robject
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.seqbias.inputs.neg.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_beforecorrection.seqbias.inputs.pos.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.before_correction.shifted_freq.txt
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimobefore.footprint.pdf

# 
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.after_buenrostromodel.shifted_freq.txt
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimoafter.footprint.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.after_shendurechr22model.shifted_freq.txt
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo.footprint.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro_neg.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro_pos.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.sb.ppm.neg.after.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.sb.ppm.pos.after.Robject
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.seqbias.corrected.neg.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postbuenrostro.seqbias.corrected.pos.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22_neg.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22_pos.kmer_freq_after_correction.pdf
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.ppm
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.sb.ppm.neg.after.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.sb.ppm.pos.after.Robject
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.seqbias.corrected.neg.Rlist
MA0098.3_ETS1.motif.meme.fimo_out.BED6.fimo_postshendurechr22.seqbias.corrected.pos.Rlist
```
