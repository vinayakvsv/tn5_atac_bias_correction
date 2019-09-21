# tn5_atac_bias_correction

This repository contains the code, models, and procedure implemented in Viswanadham et al 2019 "A Bayesian approach for correcting Tn5 transposition bias in ATAC-seq footprinting" (preprint here: https://www.biorxiv.org/content/10.1101/525808v1). 

As a test case, we will replicate the example in the paper in which we re-analyze data from Buenrostro et al 2013 (i.e. the original ATAC-seq paper) and estimate true genomewide footprints from artefactual footprints

This repository contains the following:
1. R-code that implements the Seqbias package from Jones et al 2012 (Bioconductor: https://bioconductor.org/packages/release/bioc/html/seqbias.html; Paper: https://academic.oup.com/bioinformatics/article/28/7/921/209263; Github repository: https://github.com/dcjones/seqbias) 
2. Auxilliary functions to compute the necessary input from the BAM files
3. YAML files containing the models that we constructed from ATAC-seq data and from Adey et al 2010 (Paper: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119). 
4. R-code to compute our footprint score estimates and to compute the degree of correction from before-and-after comparisons of transposition frequencies

This Github page will be updated as more and more TF footprint estimates are computed. 
