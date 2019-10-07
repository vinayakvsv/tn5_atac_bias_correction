# works for R/3.5.1
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install rtracklayer, Rsamtools, seqLogo, and seqbias
BiocManager::install(c("rtracklayer","Rsamtools","seqLogo","seqbias"))

# install ggplot2
install.packages("ggplot")