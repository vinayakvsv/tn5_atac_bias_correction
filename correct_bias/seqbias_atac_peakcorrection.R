seqbias_atac <- function(input_intervals,
                         ref_seq,
                         sample_bam,
                         n,m,
                         binarycount=TRUE,
                         k=1,
                         selectchrs,
                         L=50,R=50,
                         indet_strands='+'){
  #make sure the FIMO BED is reformatted as such: chr, start, end, name, score, strand
  #1. Require packages
  require(Rsamtools)
  require(seqbias)
  require(ggplot2)
  require(rtracklayer)
  #2. Set the reference sequence--a FASTA file--and the alignments--a BAM file
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  #3. Open the reference
  ref_f <- FaFile( ref_fn)
  open.FaFile( ref_f )
  print(ref_f)
  #4. Extract the sequences in the reference and take a random sample of n sequences that are m bases long
  ref_seqs <- scanFaIndex( ref_f )
  print(ref_seqs)
  #input_intervals is a BED file with the genomic intervals in which you would like to correct for count biases
  I <- import(con = input_intervals)
  #If the strand is indeterminate (i.e. strand is "*"), set it to the value specified in indet_strands
  print(paste("setting * strands to",indet_strands))
  star_idx <- as.logical( I@strand == '*')
  strand(I[star_idx]) <- indet_strands
  #strand(I) <- rep("+",length(I@seqnames))
  #Scan the sequences of the input intervals
  print("Scanning input intervals...")
  seqs <- scanFa( ref_f, I )
  print("finished scanning the input sequences")
  #Make sure that you get the reverse complementary sequence of the negative-stranded sample sequences
  print("Getting the reverse complement of negative-stranded sequences")
  neg_idx <- as.logical( I@strand == '-')
  seqs[neg_idx] <- reverseComplement( seqs[neg_idx] )
  #print(seqs)
  #5. Count the number of reads within the sample sequences. The counts are binary--that is, as long as one match is found, the count is 1
  print(paste("generating read counts matrix for the sequences in the intervals. Binary counts is",binarycount))
  counts <- count.reads( reads_fn, I, binary = binarycount )
  #6. Get the frequencies of particular k-mers in the sample sequences. By default, k=1
  print("generating k-mer frequencies.")
  freqs <- kmer.freq( seqs, counts ,k=k,L = L,R = R)
  #change L and R to be at most half the size of the maximum length of one of the input sequences
  #Alternatively, make sure that the motifs are slopped
  #7. Plot the frequencies of the k-mers
  print("Plotting k-mer frequencies")
  ylimit <- max(freqs$freq)
  P <- qplot( x = pos,
              y = freq,
              ylim = c(0.0,ylimit+0.01),
              color = seq,
              data = freqs,
              geom = "line" )
  P <- P + facet_grid( seq ~ . )
  print(P)
  output <- list(ref_seqs,I,seqs,counts,freqs)
  names(output) <- c("ref_seqs","I","seqs","counts","freqs")
  return(output)
}
#
seqbias_compensation <- function(ref_seq,
                                 sample_bam,
                                 L,R) {
  #fit the seqbias model. Specify an interval that is L bases to the left of the read start and R bases to the right
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  sb <- seqbias.fit( ref_fn, reads_fn, L = L, R = R )
  return(sb)
}
#
estimate_bias <- function(sb.in,
                          I.in,
                          seqs.in,
                          counts.in,
                          k.in=1,
                          L=50,R=50) {
  #1. Estimate bias
  print("estimate bias")
  bias <- seqbias.predict( sb = sb.in, I = I.in ) #You need to make sure that the "*" strands are set to "+" in the seqbias_atac function
  #2. Adjust the counts
  print("adjust counts")
  counts.adj <- mapply( FUN = '/', counts.in, bias, SIMPLIFY = F )
  #3. Recalculate the k-mer frequencies. Use the same value of k as in seqbias_atac
  print("recalculate kmer frequency")
  freqs.adj <- kmer.freq( seqs = seqs.in, counts = counts.adj, L=L, R=R,k = k.in )
  P <- qplot( x = pos,
              y = freq,
              ylim = c(0.00,0.5),
              color = seq,
              data = freqs.adj,
              geom = "line" )
  P <- P + facet_grid( seq ~ . )
  print(P)
  output <- list(bias,counts.adj,freqs.adj)
  names(output) <- c("bias","counts.adj","freqs.adj")
  return(output)
}
#
kl_divergence <- function(freqs_bkgd,
                          freqs_fgd,
                          name) {
  #' At each position in the intervals studied, the KL divergence of the k-mers found is calculated.
  #' The frequencies of the k-mers at the position represent a probability distribution, and so the KL
  #' divergence is calculated to determine how different the k-mer representation is at each genomic 
  #' position between the background (i.e. before correction) and the foreground (i.e. after correction)
  #
  if (length(freqs_bkgd) != length(freqs_fgd)) {
    print("background and foreground frequencies must be calculated over the same length of bases")
    return(NULL)
  }
  #
  #specify the Kullback-Leibler divergence function
  kl.fn <- function(f.bg,f.obs) { 
    kl_x <- f.obs*log2(f.obs/f.bg) + f.bg*log2(f.bg/f.obs) #calculate the divergence between the two distributions for a particular state in the distribution
    kl_divergence <- sum(kl_x) #sum over the instances
    return(kl_divergence)
  }
  #
  kmers <- unique(freqs_bkgd$seq)
  indices <- seq(1,length(freqs_bkgd$freq),length(kmers))
  kl <- mapply(FUN = function(x,y) {
    inds <- seq(x,x+y-1,1)
    fg <- freqs_fgd$freq[inds]
    print("fg")
    print(fg)
    bg <- freqs_bkgd$freq[inds]
    print("bg")
    print(bg)
    #freqs <- freqs_fgd$freq[inds]/freqs_bkgd$freq[inds]
    kl.fn(f.bg = bg,f.obs = fg)
  }, indices,length(kmers))
  pos <- freqs_bkgd$pos[seq(1,length(freqs_bkgd$pos),length(kmers))]
  kl.df <- cbind(pos,kl)
  pdf(paste(name,"kl_divergence","pdf",sep = "."))
  plot(kl.df,type="l")
  dev.off()
  return(kl.df)
}
#
make_counts_matrix <- function(seqbias_counts,
                               motiflen) {
  #sb.counts <- t(as.data.frame(lapply(X = seqbias_counts,function(x) {unlist(x)})))
  #sb.counts <- t(as.data.frame(lapply(X = seqbias_counts,unlist)))
  #
  #sb.atac.pos.neg <- t(mapply(function(x,y) {c(unlist(x),rev(unlist(y)))},sb.atac.ctcf.pos$counts,sb.atac.ctcf.neg$counts))
  #rbind.sb.atac.ctcf.pos <- do.call(rbind,sb.atac.ctcf.pos$counts)
  #rbind.sb.atac.ctcf.neg <- do.call(rbind,sb.atac.ctcf.neg$counts)
  #rbind.sb.atac.ctcf.pos.neg <- cbind(rbind.sb.atac.ctcf.pos,rbind.sb.atac.ctcf.neg)
  #
  sb.counts <- as.data.frame(t(as.data.frame(lapply(X = seqbias_counts,unlist))))
  sb.counts.sum <- apply(X = sb.counts,MARGIN = 2,sum)
  sb.counts.freq <- sb.counts.sum/sum(sb.counts.sum)
  #plot the footprint
  plot(sb.counts.freq,xaxt="n",type="l",ylab="Cut site probability",xlab="Position from motif")
  lim <- (length(sb.counts.sum) - motiflen)/2
  print(lim)
  axis(1,at=1:length(sb.counts.freq),labels=c(seq(-lim,-1,1),rep(0,motiflen),seq(1,lim,1)))
  abline(v=lim+1,lty=2)
  abline(v=lim+motiflen,lty=2)
  #output
  output <- list(sb.counts,sb.counts.sum,sb.counts.freq)
  names(output) <- c("sb.counts","sb.counts.sum","sb.counts.freq")
  return(output)
}
#
make_pwm <- function(freqs,
                     name) {
  #' Take the frequencies calculated above and transform them into a PWM
  library(seqLogo)
  #1. collapse the nucleotide frequency table that seqbias spits out into a proper PPM (position probability matrix)
  nts <- unique(freqs$seq)
  ppm <- sapply(nts,function(x) {
    nt.freqs <- freqs[freqs$seq == x,"freq"]
  })
  ppm <- as.data.frame(ppm)
  colnames(ppm) <- nts
  rownames(ppm) <- unique(freqs$pos)
  #2. Adjust the nucleotide frequencies such that they add up to 1 in each line
  ppm.corrected <- apply(X = ppm,MARGIN = 1,function(x) {
    x/sum(x)
  })
  freqs.pwm <- makePWM(pwm = ppm.corrected)
  #convert the ppm into TRANSFAC format for WEBLOGO
  freqs.pwm.transfac <- cbind(PO=seq(1,ncol(freqs.pwm@pwm)),t(freqs.pwm@pwm))
  write.table(freqs.pwm.transfac,file=paste(name,"ppm",sep="."),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  #fConn <- file(paste(name,"ppm",sep = "."), 'r+') 
  lines <- readLines(paste(name,"ppm",sep = "."))
  writeLines(c("ID","BF", lines), con = paste(name,"ppm",sep = "."))
  return(freqs.pwm.transfac)
}
#
#####
#MAIN
#1. Take in command-line input and set the source of the libraries
args <- commandArgs(trailingOnly=TRUE)
ref_seq <- args[1] #the reference genome
sample_bam_pos <- args[2] #the positive BAM
sample_bam_neg <- args[3] #the negative BAM
peak_intervals <- args[4] #the motif intervals
name <- args[5] #the name associated with the job
#shendurepos <- args[9]
#shendureneg <- args[10]
#buenrostropos <- args[11]
#buenrostroneg <- args[12]
#
#.libPaths( c( .libPaths(), "/n/data2/mgh/ragon/pillai/nonLSFpkgs/R/") )
.libPaths( c( "/n/data2/mgh/ragon/pillai/nonLSFpkgs/R/") )
library(Rsamtools)
library(seqbias)
library(ggplot2)
library(rtracklayer)
library(seqLogo)
#
sb.atac_pos <- seqbias_atac(input_intervals = peak_intervals,
                            ref_seq = ref_seq,
                            sample_bam = sample_bam_pos,
                            indet_strands = "+")
sb.atac_neg <- seqbias_atac(input_intervals = peak_intervals,
                            ref_seq = ref_seq,
                            sample_bam = sample_bam_neg,
                            indet_strands = "-")
save(sb.atac_pos,file = "sb.atac_pos.Rlist")
save(sb.atac_neg,file = "sb.atac_neg.Rlist")
#
sb.atac_pos.ppm <- make_pwm(freqs = sb.atac_pos$freqs,
                            name = paste(name,"sb.atac_pos",sep = "."))
sb.atac_neg.ppm <- make_pwm(freqs = sb.atac_neg$freqs,
                            name = paste(name,"sb.atac_neg",sep = "."))
save(sb.atac_pos.ppm,file = "sb.atac_pos.ppm.Robject")
save(sb.atac_neg.ppm,file = "sb.atac_neg.ppm.Robject")
#