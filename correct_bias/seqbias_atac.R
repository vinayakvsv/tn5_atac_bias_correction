#the non-example seqbias_atac needs to use rtracklayer to input a BED file
seqbias_atac <- function(input_intervals,ref_seq,sample_bam,n,m,binarycount=TRUE,k=1,selectchrs,L=50,R=50,indet_strands='+',name){
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
  #print(P)
  output <- list(ref_seqs,I,seqs,counts,freqs,P)
  saveRDS(output,file = paste(name,"rds",sep="."))
  names(output) <- c("ref_seqs","I","seqs","counts","freqs","P")
  return(output)
}
#
seqbias_compensation <- function(ref_seq,sample_bam,L,R) {
  #fit the seqbias model. Specify an interval that is L bases to the left of the read start and R bases to the right
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  sb <- seqbias.fit( ref_fn, reads_fn, L = L, R = R )
  return(sb)
}
#
estimate_bias <- function(sb.in,I.in,seqs.in,counts.in,k.in=1,L=50,R=50,name) {
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
  #print(P)
  output <- list(bias,counts.adj,freqs.adj)
  names(output) <- c("bias","counts.adj","freqs.adj")
  saveRDS(output,paste(name,"rds",sep="."))
  return(output)
}
#
kl_divergence <- function(freqs_bkgd,freqs_fgd,name) {
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
  plot(kl.df,type="l",ylim=c(0,0.5))
  dev.off()
  return(kl.df)
}
#
make_counts_matrix_old <- function(seqbias_counts,motiflen) {
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
make_counts_matrix <- function(seqbias_counts) {
  #
  sb.counts <- as.data.frame(t(as.data.frame(lapply(X = seqbias_counts,unlist))))
  sb.counts.sum <- apply(X = sb.counts,MARGIN = 2,sum)
  sb.counts.freq <- sb.counts.sum/sum(sb.counts.sum)
  #output
  output <- list(sb.counts,sb.counts.sum,sb.counts.freq)
  names(output) <- c("sb.counts","sb.counts.sum","sb.counts.freq")
  return(output)
}
#
make_pwm <- function(freqs,name) {
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
#
adjust_phase <- function(posfreq,negfreq,motiflen) {
  #Set the position indices
  indices <- c(seq(-100,-1,1),rep(0,motiflen),seq(1,100,1))
  #The negative needs to be moved up by one before being reversed
  neg.adj <- negfreq[2:length(negfreq)]
  rev.neg.adj <- rev(neg.adj)
  #Clip off the last point in positive to keep the lengths consistent
  pos.adj <- posfreq[1:(length(posfreq)-1)]
  #The positive strand needs to be shifted downstream by 4 bases. This is achieved by clipping the last four bases and adding four zeroes to the beginning
  pos.adj.shift <- c(rep(0,4),pos.adj[1:(length(pos.adj)-4)])
  #The negative strand needs to be shifted upstream by 4 bases. This is achieved by clipping the first four bases and adding four zeroes to the end
  rev.neg.adj.shift <- c(rev.neg.adj[5:(length(rev.neg.adj))],rep(0,4))
  #Sum up the two strands
  sum.adj.shift <- pos.adj.shift + rev.neg.adj.shift
  #Clip of the first and last four numbers
  sum.adj.shift.clipped <- sum.adj.shift[5:(length(sum.adj.shift)-4)]
  #Adjust the indices accordingly
  indices.new <- indices[5:(length(indices)-5)]
  #return the matrix
  cutsite.freq <- as.data.frame(cbind(pos=indices.new,freq=sum.adj.shift.clipped))
  return(cutsite.freq)
}
#
plot_footprint <- function(sb.freq,sb.freq.positions,motiflen,name) {
  #plot the footprint
  pdf(paste(name,"footprint","pdf",sep="."))
  plot(sb.freq,xaxt="n",type="l",ylab="Cut site probability",xlab="Position from motif")
  #lim <- (length(sb.counts.freq) - motiflen + 5)/2
  #print(lim)
  #print(length(sb.counts.freq))
  axis(1,at=1:length(sb.freq),labels = sb.freq.positions)
  abline(v=97,lty=2)
  abline(v=96+motiflen,lty=2)
  #axis(1,at=1:length(sb.counts.freq),labels=c(seq(-lim,-1,1),rep(0,motiflen),seq(1,(lim-5),1)))
  #abline(v=lim+1,lty=2)
  #abline(v=lim+motiflen,lty=2)
  dev.off()
}
#
pos_neg_counts <- function(sb.poscounts,sb.negcounts,motiflen,name) {
  #1. get the positive counts
  poscounts <- make_counts_matrix(sb.poscounts)
  negcounts <- make_counts_matrix(sb.negcounts)
  #####
  #   #2. reverse the read sums per position of negative and apply a five-base correction to both positive and negative. 
  #   #This is achieved by taking away the last 5 positions of positive, the first 5 bases of reversed negative, and adding the two vectors together
  #   posrange <- c(1:(length(poscounts$sb.counts.sum)-5))
  #   negrange <- c(6:length(negcounts$sb.counts.sum))
  #   print("lens")
  #   print(length(poscounts$sb.counts.sum))
  #   print(length(negcounts$sb.counts.sum))
  #   print("motiflen")
  #   print(motiflen)
  #   print("ranges")
  #   print(posrange)
  #   print(length(posrange))
  #   print(negrange)
  #   print(length(negrange))
  #   poscounts.shifted <- poscounts$sb.counts.sum[posrange]
  #   negcounts.shifted <- rev(negcounts$sb.counts.sum[negrange])
  #   #or negcounts.shifted <- rev(negcounts$sb.counts.sum)[negrange]?
  #   sum.shifted <- poscounts.shifted + negcounts.shifted
  #   print("shifted")
  #   print(length(poscounts.shifted))
  #   print(length(negcounts.shifted))
  #   print(length(sum.shifted))
  #   #3. Get the cut site frequencies at each position
  #   freq.shifted <- sum.shifted/sum(sum.shifted)
  #####
  #2. Shift the positions of the cut-site frequencies such that the positive and negative strands are in phase
  shifted.freq <- adjust_phase(posfreq = poscounts$sb.counts.freq,
                               negfreq = negcounts$sb.counts.freq,
                               motiflen = motiflen)
  #3. Plot the result
  plot_footprint(sb.freq = shifted.freq$freq,
                 sb.freq.positions = shifted.freq$pos,
                 motiflen = motiflen,
                 name = name)
  #4. Return the output
  output <- shifted.freq
  return(output)
}
#
#' for preparing a model
#' 1. 
#
#####
count.reads.edit <- function (reads_fn, I, sb = NULL, binary = FALSE, sum.counts = FALSE) {
  require(GenomicRanges)
  stopifnot(is(I, "GRanges"))
  if (!is.null(sb) & class(sb) != "seqbias") {
    stop("The 'sb' parameter of count.reads must be a seqbias class.")
  }
  bam_ptr <- .Call("seqbias_open_bam", path.expand(reads_fn), 
                   PACKAGE = "seqbias")
  counts <- tapply(I, INDEX = 1:length(I), FUN = function(x) .Call("seqbias_count_reads", 
                                                                   sb, bam_ptr, as.character(seqnames(x)), start(x), end(x), 
                                                                   as.character(strand(x)), sum.counts, PACKAGE = "seqbias"))
  #print(counts)
  print(binary)
  if (binary) 
    lapply(counts, FUN = function(c) as.integer(c > 0))
  else counts
}
#####
#OLD FUNCTIONS
seqbias_atac_example <- function(ref_seq,sample_bam,n,m,binarycount=T,k=1,selectchrs,L=50,R=50) {
  #1. Require packages
  require(Rsamtools)
  require(seqbias)
  require(ggplot2)
  #2. Set the reference sequence--a FASTA file--and the alignments--a BAM file
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  #3. Open the reference
  ref_f <- FaFile( ref_fn)
  open.FaFile( ref_f )
  #4. Extract the sequences in the reference and take a random sample of n sequences that are m bases long
  #which <- RangesList(chr22=IRanges(start=1, end=50818468))
  ref_seqs <- scanFaIndex( ref_f )
  print(ref_seqs)
  #ref_seqs[38] contains chr22
  #ref_seqs <- seqbias_beforemodel$ref_seqs[seqnames(seqbias_beforemodel$ref_seqs) %in% selectchrs]
  ref_seqs <- ref_seqs[seqnames(ref_seqs) %in% selectchrs]
  I <- random.intervals( ref_seqs, n = n, m = m ) #This can be replaced with an import of the sequences in an ATAC-seq interval
  seqs <- scanFa( ref_f, I )
  
  
  #Make sure that you get the reverse complementary sequence of the negative-stranded sample sequences
  neg_idx <- as.logical( I@strand == '-')
  seqs[neg_idx] <- reverseComplement( seqs[neg_idx] )
  print(seqs)
  
  
  #5. Count the number of reads within the sample sequences. The counts are binary--that is, as long as one match is found, the count is 1
  counts <- count.reads( reads_fn, I, binary = binarycount )
  #6. Get the frequencies of particular k-mers in the sample sequences. By default, k=1
  freqs <- kmer.freq( seqs, counts ,k=k,L = L,R = R)
  #7. Plot the frequencies of the k-mers
  P <- qplot( x = pos,
              y = freq,
              ylim = c(0.15,0.4),
              color = seq,
              data = freqs,
              geom = "line" )
  P <- P + facet_grid( seq ~ . )
  print(P)
  output <- list(ref_seqs,I,seqs,counts,freqs)
  names(output) <- c("ref_seqs","I","seqs","counts","freqs")
  return(output)
}
make_pwm_old <- function(freqs) {
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
  return(freqs.pwm)
}
#Suggested inputs
#input <- "~/Desktop/mnt.transfer/analysis/bcell_subsets/gc/gc.atac-seq/Run_195/human.samples/enhancer.clusters/en-tss/window-500kb-wtss/up-gc-vs-nongc/up-gc-vs-nongc.non-tss_peaks.bed"
#hg38_ref <- "~/Desktop/mnt.transfer/referenceGenomes/hg38_ucsc/hg38.fa"
#test_bam <- "~/Desktop/mnt.transfer/analysis/bcell_subsets/gc/gc.atac-seq/Run_195/human.samples/alignments.q10filtered.sorted.hg38-autosom-sex-chr/S0004.q10filtered.sorted.hg38-autosom-sex-chr.bam"
#
#sb.atac <- seqbias_atac(input_intervals = input,ref_seq = hg38_ref,sample_bam = test_bam)
#sb.gc_004.model <- seqbias_compensation(ref_seq = hg38_ref,sample_bam = test_bam,L = 15,R = 15)
#bias.est.gc_s004 <- estimate_bias(sb = sb.gc_004.model,I = sb.atac$I,k = 1,seqs = sb.atac$seqs,counts = sb.atac$counts)
#