#the non-example seqbias_atac needs to use rtracklayer to input a BED file
seqbias_atac_old <- function(input_intervals,ref_seq,sample_bam,n,m,binarycount=T,k=1,selectchrs,L=100,R=100) {
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
  #4. Extract the sequences in the reference and take a random sample of n sequences that are m bases long
  #which <- RangesList(chr22=IRanges(start=1, end=50818468))
  ref_seqs <- scanFaIndex( ref_f )
  print(ref_seqs)
  #ref_seqs[38] contains chr22
  ref_seqs <- seqbias_beforemodel$ref_seqs[seqnames(seqbias_beforemodel$ref_seqs) %in% selectchrs]
  #I <- random.intervals( ref_seqs, n = n, m = m ) #This can be replaced with an import of the sequences in an ATAC-seq interval
  #input_intervals is a BED file with the genomic intervals in which you would like to correct for count biases
  I <- import(con = input_intervals)
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
              ylim = c(0.0,0.5),
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
seqbias_atac <- function(input_intervals,ref_seq,sample_bam,n,m,binarycount=TRUE,k=1,selectchrs,L=100,R=100,indet_strands='+',name){
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
  pdf(paste(name,"kmer_freq_before_correction","pdf",sep="."))
  print(P)
  dev.off()
  output <- list(ref_seqs,I,seqs,counts,freqs,P)
  names(output) <- c("ref_seqs","I","seqs","counts","freqs",P)
  return(output)
}
#
seqbias_compensation <- function(ref_seq,sample_bam,L,R) {
  #fit the seqbias model. Specify an interval that is L bases to the left of the read start and R bases to the right
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  set.seed(0)
  sb <- seqbias.fit( ref_fn, reads_fn, L = L, R = R )
  return(sb)
}
#
estimate_bias <- function(sb,I,k,seqs,counts,L=100,R=100,name) {
  #1. Estimate bias
  bias <- seqbias.predict( sb, I )
  #2. Adjust the counts
  counts.adj <- mapply( FUN = '/', counts, bias, SIMPLIFY = F )
  #3. Recalculate the k-mer frequencies. Use the same value of k as in seqbias_atac
  freqs.adj <- kmer.freq( seqs = seqs, counts = counts.adj,L=L,R=R )
  P <- qplot( x = pos,
              y = freq,
              ylim = c(0.0,0.5),
              color = seq,
              data = freqs.adj,
              geom = "line" )
  P <- P + facet_grid( seq ~ . )
  pdf(paste(name,"kmer_freq_after_correction","pdf",sep="."))
  print(P)
  dev.off()
  output <- list(bias,counts.adj,freqs.adj)
  names(output) <- c("bias","counts.adj","freqs.adj")
  return(output)
}
#
make_counts_matrix_old <- function(seqbias_counts) {
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
adjust_phase_old <- function(posfreq,negfreq,motiflen) {
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
adjust_phase <- function(posfreq,
                         negfreq,
                         motiflen,
                         clipbp=4,
                         ...) {
  #Set the position indices
  indices <- c(seq(-100,-1,1),rep(0,motiflen),seq(1,99,1)) #the 99 accounts for the fact that SeqBias neglects to count the last base of the motif from the zero-based coordinate system of the intervals
  #The negative needs to be moved up by one before being reversed
  print(length(posfreq))
  print(length(negfreq))
  neg.adj <- negfreq[2:length(negfreq)]
  rev.neg.adj <- rev(neg.adj)
  #Clip off the last point in positive to keep the lengths consistent. Do so for the indices as well
  pos.adj <- posfreq[1:(length(posfreq)-1)]
  indices <- indices[1:(length(posfreq)-1)]
  #Are the above necessary after the new adjust_phase in seqbias_atac_orchestra_setmodels.R?
  
  #The positive strand needs to be shifted downstream by 4 bases. This is achieved by clipping the last four bases and adding four zeroes to the beginning
  pos.adj.shift <- c(rep(0,4),pos.adj[1:(length(pos.adj)-4)])
  #The negative strand needs to be shifted upstream by 4 bases. This is achieved by clipping the first four bases and adding four zeroes to the end
  rev.neg.adj.shift <- c(rev.neg.adj[5:(length(rev.neg.adj))],rep(0,4))
  #Sum up the two strands
  sum.adj.shift <- pos.adj.shift + rev.neg.adj.shift
  #Clip of the first and last 'clipbp' bases. 
  #The mininum possible value for clipbp is 4 bases (since we did enforce a 4-base correction to account for Tn5's staggered cut)
  clipbp[clipbp < 4] <- 4 
  print(paste("clip by",clipbp,"bases"))
  sum.adj.shift.clipped <- sum.adj.shift[(clipbp+1):(length(sum.adj.shift)-clipbp)]
  #Adjust the indices accordingly
  #indices.new <- indices[(clipbp+1):(length(indices)-5)] #If we are clipping the first 4 bases from the beginning, then the length should be cut by -4, not -5. This mistake made the right non-motif sequence 95 instead of 96.
  indices.new <- indices[(clipbp+1+2):(length(indices)-clipbp)]
  #indices.new <- indices[(clipbp+1):(length(indices)-clipbp)]
  #return the matrix
  print(paste("length of indices",length(indices)))
  print(paste("length of indices.new",length(indices.new)))
  print(indices.new)
  print(paste("length of sum.adj.shift",length(sum.adj.shift)))
  print(paste("length of sum.adj.shift.clipped",length(sum.adj.shift.clipped)))
  #print()
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
pos_neg_counts <- function(sb.poscounts,
                           sb.negcounts,
                           motiflen,
                           name) {
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
  print(shifted.freq)
  #3. Plot the result
  plot_footprint(sb.freq = shifted.freq$freq,
                 sb.freq.positions = shifted.freq$pos,
                 motiflen = motiflen,
                 name = name)
  #4. Return the output
  output <- shifted.freq
  return(output)
}
make_pwm <- function(freqs,name) {
  #' Take the frequencies calculated above and transform them into a PWM
  require(seqLogo)
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
  ppmfile <- paste(name,"ppm",sep=".")
  write.table(freqs.pwm.transfac,file=ppmfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  lines <- readLines(ppmfile)
  writeLines(c("ID","BF", lines), con = ppmfile)
  return(freqs.pwm.transfac)
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
  plot(kl.df,type="l")
  dev.off()
  return(kl.df)
}
#
seqbias_before <- function(ref_seq,posbam,negbam,input_intervals,bincount,name) {
  #positive strand
  seqbias.inputs.pos <- seqbias_atac(ref_seq=ref_seq,
                                     input_intervals=input_intervals,
                                     sample_bam=posbam,
                                     binarycount = bincount,
                                     indet_strands = "+",
                                     name = paste(name,"pos",sep="_"))
  sb.ppm.pos.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
  #negative strand
  seqbias.inputs.neg <- seqbias_atac(ref_seq=ref_seq,
                                     input_intervals=input_intervals,
                                     sample_bam=sample_bam_neg,
                                     binarycount = bincount,
                                     indet_strands = "-",
                                     name = paste(name,"neg",sep="_"))
  sb.ppm.neg.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
  #save the lists of objects
  save(seqbias.inputs.pos,file = paste(name,"seqbias.inputs.pos","Rlist",sep="."))
  save(seqbias.inputs.neg,file = paste(name,"seqbias.inputs.neg","Rlist",sep="."))
  #save the PPM's
  save(sb.ppm.pos.before,file = paste(name,"sb.ppm.pos.before","Robject",sep="."))
  save(sb.ppm.neg.before,file = paste(name,"sb.ppm.neg.before","Rlist",sep="."))
  #prepare the output
  posout <- list(seqbias.inputs.pos,sb.ppm.pos.before)
  names(posout) <- c("seqbias","ppm")
  negout <- list(seqbias.inputs.neg,sb.ppm.neg.before)
  names(negout) <- c("seqbias","ppm")
  output <- list(posout,negout)
  names(output) <- c("pos","neg")
  #return
  return(output)
}
#
seqbias_after <- function(seqbias.inputs.pos,seqbias.inputs.neg,posmodel,negmodel,name) {
  #pos
  seqbias.corrected.pos <- estimate_bias(sb=posmodel,
                                         I=seqbias.inputs.pos$"I",
                                         seqs=seqbias.inputs.pos$"seqs",
                                         counts=seqbias.inputs.pos$"counts",
                                         name = paste(name,"pos",sep="_"))
  sb.ppm.pos.after <- make_pwm(freqs = seqbias.corrected.pos$freqs.adj,name = name)
  #neg
  seqbias.corrected.neg <- estimate_bias(sb=negmodel,
                                         I=seqbias.inputs.neg$"I",
                                         seqs=seqbias.inputs.neg$"seqs",
                                         counts=seqbias.inputs.neg$"counts",
                                         name = paste(name,"neg",sep = "_"))
  sb.ppm.neg.after <- make_pwm(freqs = seqbias.corrected.neg$freqs.adj,name = name)
  #save outputs
  save(seqbias.corrected.pos,file = paste(name,"seqbias.corrected.pos","Rlist",sep="."))
  save(seqbias.corrected.neg,file = paste(name,"seqbias.corrected.neg","Rlist",sep="."))
  #save ppms
  save(sb.ppm.pos.after,file = paste(name,"sb.ppm.pos.after","Robject",sep="."))
  save(sb.ppm.neg.after,file = paste(name,"sb.ppm.neg.after","Rlist",sep="."))
  #prepare the output
  posout <- list(seqbias.corrected.pos,sb.ppm.pos.after)
  names(posout) <- c("seqbias_atac","ppm")
  negout <- list(seqbias.corrected.neg,sb.ppm.neg.after)
  names(negout) <- c("seqbias_atac","ppm")
  output <- list(posout,negout)
  names(output) <- c("pos","neg")
  #return
  return(output)
}
#
#####
#MAIN
#1. Take in command-line input and set the source of the libraries
args <- commandArgs(trailingOnly=TRUE)
ref_seq <- args[1] #the reference genome
sample_bam_pos <- args[2] #the positive BAM
sample_bam_neg <- args[3] #the negative BAM
input_intervals <- args[4] #the motif intervals
motiflen <- args[5] #the motif length
motiflen <- ceiling(as.numeric(motiflen)) #round the motif length up
print(motiflen)
name <- args[6] #the name associated with the job
posmodel <- args[7]
negmodel <- args[8]
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
#2. Collect the models to be respectively applied to the positive-stranded reads and the negative-stranded reads
print(ref_seq)
#posmodel.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.pos.dedupbam.L20.R20"
#negmodel.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.neg.dedupbam.L20.R20"
#
posmodel.fn <- posmodel
negmodel.fn <- negmodel
print(posmodel.fn)
print(negmodel.fn)
posmodel <- seqbias.load(ref_fn = ref_seq, model_fn = posmodel.fn)
negmodel <- seqbias.load(ref_fn = ref_seq, model_fn = negmodel.fn)
#
shendure_chr22_model.pos.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.pos.dedupbam.L20.R20"
#shendure_chr22_model.pos.fn <- buenrostropos
shendure_chr22_model.pos <- seqbias.load(ref_fn = ref_seq,
                                         model_fn = shendure_chr22_model.pos.fn)
print("Loaded /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.pos.dedupbam.L20.R2")
#
shendure_chr22_model.neg.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.neg.dedupbam.L20.R20"
#shendure_chr22_model.neg.fn <- buenrostroneg
shendure_chr22_model.neg <- seqbias.load(ref_fn = ref_seq,
                                         model_fn = shendure_chr22_model.neg.fn)
print("Loaded /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/shendure.chr22.neg.dedupbam.L20.R20")
#
buenrostro_atac_model.pos.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/sb.buenrostro_SRR891268_atac_dedup_posbam.posmodel.L20.R20"
#buenrostro_atac_model.pos.fn <- buenrostropos
buenrostro_atac_model.pos <- seqbias.load(ref_fn = ref_seq,
                                         model_fn = buenrostro_atac_model.pos.fn)
print("Loaded /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/sb.buenrostro_SRR891268_atac_dedup_posbam.posmodel.L20.R20")
#
buenrostro_atac_model.neg.fn <- "/n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/sb.buenrostro_SRR891268_atac_dedup_negbam.negmodel.L20.R20"
#buenrostro_atac_model.neg.fn <- buenrostroneg
buenrostro_atac_model.neg <- seqbias.load(ref_fn = ref_seq,
                                          model_fn = buenrostro_atac_model.neg.fn)
print("Loaded /n/data2/mgh/ragon/pillai/pipelines/quality-control_pipeline/seqbias/sb.buenrostro_SRR891268_atac_dedup_negbam.negmodel.L20.R20")
#
bincount=TRUE
#3. Get the important seqbias objects prior to any correction
print("Obtain seqbias counts and frequencies...")
before <- seqbias_before(ref_seq = ref_seq,
                         posbam = sample_bam_pos,
                         negbam = sample_bam_neg,
                         bincount = bincount,
                         input_intervals = input_intervals,
                         name = paste(name,"beforecorrection",sep="_"))

#4. Impose the corrections from the inputted models and the two hard-coded models: the Shendure chr22 naked DNA model, and the Buenrostro ATAC-seq model
#The inputted models
print("Impose inputted models")
after.inputmodel <- seqbias_after(seqbias.inputs.pos = before$pos$seqbias,
                       seqbias.inputs.neg = before$neg$seqbias,
                       posmodel = posmodel,
                       negmodel = negmodel,
                       name = paste(name,"postinputted",sep = "_"))
#The Shendure model
print("Impose Shendure chr22 model")
after.shendurechr22model <- seqbias_after(seqbias.inputs.pos = before$pos$seqbias,
                                  seqbias.inputs.neg = before$neg$seqbias,
                                  posmodel = shendure_chr22_model.pos,
                                  negmodel = shendure_chr22_model.neg,
                                  name = paste(name,"postshendurechr22",sep = "_"))
#The Buenrostro ATAC-seq model
print("Impose Buenrostro ATAC-seq model")
after.buenrostromodel <- seqbias_after(seqbias.inputs.pos = before$pos$seqbias,
                                       seqbias.inputs.neg = before$neg$seqbias,
                                       posmodel = buenrostro_atac_model.pos,
                                       negmodel = buenrostro_atac_model.neg,
                                       name = paste(name,"postbuenrostro",sep = "_"))
#####
# #3. Positive strand
# # Starting with the positive strand reads, collect the input interval sequences, the nucleotide frequencies, and the read start counts in each interval. Generate a PPM.
# seqbias.inputs.pos <- seqbias_atac(ref_seq=ref_seq,
#                                    input_intervals=input_intervals,
#                                    sample_bam=sample_bam_pos,
#                                    binarycount = bincount,
#                                    indet_strands = "+",
#                                    name = paste(name,"pos",sep="_"))
# sb.ppm.pos.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
# save(seqbias.inputs.pos,file = paste(name,"seqbias.inputs.pos","Rlist",sep="."))
# save(sb.ppm.pos.before,file = paste(name,"sb.ppm.pos.before","Robject",sep="."))
# # Apply the positive-strand correction and generate a PPM of the result
# seqbias.corrected.pos <- estimate_bias(sb=posmodel,
#                                        I=seqbias.inputs.pos$"I",
#                                        seqs=seqbias.inputs.pos$"seqs",
#                                        counts=seqbias.inputs.pos$"counts",
#                                        name = paste(name,"pos",sep="_"))
# sb.ppm.pos.after <- make_pwm(freqs = seqbias.corrected.pos$freqs.adj,name = name)
# save(seqbias.corrected.pos,file = paste(name,"seqbias.corrected.pos","Rlist",sep="."))
# save(sb.ppm.pos.after,file = paste(name,"sb.ppm.pos.after","Robject",sep="."))
# #
# #4. Negative strand
# seqbias.inputs.neg <- seqbias_atac(ref_seq=ref_seq,
#                                    input_intervals=input_intervals,
#                                    sample_bam=sample_bam_neg,
#                                    binarycount = bincount,
#                                    indet_strands = "-",
#                                    name = paste(name,"neg",sep="_"))
# sb.ppm.neg.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
# save(seqbias.inputs.neg,file = paste(name,"seqbias.inputs.neg","Rlist",sep="."))
# save(sb.ppm.neg.before,file = paste(name,"sb.ppm.neg.before","Rlist",sep="."))
# # Apply the negative-strand correction and generate a PPM of the result
# seqbias.corrected.neg <- estimate_bias(sb=negmodel,
#                                        I=seqbias.inputs.neg$"I",
#                                        seqs=seqbias.inputs.neg$"seqs",
#                                        counts=seqbias.inputs.neg$"counts",
#                                        name = paste(name,"neg",sep = "_"))
# sb.ppm.neg.after <- make_pwm(freqs = seqbias.corrected.neg$freqs.adj,name = name)
# save(seqbias.corrected.neg,file = paste(name,"seqbias.corrected.neg","Rlist",sep="."))
# save(sb.ppm.neg.after,file = paste(name,"sb.ppm.neg.after","Rlist",sep="."))
# #
#####
#5. Generate counts matrices for the positive and negative stranded objects of the pre- and post-corrected runs. Produce a footprint plot for each
#Before correction
print("Generate footprints prior to correction")
before.footprint <- pos_neg_counts(sb.poscounts = before$pos$seqbias$counts,
                                  sb.negcounts = before$neg$seqbias$counts,
                                  motiflen = motiflen,
                                  name = paste(name,"before",sep=""))
write.table(x = before.footprint,
            file = paste(name,"before_correction","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
#After correction with the inputted models
print("Generate footprints after correction with inputted models")
after.inputmodel.footprint <- pos_neg_counts(sb.poscounts = after.inputmodel$pos$seqbias$counts,
                                             sb.negcounts = after.inputmodel$neg$seqbias$counts,
                                             motiflen = motiflen,
                                             name = paste(name,"after",sep=""))
write.table(x = after.inputmodel.footprint,
            file = paste(name,"after_inputmodel","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
#After correction with the Buenrostro models
print("Generate footprints after correction with Buenrostro ATAC-seq models")
after.buenrostromodel.footprint <- pos_neg_counts(sb.poscounts = after.buenrostromodel$pos$seqbias$counts,
                                             sb.negcounts = after.buenrostromodel$neg$seqbias$counts,
                                             motiflen = motiflen,
                                             name = paste(name,"after",sep=""))
write.table(x = after.buenrostromodel.footprint,
            file = paste(name,"after_buenrostromodel","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
#After correction with the Shendure chr22 model
print("Generate footprints after correction with Shendure chr22 models")
after.shendurechr22model.footprint <- pos_neg_counts(sb.poscounts = after.shendurechr22model$pos$seqbias$counts,
                                             sb.negcounts = after.shendurechr22model$neg$seqbias$counts,
                                             motiflen = motiflen,
                                             name = paste(name,"after",sep=""))
write.table(x = after.shendurechr22model.footprint,
            file = paste(name,"after_shendurechr22model","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
# after.footprint <- pos_neg_counts(sb.poscounts = after.shendurechr22model$seqbias$counts,
#                                   sb.negcounts = seqbias.corrected.neg$counts,
#                                   motiflen = motiflen,
#                                   name = paste(name,"after",sep=""))
# write.table(x = after.footprint,
#             file = paste(name,"after_correction","shifted_freq","txt",sep = "."),
#             quote = FALSE,sep = "\t",
#             row.names = TRUE,
#             col.names = TRUE)
#Plot
print("Generate plots")
pdf(paste(name,"footprint","pdf",sep = "."),width = 28,height = 7)
par(mfrow=c(1,4))
#Before correction
plot(before.footprint$freq,
     xaxt="n",
     type="l",
     ylab="Cut site probability",
     xlab="Position from motif",
     ylim=c(0,0.03),
     col="red",
     main="Before Correction")
##These lines are off...
#abline(v = 97,lty=2)
#abline(v = 96+motiflen,lty=2)
#These are correct...NOTE: 
abline(v = 96,lty=2)
abline(v = 95+motiflen,lty=2)
legend("topright",legend = c("before","after"),col = c("red","blue"),lty = 2)
#After correction with inputted models
plot(after.inputmodel.footprint$freq,
     xaxt="n",
     type="l",
     ylim=c(0,0.03),
     col="blue",
     main="Shendure SRR1554094\nmotif interval coverage model")
# abline(v = 97,lty=2)
# abline(v = 96+motiflen,lty=2)
#These are correct...
abline(v = 96,lty=2)
abline(v = 95+motiflen,lty=2)
#After correction with Shendure models
plot(after.shendurechr22model.footprint$freq,
     xaxt="n",
     type="l",
     ylim=c(0,0.03),
     col="blue",
     main="Shendure SRR1554094\nchr22 WGS model")
# abline(v = 97,lty=2)
# abline(v = 96+motiflen,lty=2)
#These are correct...
abline(v = 96,lty=2)
abline(v = 95+motiflen,lty=2)
#After correction with Buenrostro models
plot(after.buenrostromodel.footprint$freq,
     xaxt="n",
     type="l",
     ylim=c(0,0.03),
     col="blue",
     main="Buenrostro SRR891268\nATAC-seq model ")
# abline(v = 97,lty=2)
# abline(v = 96+motiflen,lty=2)
#These are correct...
abline(v = 96,lty=2)
abline(v = 95+motiflen,lty=2)
#turn off the device
dev.off()
print("Done")
#
# save(seqbias.inputs,file = "seqbias.inputs.Robject")
# seqbias.save(sb = sb.model,fn = "sb.model.sbobject")
# save(seqbias.corrected,file = "seqbias.corrected.Robject")
#plot(sb.buenrostro_SRR891268_bam_pos.mef2c_genomewide.binT.counts$sb.counts.sum[c(1:106)]+rev(sb.buenrostro_SRR891268_bam_neg.mef2c_genomewide.binT.counts$sb.counts.sum)[c(6:111)],type="l")