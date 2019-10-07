#
seqbias_atac <- function(input_intervals,
                         ref_seq,
                         sample_bam,
                         n,
                         m,
                         binarycount=TRUE,
                         k=1,
                         selectchrs,
                         L=100,
                         R=100,
                         indet_strands='+',
                         name,
                         ...){
  #How to make sure the FIMO BED is reformatted as [chr, start, end, name, score, strand]?
  
  #1. Require packages
  require(Rsamtools)
  require(seqbias)
  require(ggplot2)
  require(rtracklayer)
  require(GenomicRanges)
  
  #2. Set the reference sequence--a FASTA file--and the alignments--a BAM file
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  
  #3. Open the reference
  print(ref_fn)
  ref_f <- FaFile( ref_fn)
  # ref_f <- scanFa( ref_fn)
  open.FaFile( ref_f )
  print(ref_f)
  
  #4. Extract the sequences in the reference and take a random sample of n sequences that are m bases long
  # ref_seqs <- scanFaIndex( ref_f ,as=c("GRangesList", "GRanges")) # this is throwing an error...
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

  #7. Prepare and submit the output
  output <- list(ref_seqs,I,seqs,counts,freqs)
  names(output) <- c("ref_seqs","I","seqs","counts","freqs")
  return(output)
}
#
kmer_freq_plot <- function(freqs,
                           name) {
  print("Plotting k-mer frequencies")
  ylimit <- max(freqs$freq)
  P <- qplot( x = pos,
              y = freq,
              ylim = c(0.0,ylimit+0.01),
              color = seq,
              data = freqs,
              geom = "line" )
  P <- P + facet_grid( seq ~ . )
  return(P)
}
#
seqbias_compensation <- function(ref_seq,
                                 sample_bam,
                                 L,
                                 R,
                                 ...) {
  #fit the seqbias model. Specify an interval that is L bases to the left of the read start and R bases to the right
  ref_fn <- ref_seq
  reads_fn <- sample_bam
  set.seed(0)
  sb <- seqbias.fit( ref_fn, reads_fn, L = L, R = R )
  return(sb)
}
#
estimate_bias <- function(sb,
                          I,
                          k,
                          seqs,
                          counts,
                          L=100,
                          R=100,
                          name,
                          ...) {
  #1. Estimate bias
  bias <- seqbias.predict( sb, I )
  
  #2. Adjust the counts
  counts <- mapply( FUN = '/', counts, bias, SIMPLIFY = F )
  
  #3. Recalculate the k-mer frequencies. Use the same value of k as in seqbias_atac
  freqs <- kmer.freq( seqs = seqs, counts = counts,L=L,R=R )
  
  #4. Prepare output
  output <- list(bias,counts,freqs)
  names(output) <- c("bias","counts","freqs")
  return(output)
}
#
make_counts_matrix <- function(seqbias_counts,
                               ...) {
  #make the counts matrix
  sb.counts <- as.data.frame(t(as.data.frame(lapply(X = seqbias_counts,unlist))))
  sb.counts.sum <- apply(X = sb.counts,MARGIN = 2,sum)
  sb.counts.freq <- sb.counts.sum/sum(sb.counts.sum)
  
  #output
  output <- list(sb.counts,sb.counts.sum,sb.counts.freq)
  names(output) <- c("sb.counts","sb.counts.sum","sb.counts.freq")
  return(output)
}
#
adjust_phase <- function(posfreq,
                         negfreq,
                         motiflen,
                         clipbp=5,
                         ...) {
  #Set the position indices
  indices <- c(seq(-100,-1,1),rep(0,motiflen),seq(1,99,1)) #the 99 accounts for the fact that SeqBias neglects to count the last base of the motif from the zero-based coordinate system of the intervals
  
  #The negative needs to be moved up by one before being reversed
  neg.adj <- negfreq[2:length(negfreq)]
  rev.neg.adj <- rev(neg.adj)
  
  #Clip off the last point in positive to keep the lengths consistent. Do so for the indices as well
  pos.adj <- posfreq[1:(length(posfreq)-1)]
  indices <- indices[1:(length(posfreq)-1)]
  
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
  sum.adj.shift.clipped <- sum.adj.shift[(clipbp+1+2):(length(sum.adj.shift)-clipbp)]
  
  #Adjust the indices accordingly
  #indices.new <- indices[(clipbp+1):(length(indices)-5)] #If we are clipping the first 4 bases from the beginning, then the length should be cut by -4, not -5. This mistake made the right non-motif sequence 95 instead of 96.
  indices.new <- indices[(clipbp+1+2):(length(indices)-clipbp)]
  
  #return the matrix
  cutsite.freq <- as.data.frame(cbind(pos=indices.new,freq=sum.adj.shift.clipped))
  output <- list(cutsite.freq,clipbp)
  names(output) <- c("cutsite.freq","clipbp")
  return(output)
}
#
pos_neg_counts <- function(sb.poscounts,
                           sb.negcounts,
                           motiflen,
                           name,
                           ...) {
  #1. get the positive counts
  poscounts <- make_counts_matrix(sb.poscounts)
  negcounts <- make_counts_matrix(sb.negcounts)
  
  #2. Shift the positions of the cut-site frequencies such that the positive and negative strands are in phase
  shifted.freq <- adjust_phase(posfreq = poscounts$sb.counts.freq,
                               negfreq = negcounts$sb.counts.freq,
                               motiflen = motiflen,
                               ...)
  shifted.freq.list <- shifted.freq$cutsite.freq
  clipbp.list <- shifted.freq$clipbp
  
  # #3. Plot the result -- uncomment this if you wish, but this is not necessary as the footprint is plotted below
  # pdf(paste(name,"footprint","pdf",sep="."))
  # plot(shifted.freq.list$freq,xaxt="n",type="l",ylab="Cut site probability",xlab="Position from motif")
  # axis(1,at=1:length(shifted.freq.list$freq),labels = shifted.freq.list$pos)
  # print(100-clipbp.list+1)
  # print(100-clipbp.list+motiflen)
  # abline(v=100-clipbp.list+1-2,lty=2)
  # abline(v=100-clipbp.list+motiflen-2,lty=2)
  # dev.off()
  
  #4. Return the output
  output <- shifted.freq
  return(output)
}
#
make_pwm <- function(freqs,
                     name,
                     ...) {
  #This function takes the frequencies calculated above and transform them into a PWM
  
  #1. collapse the nucleotide frequency table that seqbias spits out into a proper PPM (position probability matrix)
  require(seqLogo)
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
  
  #3. convert the ppm into TRANSFAC format for WEBLOGO
  freqs.pwm.transfac <- cbind(PO=seq(1,ncol(freqs.pwm@pwm)),t(freqs.pwm@pwm))
  ppmfile <- paste(name,"ppm",sep=".")
  
  #4. save output
  write.table(freqs.pwm.transfac,file=ppmfile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  lines <- readLines(ppmfile)
  writeLines(c("ID","BF", lines), con = ppmfile)
  return(freqs.pwm.transfac)
}
#
kl_divergence <- function(freqs_bkgd,
                          freqs_fgd,
                          name,...) {
  #' At each position in the intervals studied, the KL divergence of the k-mers found is calculated.
  #' The frequencies of the k-mers at the position represent a probability distribution, and so the KL
  #' divergence is calculated to determine how different the k-mer representation is at each genomic 
  #' position between the background (i.e. before correction) and the foreground (i.e. after correction)
  
  #1. make sure that the foreground and background frequency lengths agree
  if (length(freqs_bkgd) != length(freqs_fgd)) {
    print("background and foreground frequencies must be calculated over the same length of bases")
    return(NULL)
  }
  
  #2. specify the Kullback-Leibler divergence function
  kl.fn <- function(f.bg,f.obs) { 
    kl_x <- f.obs*log2(f.obs/f.bg) + f.bg*log2(f.bg/f.obs) #calculate the divergence between the two distributions for a particular state in the distribution
    kl_divergence <- sum(kl_x) #sum over the instances
    return(kl_divergence)
  }
  
  #3. Apply the KL divergence to each genomic positional unit where a set of k-mers occur. 
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
    kl.fn(f.bg = bg,f.obs = fg)
  }, indices,length(kmers))
  pos <- freqs_bkgd$pos[seq(1,length(freqs_bkgd$pos),length(kmers))]
  kl.df <- cbind(pos,kl)
  
  #4. Plot the result
  pdf(paste(name,"kl_divergence","pdf",sep = "."))
  plot(kl.df,type="l")
  dev.off()
  return(kl.df)
}
#
seqbias_before <- function(ref_seq,
                           posbam,
                           negbam,
                           input_intervals,
                           bincount,
                           name,
                           ...) {
  #for this one, should I consolidate all of these into a single object?
  
  #1. positive strand
  seqbias.inputs.pos <- seqbias_atac(ref_seq=ref_seq,
                                     input_intervals=input_intervals,
                                     sample_bam=posbam,
                                     binarycount = bincount,
                                     indet_strands = "+",
                                     name = paste(name,"pos",sep="_"))
  sb.ppm.pos.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
  
  #2. negative strand
  seqbias.inputs.neg <- seqbias_atac(ref_seq=ref_seq,
                                     input_intervals=input_intervals,
                                     sample_bam=sample_bam_neg,
                                     binarycount = bincount,
                                     indet_strands = "-",
                                     name = paste(name,"neg",sep="_"))
  sb.ppm.neg.before <- make_pwm(freqs = seqbias.inputs.pos$freqs,name = name)
  
  #3. save the lists of objects
  saveRDS(seqbias.inputs.pos,file = paste(name,"seqbias.inputs.pos","rds",sep="."))
  saveRDS(seqbias.inputs.neg,file = paste(name,"seqbias.inputs.neg","rds",sep="."))

  #4. save the PPM's
  saveRDS(sb.ppm.pos.before,file = paste(name,"sb.ppm.pos.before","rds",sep="."))
  saveRDS(sb.ppm.neg.before,file = paste(name,"sb.ppm.neg.before","rds",sep="."))
  
  #5. prepare the output
  posout <- list(seqbias.inputs.pos,sb.ppm.pos.before)
  names(posout) <- c("seqbias","ppm")
  negout <- list(seqbias.inputs.neg,sb.ppm.neg.before)
  names(negout) <- c("seqbias","ppm")
  output <- list(posout,negout)
  names(output) <- c("pos","neg")
  #saveRDS(object = output,file = paste(name,"sb.before","rds",sep="."))
  return(output)
}
#
seqbias_after <- function(seqbias.inputs.pos,
                          seqbias.inputs.neg,
                          posmodel,
                          negmodel,
                          name,
                          ...) {
  #1. positive strand
  seqbias.corrected.pos <- estimate_bias(sb=posmodel,
                                         I=seqbias.inputs.pos$"I",
                                         seqs=seqbias.inputs.pos$"seqs",
                                         counts=seqbias.inputs.pos$"counts",
                                         name = paste(name,"pos",sep="_"))
  sb.ppm.pos.after <- make_pwm(freqs = seqbias.corrected.pos$freqs,name = name)
  
  #2. negative strand
  seqbias.corrected.neg <- estimate_bias(sb=negmodel,
                                         I=seqbias.inputs.neg$"I",
                                         seqs=seqbias.inputs.neg$"seqs",
                                         counts=seqbias.inputs.neg$"counts",
                                         name = paste(name,"neg",sep = "_"))
  sb.ppm.neg.after <- make_pwm(freqs = seqbias.corrected.neg$freqs,name = name)
  
  #3. save the lists of objects
  saveRDS(seqbias.corrected.pos,file = paste(name,"seqbias.corrected.pos","rds",sep="."))
  saveRDS(seqbias.corrected.neg,file = paste(name,"seqbias.corrected.neg","rds",sep="."))

  #4. save the PPM's
  saveRDS(sb.ppm.pos.after,file = paste(name,"sb.ppm.pos.after","rds",sep="."))
  saveRDS(sb.ppm.neg.after,file = paste(name,"sb.ppm.neg.after","rds",sep="."))
  
  #5. prepare the output
  posout <- list(seqbias.corrected.pos,sb.ppm.pos.after)
  names(posout) <- c("seqbias_atac","ppm")
  negout <- list(seqbias.corrected.neg,sb.ppm.neg.after)
  names(negout) <- c("seqbias_atac","ppm")
  output <- list(posout,negout)
  names(output) <- c("pos","neg")
  #saveRDS(object = output,file = paste(name,"sb.after","rds",sep="."))
  return(output)
}
#' Might it be a good idea to consolidate the seqbias_before and seqbias_after functions into one? 
#' if a seqbias model has been inputted, then the function would be able to run the correction with the model.
#
measure_deflection.spline <- function(footprint,
                                      name,
                                      motiflen,
                                      nonmotiflen=95,
                                      extralen=19,
                                      convolute=TRUE,
                                      confidence.level=0.99) {
  #' The strategy: define a natural spline on the predictors of the footprint (in this case, the serial position of bases). For knots, use five positions:
  #' 1. the start of the 19 bases to the left of the motif
  #' 2. the left of the motif
  #' 3. the low-end estimate of the motif midpoint
  #' 4. the right of the motif
  #' 5. the end of the 19 bases to the right of the motif
  #' Fit the relationship between the footprint transposition frequencies and the natural spline
  library(splines)
  
  #1. Add the footprint to the reverse of itself and divide the result by 2. This gives a symmetric footprint with which we can estimate
  if (convolute) {
    footprint.symm <- (footprint$freq + rev(footprint$freq))/2
  } else {
    footprint.symm <- footprint$freq
  }
  index <- seq_along(footprint.symm)
  print(index)
  print(footprint.symm)
  
  #2. Set the knots for the spline model
  set.knots <- c(nonmotiflen,
                 nonmotiflen+motiflen)
  if (extralen) {
    set.knots <- c(nonmotiflen - extralen,
                   set.knots,
                   set.knots[length(set.knots)] + extralen)
  }
  mdptpos <- floor(median(seq_along(set.knots))+1)
  mdptknot <- ceiling((set.knots[mdptpos-1] + set.knots[mdptpos])/2)
  set.knots <- c(set.knots[c(1:(mdptpos-1))],
                 mdptknot,
                 set.knots[c(mdptpos:length(set.knots))])
  #for the knots, we do not want to include the last index of the footprint. Otherwise, the model will attempt to estimate the curve for values that do not exist (as they would go off the end of the dataset)
  set.knots <- set.knots[!set.knots %in% c(max(index))]
  print(motiflen)
  print(set.knots)
  
  #3. Set the spline
  spline <- ns(index,knots = set.knots)
  #fit the model
  spline.model <- lm(footprint.symm ~ spline)
  #generate the prediction
  spline.model.predict <- predict(spline.model,
                                  data.frame(x=index),
                                  interval='confidence',
                                  level=confidence.level)
  #4. generate the plot
  pdf(paste(name,"pdf",sep = "."))
  #first, plot the footprint
  plot(index,footprint$freq,type="l")
  #now plot the different parts of the spline model
  lines(index,spline.model.predict[,1],col="red")
  lines(index,spline.model.predict[,2],col="red",lty=2)
  lines(index,spline.model.predict[,3],col="red",lty=2)
  abline(v = set.knots,col="blue",lty=2)
  
  #5. calculate deflection
  nonmotif.est <- min(spline.model.predict[,1][nonmotiflen-extralen],spline.model.predict[,1][nonmotiflen+motiflen+1+extralen])
  abline(h = nonmotif.est,lty=1,col="dark green")
  #motif.est <- spline.model.predict[,1][mdptknot]
  motif.est <- as.numeric(spline.model.predict[,1][mdptknot])
  abline(h = motif.est,lty=1,col="dark green")
  dev.off()
  deflection <- nonmotif.est - motif.est
  
  #6. Prepare output
  deflection_est <- list(nonmotif.est,motif.est,deflection)
  names(deflection_est) <- c("nonmotif.est","motif.est","deflection")
  output <- list(footprint.symm = footprint.symm,
                 footprint = footprint,
                 motiflen = motiflen,
                 set.knots = set.knots,
                 spline = spline,
                 spline.model = spline.model,
                 spline.model.predict = spline.model.predict,
                 deflection_est = deflection_est)
  saveRDS(output,file = paste(name,"deflection_model","rds",sep="."))
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
posmodel.fn <- args[7]
negmodel.fn <- args[8]
# cliplen <- args[7] # this was initially an argument, but 
cliplen <- 4
cliplen <- ceiling(as.numeric(cliplen)) #convert the cliplen to a numeric


# .libPaths( c( "/n/data2/mgh/ragon/pillai/nonLSFpkgs/R/") )
library(Rsamtools)
library(seqbias)
library(ggplot2)
library(rtracklayer)
library(seqLogo)

#2. Collect the models to be respectively applied to the positive-stranded reads and the negative-stranded reads
posmodel <- seqbias.load(ref_fn = ref_seq,
                         model_fn = posmodel.fn)
print(paste("Loaded positive-strand model at ",posmodel.fn))

#
negmodel <- seqbias.load(ref_fn = ref_seq,
                         model_fn = negmodel.fn)
print(paste("Loaded negative-strand model at ",negmodel.fn))

#
#3. Get the important seqbias objects prior to any correction
print("Obtain seqbias counts and frequencies...")
bincount=TRUE
before <- seqbias_before(ref_seq = ref_seq,
                         posbam = sample_bam_pos,
                         negbam = sample_bam_neg,
                         bincount = bincount,
                         input_intervals = input_intervals,
                         name = paste(name,"beforecorrection",sep="_"))

#4. Impose the corrections from the inputted models and the two hard-coded models: the Shendure chr22 naked DNA model, and the Buenrostro ATAC-seq model
print("Imposing trained Seqbias model...")
after.model <- seqbias_after(seqbias.inputs.pos = before$pos$seqbias,
                                          seqbias.inputs.neg = before$neg$seqbias,
                                          posmodel = posmodel,
                                          negmodel = negmodel,
                                          name = paste(name,"seqbiasmodel",sep = "_"))
print("...model imposed")

#5. Generate counts matrices for the positive and negative stranded objects of the pre- and post-corrected runs. Produce a footprint plot for each. Generate a spline model
#When adjusting phase, we will need to trim off 4 bases from the positive and negative strands (the latter will also require an additional base trim to account for 9-base overhang that Tn5 imposes)
clip <- cliplen
print(paste("input clip:",clip))
#Before correction
print("Generate footprints prior to correction")
#generate the counts track
before.footprint <- pos_neg_counts(sb.poscounts = before$pos$seqbias$counts,
                                   sb.negcounts = before$neg$seqbias$counts,
                                   motiflen = motiflen,
                                   name = paste(name,"before",sep="."),
                                   clipbp=clip)
#create a spline model
print(paste("clip from function",before.footprint$clipbp))
nonmotiflen <- 100-before.footprint$clipbp-2
before.footprint.spline <- measure_deflection.spline(footprint = before.footprint$cutsite.freq,
                                                     name = paste(name,"before",sep="."),
                                                     motiflen = motiflen,
                                                     nonmotiflen = nonmotiflen)
#write the counts track to an output
write.table(x = before.footprint$cutsite.freq,
            file = paste(name,"before_correction","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
before.sb <- c(before,before.footprint,before.footprint.spline)
saveRDS(before.sb,file = paste(name,"before_correction","rds",sep = "."))

#After correction with the inputed models
#generate the counts track
print("Generate footprints after correction with inputted models")
after.model.footprint <- pos_neg_counts(sb.poscounts = after.model$pos$seqbias$counts,
                                        sb.negcounts = after.model$neg$seqbias$counts,
                                        motiflen = motiflen,
                                        name = paste(name,"after_model",sep="."),
                                        clipbp=clip)
#create a spline model
after.model.spline <- measure_deflection.spline(footprint = after.model.footprint$cutsite.freq,
                                                name = paste(name,"after_model",sep="."),
                                                motiflen = motiflen,
                                                nonmotiflen = nonmotiflen)
#write the counts track to an output
write.table(x = after.model.footprint$cutsite.freq,
            file = paste(name,"after_model","shifted_freq","txt",sep = "."),
            quote = FALSE,sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
after.model.sb <- c(after.buenrostromodel,after.model.footprint,after.model.spline)
saveRDS(after.model.sb,file = paste(name,"after_model","rds",sep = "."))

#6. Plot
#indicate the motif boundaries

#
vline1 <- 100-before.footprint$clipbp-2+1
vline2 <- 100-before.footprint$clipbp-2+motiflen

#
print("Generate plots")
pdf(paste(name,"footprint","pdf",sep = "."),width = 28,height = 7)
par(mfrow=c(1,2))

#Before correction
plot(before.footprint$cutsite.freq$freq,
     xaxt="n",
     type="l",
     ylab="Cut site probability",
     xlab="Position from motif",
     ylim=c(0,0.03),
     col="red",
     main="Before Correction")
abline(v = c(vline1,vline2),lty=2)
legend("topright",legend = c("before","after"),col = c("red","blue"),lty = 2)

#After correction with inputted models
plot(after.model.footprint$cutsite.freq$freq,
     xaxt="n",
     type="l",
     ylim=c(0,0.03),
     col="blue",
     main=paste("Inputted seqbias model\nfor",posmodel.fn,"and\n",negmodel.fn,))
abline(v = c(vline1,vline2),lty=2)

#turn off the device
dev.off()
print("Done")

#Done