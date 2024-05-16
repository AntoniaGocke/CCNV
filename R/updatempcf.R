#' Segment multiple data together using the partial least squares regression
#' as presented in XXX
#'
#' @param data dataframe of binned copy number data in the format (samples,
#'  values) with two additional columns. The first representing the chromosome,
#'  the second the position on the chromosome
#' @param gamma value representing a threshhold for the segmentation
#' @return A dataframe representing segmentation data, containing columns for 
#' the chromosome ,segmentation start ,segmentation end ,number of bins in that 
#' segment, mean value for that segment, each segmentation value for each sample
updatempcf <- function(data, gamma=5, frac1=0.15, frac2=0.15){
  
  chrom <- data[,1]
  position <- data[,2]
  sampleid <- colnames(data)[-c(1:2)]
  nSample <- length(sampleid)
  
  mpcf.names <- c("chrom","pos",sampleid)
  segments <- data.frame(matrix(nrow=0,ncol=nSample+5))
  seg.names <- c("chrom","start.pos","end.pos","n.probes",sampleid)
  colnames(segments) <- seg.names
  
  #Scale gamma according to the number of samples:
  gamma <- gamma*nSample
  chrom.list <- unique(chrom)
  
  #run multiPCF separately on each chromosome:
  for(c in 1:length(chrom.list)){
    
    probe.c <- which(chrom==chrom.list[c])
    pos.c <- position[probe.c]
    
    #get data for this chrom
    chrom.data <- data[which(data$Chrom == c),-c(1:2)]
    
    #Run multipcf:
    mpcf <- runFastMultiPCF(as.matrix(chrom.data),gamma=gamma, frac1, frac2)    #requires samples in columns, probes in rows
    
    #Information about segments:
    nSeg <- mpcf$nIntervals
    start0 <- mpcf$start0
    n.pos <- mpcf$length
    seg.mean <- t(mpcf$mean)  #get samples in columns
    posStart <- pos.c[start0]
    posEnd <- c(pos.c[start0-1],pos.c[length(probe.c)])
    
    #Chromosome number :
    chrid <- rep(unique(chrom[probe.c]),times=nSeg)
    
    #Round
    seg.mean <- round(seg.mean,digits=3)
    
    #Data frame:
    segments.c <- data.frame(chrid,posStart,posEnd,n.pos,seg.mean,stringsAsFactors=FALSE)
    colnames(segments.c) <- seg.names
    
    #Append results for this arm:
    segments <- rbind(segments,segments.c)
    
  }
  return(segments)
}


# Fast version 
runFastMultiPCF <- function (x, gamma, frac1, frac2) {
  mark <- rep(0, nrow(x))
  mark<-sawMarkM(x,frac1,frac2)
  dense <- compactMulti(t(x), mark)
  compPotts <- multiPCFcompact(dense$Nr, dense$Sum, gamma)
  return(list(length = compPotts$Lengde, start0 = compPotts$sta,
              mean = compPotts$mean, nIntervals = compPotts$nIntervals))
}

define_sawtooth <- function(half_length){
    filter <- rep(0, 2*half_length)
    for (k in 1:half_length) {
        filter[k] <- k / half_length
        filter[2 * half_length + 1 - k] <- -k / half_length
    }
    return(filter)
}

apply_filter <- function(x, filter){
    sawValue <- rep(0, nrow(x))
    L <- length(filter)/2 # filter size assumed to be even
    for (l in 1:(nrow(x) - 2 * L + 1)) {
        for (m in 1:ncol(x)) {
            diff <- crossprod(filter, x[l:(l + 2 * L - 1), m])
            sawValue[l + L - 1] <- sawValue[l + L - 1] + abs(diff)
        }
    }
    return(sawValue)
}

mark_points <- function(mark, sawValue, halflength, frac){
    limit <- quantile(sawValue, (1 - frac))
    for (l in 1:(length(sawValue) - 2 * halflength)) {
        if (sawValue[l + halflength - 1] > limit) {
            mark[l + halflength - 1] <- 1
        }
    }
    return(mark)
}

## sawtooth-filter for multiPCF - marks potential breakpoints. Uses two
## sawtoothfilters, one lang (length L) and one short (fixed length 6)
sawMarkM <- function(x, frac1, frac2){
    L <- 15
    filter <- define_sawtooth(L)
    filter2 <- define_sawtooth(3)
    sawValue <- apply_filter(x, filter)
    sawValue2 <- apply_filter(x, filter2)
    
    
    nrProbes <- nrow(x)
    nrSample <- ncol(x)
    mark <- rep(0, nrow(x))
    
    mark <- mark_points(mark, sawValue, L, frac1)
    mark <- mark_points(mark, sawValue2, L, frac2)
    
    mark[1:L] <- 1
    mark[nrow(x):nrow(x)+L-1] <- 1
    
    return(mark)
}

# function that accumulates numbers of observations and sums between potential breakpoints
compactMulti <- function(y, mark) {
  antGen <- ncol(y)
  antSample <- nrow(y)
  antMark <- sum(mark)
  ant <- rep(0, antMark)
  sum <- rep(0, antMark * antSample)
  dim(sum) <- c(antSample, antMark)
  pos <- 1
  oldPos <- 0
  count <- 1
  delSum <- rep(0, antSample)
  while (pos <= antGen) {
    delSum <- 0
    while (mark[pos] < 1) {
      delSum <- delSum + y[, pos]
      pos <- pos + 1
    }
    ant[count] <- pos - oldPos
    sum[, count] <- delSum + y[, pos]
    oldPos <- pos
    pos <- pos + 1
    count <- count + 1
  }
  list(Nr = ant, Sum = sum)
}

multiPCFcompact <- function(nr,sum,gamma) {
  ## nr,sum : numbers and sums for one analysis unit,
  ## typically one chromosomal arm. Samples assumed to be in rows.
  ## gamma: penalty for discontinuities
  N <- length(nr)
  nSamples <- nrow(sum)
  ## initialisations
  yhat <- rep(0,N*nSamples);
  dim(yhat) <- c(nSamples,N)
  bestCost <- rep(0,N)
  bestSplit <- rep(0,N+1)
  bestAver <- rep(0,N*nSamples)
  dim(bestAver) <- c(nSamples,N)
  Sum <- rep(0,N*nSamples)
  dim(Sum) <- c(nSamples,N)
  Nevner <- rep(0,N*nSamples)
  dim(Nevner) <- c(nSamples,N)
  eachCost <- rep(0,N*nSamples)
  dim(eachCost) <- c(nSamples,N)
  Cost <- rep(0,N)
  ## Filling of first elements
  Sum[ ,1]<-sum[,1]
  Nevner[,1]<-nr[1]
  bestSplit[1]<-0
  bestAver[,1] <- sum[,1]/nr[1]
  helper <- rep(1, nSamples)
  bestCost[1]<-helper%*%(-Sum[,1]*bestAver[,1])
  lengde <- rep(0,N)
  
  ## Solving for gradually longer arrays. Sum accumulates
  ## error values for righthand plateau downward from n;
  ## this error added to gamma and the stored cost in bestCost
  ## give the total cost. Cost stores the total cost for breaks
  ## at any position below n, and which.min finds the position
  ## with lowest cost (best split). Aver is the average of the
  ## righthand plateau.
  for (n in 2:N) {
    Sum[,1:n] <- Sum[ ,1:n]+sum[,n]
    Nevner[,1:n] <- Nevner[,1:n]+nr[n]
    eachCost[,1:n] <- -(Sum[ ,1:n]^2)/Nevner[ ,1:n]
    Cost[1:n] <- helper %*% eachCost[, 1:n]
    Cost[2:n] <- Cost[2:n]+bestCost[1:(n-1)]+gamma
    Pos <- which.min(Cost[1:n])
    cost <- Cost[Pos]
    aver <- Sum[ ,Pos]/Nevner[,Pos]
    bestCost[n] <- cost
    bestAver[ ,n] <- aver
    bestSplit[n] <- Pos-1
    
  }
  ## The final solution is found iteratively from the sequence
  ## of split positions stored in bestSplit and the averages
  ## for each plateau stored in bestAver
  n <- N
  antInt <- 0
  while (n > 0) {
    yhat[ ,(bestSplit[n]+1):n] <- bestAver[ ,n]
    antInt <- antInt+1
    lengde[antInt] <- sum(nr[(bestSplit[n]+1):n])
    n <- bestSplit[n]
  }
  lengdeRev <- lengde[antInt:1]
  init <- rep(0,antInt)
  init[1]<-1
  if(antInt>=2){
    for(k in 2:antInt){
      init[k]<-init[k-1]+lengdeRev[k-1]
    }
  }
  
  n <- N
  verdi <- rep(0,antInt*nSamples)
  dim(verdi) <- c(nSamples,antInt)
  bestSplit[n+1] <- n
  antall <- antInt
  while (n > 0) {
    verdi[ ,antall] <- bestAver[ ,n]
    n <- bestSplit[n]
    antall <- antall-1
  }
  
  list(Lengde = lengdeRev, sta = init, mean = verdi, nIntervals=antInt)
}