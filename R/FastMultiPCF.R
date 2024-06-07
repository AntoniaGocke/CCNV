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
FastMultiPCF <- function(data, gamma=5, frac1=0.15, frac2=0.15){
  
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
    nSeg <- mpcf$nBins
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
runFastMultiPCF <- function (intensities, gamma, frac1, frac2) {
  bpoints<-findBreakpoints(intensities,frac1,frac2)
  aggr <- aggregateBins(t(intensities), bpoints)
  segs <- calcScoresSegs(aggr$Dist, aggr$Sum, gamma)
  return(list(length = segs$segLen, start0 = segs$sta,
              mean = segs$mean, nBins = segs$nBins))
}


#' Defines a highpass filter with double the specified length
#' @param half_length Parameter for the filter length, which will be 2*half_length
#' @return vector of numeric with the filter elements
define_cutoff <- function(half_length){
    filter <- rep(0, 2*half_length)
    for (k in 1:half_length) {
        filter[k] <- k / half_length
        filter[2 * half_length + 1 - k] <- -k / half_length
    }
    return(filter)
}

#' Applies a given filter over each column of a matrix and
#' computes the respective sums
#' @param x The matrix
#' @param filter the filter (vector of numeric)
#' @param start left-spacing of each row
#' @param end right-spacing of each row
#' @return vector of numeric with the filter elements
apply_filter <- function(x, filter, start, end){
  hpfValue <- rep(0, nrow(x))
    L <- length(filter)/2 # filter size assumed to be even
    for (l in start:(nrow(x) - end)) {
        for (m in 1:ncol(x)) {
            diff <- crossprod(filter, x[l:(l + 2 * L - 1), m])
            hpfValue[l + L - 1] <- hpfValue[l + L - 1] + abs(diff)
        }
    }
    return(hpfValue)
}

#' Takes a vector of values and marks positions at which the respective elements
#' are larger than a given quantile
#' @param bpoints Vector to write the marks to. Smaller elements are marked with
#' a zero and larger elements are marked by one
#' @param hpfValue Vector of numeric. Must have same length as bpoints
#' @param halflength Half the length of the filter used to create hpfValue
#' @param frac Threshold at 1-frac quantile
#' @param start left-spacing of each row
#' @param end right-spacing of each row
#' @return the modified vector of potential breakpoints
mark_points <- function(bpoints, hpfValue, halflength, frac, start, end){
    limit <- quantile(hpfValue, (1 - frac))
    
    for (l in start:(length(hpfValue) - end)) {
        if (hpfValue[l + halflength - 1] > limit) {
          bpoints[l + halflength - 1] <- 1
        }
    }
    return(bpoints)
}

## Highpass filter for heuristically determining potential breakpoints for multiPCF based on implementation in XXX
findBreakpoints <- function(x, frac1, frac2){
    L <- 15
    bpoints <- rep(0, nrow(x))
    hpfValue <- rep(0, nrow(x))
    filter <- define_cutoff(L)
    hpfValue2 <- rep(0, nrow(x))
    filter2 <- define_cutoff(3)
    hpfValue <- apply_filter(x, filter, 1, (2*L-1))
    hpfValue2 <- apply_filter(x, filter2, L-1, L+2)
    bpoints <- mark_points(bpoints, hpfValue, L, frac1, 1, 2*L)
    bpoints <- mark_points(bpoints, hpfValue2, 3, frac2, L-1, L+2)
    bpoints[1:L] <- 1
    for (l in 1:L) {
      bpoints[nrow(x) + 1 - l] <- 1
    }
    return(bpoints)
}


# function that accumulates numbers of observations and sums between potential breakpoints
#delsum ist summe aller bin intensitäten bis zum heuristischen breakpoint
#returns a list nr and sum. nr is a vector, which tracks the number of bins to the last breakpoint, sum is a dataframe, that 
#has the dimensions nof_breakpoints x samples and adds the intensities of each bin for each sample until the breakpoint
aggregateBins <- function(intensities, bpoints) {
  nof_bpoints <- sum(bpoints) #anzahl gesetzte breakpoints
  bDist <- rep(0, nof_bpoints)
  sums <- data.frame(matrix(0, nrow(intensities), nof_bpoints)) #probenanzal x breakpoints
  intSum <- 0
  oldPos <- 0
  counter <- 1
  
  #iterate over each bin. calculate distances and sums between bpoints 
  foreach::foreach(i = 1:ncol(intensities)) %do% {
    if(bpoints[i] != 1) {
      intSum <- intSum + intensities[, i] 
    } else {
      #potential breakpoint
      bDist[counter] <- i - oldPos 
      sums[, counter] <- intSum + intensities[, i]
      oldPos <- i
      counter <- counter + 1
      intSum <- 0
    }
  }
  list(Dist = bDist, Sum = sums)
}


#' Takes a vector of values and marks positions at which the respective elements
#' are larger than a given quantile
#' @param dist Vector with distances between breakpoints
#' @param sum matrix with summed intensities between each breakpoint
#' @param gamma threshold
#' @return the modified vector of potential breakpoints
calcScoresSegs <- function(dist, sum, gamma) {
  N <- length(dist)
  nSamples <- nrow(sum)
  
  # initialise variables and matrices
  bestCost <- rep(0, N)
  Cost <- rep(0, N)  
  segLen <- rep(0, N)
  bestSeg <- rep(0, N+1)
  Aver <- matrix(0, nSamples, N)
  Sum <- matrix(0, nSamples, N)
  Dist <- matrix(0, nSamples, N)
  eachCost <- matrix(0, nSamples, N)
  
  # initialise matrices
  Sum[, 1] <- sum[, 1]
  Dist[, 1] <- dist[1]
  Aver[, 1] <- sum[, 1] / dist[1]
  bestCost[1] <- -sum(sum[, 1] * Aver[, 1])
  
  for (n in 2:N) {
    Sum[,1:n] <- Sum[,1:n] + sum[,n]
    Dist[,1:n] <- Dist[,1:n] + dist[n]
    eachCost[,1:n] <- -(Sum[,1:n]^2) / Dist[,1:n]
    Cost <- colSums(eachCost[, 1:n])
    Cost[2:n] <- Cost[2:n] + bestCost[1:(n-1)] + gamma
    bestCost[n] <- Cost[which.min(Cost[1:n])]
    Aver[,n] <- Sum[,which.min(Cost[1:n])] / Dist[,which.min(Cost[1:n])]
    bestSeg[n] <- which.min(Cost[1:n]) - 1
  }
  
  # Calculate best segments
  n <- N
  antInt <- 0
  while (n > 0) {
    antInt <- antInt + 1
    segLen[antInt] <- sum(dist[(bestSeg[n] + 1):n])
    n <- bestSeg[n]
  }
  #move from right to left
  segLenRev <- rev(segLen[1:antInt])
  SegSta <- cumsum(c(1, segLenRev[-antInt]))
  if (antInt >= 2) {
    for (j in 2:antInt) {
      SegSta[j] <- SegSta[j - 1] + segLenRev[j - 1]
    }
  }
  
  # calculate mean intensities per segment
  segMean <- matrix(0, nSamples, antInt)
  n <- N
  bestSeg[n + 1] <- n
  antall <- antInt
  while (n > 0) {
    segMean[, antall] <- Aver[, n]
    n <- bestSeg[n]
    antall <- antall - 1
  }
  
  # Return results as a list
  list(segLen = segLenRev, sta = SegSta, mean = segMean, nBins = antInt)
}