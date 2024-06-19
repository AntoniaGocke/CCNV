#' Segments the data using the multipcf algorithm and visualizes DNA methylation data and generative cumulative plots using conumee
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param array_type Type of methylation array
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param showPlot boolean that determines, if the plots will be displayed
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
multiSampleSeg <- function(mSetsAnno, thresh, array_type, colour.amplification, colour.loss, detail.regions, showPlot){
  target_mset_loaded <- mSetsAnno$target_mset_loaded
  control_mset_loaded <- mSetsAnno$control_mset_loaded
  anno_targets <- mSetsAnno$anno_targets
  #load and bin each sample in conumee
  foreach(i=1:ncol(mSetsAnno$target_mset_loaded@intensity)) %do%
    {
      if (i == 1) {
        x <- conumee::CNV.bin(conumee::CNV.fit(query = target_mset_loaded[names(target_mset_loaded@intensity[i])], ref = control_mset_loaded,  anno_targets))
        target_ratios <- cbind(x@anno@bins@ranges@NAMES, x@anno@bins@ranges@start, x@bin$ratio)
      }
      else {
        x <- conumee::CNV.bin(conumee::CNV.fit(query = target_mset_loaded[names(target_mset_loaded@intensity[i])], ref = control_mset_loaded,  anno_targets))
        target_ratios <- as.data.frame(cbind(target_ratios, x@bin$ratio))
      }
    }
    
  names(target_ratios) <- c("Chrom", "Median.bp", names(mSetsAnno$target_mset_loaded@intensity))
  
  target_ratios$Chrom <- gsub("chr","",target_ratios$Chrom)
  target_ratios$Chrom <- substr(target_ratios$Chrom,1,nchar(target_ratios$Chrom)-4)
  target_ratios$Chrom <- sub("\\-", "", target_ratios$Chrom)
  
  target_ratios <- na.omit(as.data.frame(sapply(target_ratios, as.numeric)))
  
  #################### Segmentation #################################
  start <- Sys.time()
  seg_mpcf <- FastMultiPCF(target_ratios, gamma = 5)
  end <- Sys.time()
  execution_time_multi <- as.numeric(as.POSIXct(start,origin = "1970-01-01")) - as.numeric(as.POSIXct(end,origin = "1970-01-01"))
  print(paste("Runtime multi:", execution_time_multi))
  suppressPackageStartupMessages(require(dplyr))
  cumCNV <- CCNV(mSetsAnno, seg_mpcf, target_ratios, array_type, colour.amplification, colour.loss, detail.regions, array_type)
  cumFreq <- cumFreq(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, thresh, array_type)
  
  
  if (showPlot == "TRUE") {
    #draw plots
    suppressWarnings(print(cumCNV))
    suppressWarnings(print(cumFreq))
  }
  
  return(seg_mpcf)
  
}