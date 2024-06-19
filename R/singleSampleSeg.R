#' Segments the data using the conumee package and visualizes DNA methylation data and generative cumulative plots
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param showPlot boolean that determines, if the plots will be displayed
#' 
#' @return Nothing. Will print the figures to the default plotting terminal.
singleSampleSeg<- function(mSetsAnno, thresh, colour.amplification, colour.loss, array_type, showPlot){
  
  #load bin and segment each sample in conumee
  #load bin and segment each sample in conumee
  foreach(i = 1:ncol(mSetsAnno$target_mset_loaded@intensity)) %do%
    {
      if(i==1) {
        x <- conumee::CNV.bin(conumee::CNV.fit( query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets))
        start <- Sys.time()
        x <- conumee::CNV.segment(x)
        end <- Sys.time()
        execution_time_single <- as.numeric(as.POSIXct(start,origin = "1970-01-01")) - as.numeric(as.POSIXct(end,origin = "1970-01-01"))
        segmentation_data <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i]), x@fit[["noise"]]))
      }
      else {
        x <- conumee::CNV.bin(conumee::CNV.fit( query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets))
        start <- Sys.time()
        x <- conumee::CNV.segment(x)
        end <- Sys.time()
        execution_time_singleSeg <- as.numeric(as.POSIXct(start,origin = "1970-01-01")) - as.numeric(as.POSIXct(end,origin = "1970-01-01"))
        execution_time_single <- execution_time_single + execution_time_singleSeg
        target_segmentation <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i]), x@fit[["noise"]]))
        names(target_segmentation) <- names(segmentation_data)
        segmentation_data <- rbind(segmentation_data, target_segmentation)
      }
    }
  print(paste("Runtime single:", execution_time_single))
  
  names(segmentation_data) <- c("chromosome", "start", "end","segmean", "sample", "noise")
  
  CNV_dynamic_cutoff <- thresh
  
  segmentation_data$chromosome <- gsub("chr","",segmentation_data$chromosome)
  segmentation_data$chromosome <- as.numeric(segmentation_data$chromosome)
  sample_no <- as.numeric(length(unique(segmentation_data$sample)))
  segmentation_data$start <- as.numeric(segmentation_data$start)
  segmentation_data$end <- as.numeric(segmentation_data$end)
  segmentation_data$segmean <- as.numeric(segmentation_data$segmean)
  segmentation_data$noise <- as.numeric(segmentation_data$noise)
  segmentation_data2 <- as.data.frame(segmentation_data)
  
  overlayPlot <- overlayPlot(mSetsAnno, segmentation_data, colour.amplification, colour.loss, array_type)
  
  #redefine segmean as gain (1), loss (-1), or neutral (0) in new segmentationdata dataframe
  segmentation_data2$segmean_OLD <- segmentation_data2$segmean
  foreach::foreach(j = 1:nrow(segmentation_data2)) %do% {
    if (segmentation_data2$segmean_OLD[j] > (segmentation_data2$noise[j] * CNV_dynamic_cutoff)) {
      segmentation_data2$segmean[j] <- 1
    } else if (segmentation_data2$segmean_OLD[j] < -(segmentation_data2$noise[j] * CNV_dynamic_cutoff)) {
      segmentation_data2$segmean[j] <- -1
    } else {
      segmentation_data2$segmean[j] <- 0
    }
     
  }
  
  singleFreqPlot <- singleFrequencyPlot(mSetsAnno, segmentation_data2, colour.amplification, colour.loss, thresh, array_type)
  
  if (showPlot == "TRUE"){
    #draw plots
    suppressWarnings(print(overlayPlot))
    suppressWarnings(print(singleFreqPlot))  
  }
  return(segmentation_data)
}