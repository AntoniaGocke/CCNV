#' Segments the data using the conumee package and visualizes DNA methylation data and generative cumulative plots
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
singleSampleSeg<- function(mSetsAnno, thresh, colour.amplification, colour.loss){
  
  #load bin and segment each sample in conumee
  foreach(i=1:ncol(mSetsAnno$target_mset_loaded@intensity)) %do%
    {
      if(i==1) {
        x <- conumee::CNV.segment(conumee::CNV.bin(conumee::CNV.fit(query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_set_loaded, anno = mSetsAnno$anno_targets)))
        segmentation_data <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i])))
      }
      else {
        x <- conumee::CNV.segment(conumee::CNV.bin(conumee::CNV.fit(query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_set_loaded, anno = mSetsAnno$anno_targets)))
        target_segmentation <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i])))
        names(target_segmentation) <- names(segmentation_data)
        segmentation_data <- rbind(segmentation_data, target_segmentation)
      }
    }
  
  names(segmentation_data) <- c("chromosome", "start", "end","segmean", "sample")
  
  segmentation_data$chromosome <- gsub("chr","",segmentation_data$chromosome)
  segmentation_data$chromosome <- as.numeric(segmentation_data$chromosome)
  sample_no <- as.numeric(length(unique(segmentation_data$sample)))
  segmentation_data$start <- as.numeric(segmentation_data$start)
  segmentation_data$end <- as.numeric(segmentation_data$end)
  segmentation_data$segmean <- as.numeric(segmentation_data$segmean)
  segmentation_data2 <- as.data.frame(segmentation_data)
  
  overlayPlot <- overlayPlot(mSetsAnno, segmentation_data2, colour.amplification, colour.loss)
  singleFreqPlot <- singleFrequencyPlot(mSetsAnno, segmentation_data2, colour.amplification, colour.loss, thresh)
  
  #draw plots
  print(overlayPlot)
  print(singleFreqPlot)
}