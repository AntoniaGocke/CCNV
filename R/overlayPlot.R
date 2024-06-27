#' Generates the single sample segmentation aberration frequency plot
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param segmentation_data a dataframe of the segmentation results of the single sample segmentation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#'
#' @return returns an overlay plot of all segmentations resulting from the single sample segmentation algorithm
overlayPlot <- function(mSetsAnno, segmentation_data, colour.amplification, colour.loss , array_type) {
  if (array_type == "mouse") {
    #genome data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    
    #calculate position on genome
    seg_start <- c()
    seg_end <- c()
    foreach(i=1:length(segmentation_data$chromosome)) %do%
      {
        seg_start[i] <- segmentation_data$start[i] + addChr[segmentation_data$chromosome[i]]
        seg_end[i] <- segmentation_data$end[i] + addChr[segmentation_data$chromosome[i]] 
        if ((segmentation_data$end[i] - segmentation_data$start[i]) <5000000) {
          segmentation_data$end[i] <- segmentation_data$start[i] + 5000000
        }
      }
    
    segmentation_data$start <- seg_start
    segmentation_data$end <- seg_end
    sample_no <- as.numeric(length(unique(segmentation_data$sample)))
    
    axis_break <- genome_chr
    axis_label <- c(1:length(genome_chr))
    
    segcount <- c(1:length(segmentation_data$segmean))
    segmentation_data$col <- "amplification"
    for (i in segcount) {
      if (segmentation_data$segmean[i] < 0)
      {
        segmentation_data$col[i] <- "deletion"
      }
    }
    
    overlayPlot <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data = segmentation_data,  aes(xmin = start, xmax = end, ymax = segmean, ymin = 0 , fill = col))  +
      ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") +
      ggplot2::ylim(-2.5, 2.5) + 
      ggplot2::labs(y = "Intensity") +
      ggplot2::geom_hline( yintercept = 0, colour ="darkgrey") +
      ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
      ggplot2::theme_classic(base_size = 15) +
      ggplot2::guides(x = guide_axis(angle = 40)) +
      ggplot2::scale_fill_discrete(guide="none") +
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
      ggplot2::guides(x = guide_axis(n.dodge = 2))+
      ggplot2::scale_fill_manual(values = alpha(c("amplification" = colour.amplification,"deletion" = colour.loss), 1/sample_no), guide="none")
    print(overlayPlot)
  }
  else {
    #genome data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    chr_centr <- mSetsAnno$anno_targets@genome$pq
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    genome_centr <- chr_centr + addChr
    
    #calculate position on genome
    seg_start <- c()
    seg_end <- c()
    foreach(i=1:length(segmentation_data$chromosome)) %do%
      {
        seg_start[i] <- segmentation_data$start[i] + addChr[segmentation_data$chromosome[i]]
        seg_end[i] <- segmentation_data$end[i] + addChr[segmentation_data$chromosome[i]] 
        if ((segmentation_data$end[i] - segmentation_data$start[i]) <5000000) {
          segmentation_data$end[i] <- segmentation_data$start[i] + 5000000
        }
      }
    
    segmentation_data$start <- seg_start
    segmentation_data$end <- seg_end
    sample_no <- as.numeric(length(unique(segmentation_data$sample)))
    
    axis_break <- genome_centr
    axis_label <- c(1:length(genome_centr))
    
    segcount <- c(1:length(segmentation_data$segmean))
    segmentation_data$col <- "amplification"
    for (i in segcount) {
      if (segmentation_data$segmean[i] < 0)
      {
        segmentation_data$col[i] <- "deletion"
      }
    }
    
    overlayPlot <- ggplot2::ggplot() + 
      ggplot2::geom_rect(data = segmentation_data,  aes(xmin = start, xmax = end, ymax = segmean, ymin = 0 , fill = col))  +
      ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") +
      ggplot2::ylim(-2.5, 2.5) + 
      ggplot2::labs(y = "Intensity") +
      ggplot2::geom_vline(xintercept = genome_centr, linetype="dotted", colour = "grey") +
      ggplot2::geom_hline( yintercept = 0, colour ="darkgrey") +
      ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
      ggplot2::theme_classic(base_size = 15) +
      ggplot2::guides(x = guide_axis(angle = 40)) +
      ggplot2::scale_fill_discrete(guide="none") +
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
      ggplot2::guides(x = guide_axis(n.dodge = 2))+
      ggplot2::scale_fill_manual(values = alpha(c("amplification" = colour.amplification,"deletion" = colour.loss), 1/sample_no), guide="none")
    
  }

  return(overlayPlot)
  
}