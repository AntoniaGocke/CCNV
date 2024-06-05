#' Generates the multi sample segmentation aberration frequency plot
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param seg_mpcf a dataframe of the segmentation results of the multi sample segmentation.
#' @param target_ratios a dataframe of the intesity of all bins for each sample.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param thresh a float specifying the threshold when a segment is called as an aberration
#'
#' @return returns the frequency plot of the samples that were segmented togethher using the multi sample segmentation algorithm
#' 
cumFreq <- function(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, thresh){
  chr_name <- mSetsAnno$anno_targets@genome$chr
  chr_size <- mSetsAnno$anno_targets@genome$size
  genome_chr <- cumsum(as.numeric(chr_size))
  chr_centr <- mSetsAnno$anno_targets@genome$pq
  
  #just for getting pq - lines
  addChr <- c(0,genome_chr)
  addChr <- addChr[-length(addChr)]
  genome_centr <- chr_centr + addChr
  
  #defining threshholds
  broad_min_probes <- 300
  focal_min_probes <- 2
  
  #x-axis
  axis_break <- genome_centr
  axis_label <- c(1:22)
  b <- c(-1, -0.5, 0, 0.5, 1)
  
  #segment
  segstart <- seg_mpcf$start.pos
  segend <- seg_mpcf$end.pos
  segnames <- as.numeric(seg_mpcf$chrom)
  
  segment_data <- as.data.frame(cbind(segstart, segend, segnames))
  colnames(segment_data) <- c("seg_start", "seg_end" , "chr")
  
  #calculate position on genome
  seg_start <- c()
  seg_end <- c()
  foreach(i=1:length(segment_data$chr)) %do%
    {
      seg_start[i] <- segment_data$seg_start[i] + addChr[segment_data$chr[i]]
      seg_end[i] <- segment_data$seg_end[i] + addChr[segment_data$chr[i]]
    }
  
  seg_data_samples <- seg_mpcf[,-(which(names(seg_mpcf) %in% c("chrom", "arm", "start.pos", "end.pos", "n.probes")))]
  df_seg_val <- as.data.frame(cbind(seg_start, seg_end,  rowMeans(seg_data_samples)))
  df_abberation_candidates <- as.data.frame(cbind(seg_start, seg_end, seg_mpcf$n.probes, seg_data_samples))
  
  #focal candidates
  focal_candidates <- df_abberation_candidates
  #focal_candidates$seg_end <- focal_candidates$seg_end + 10000000
  for (i in 1:length(focal_candidates$seg_end)) {
    if ((focal_candidates$seg_end[i] - focal_candidates$seg_start[i]) <5000000) {
      focal_candidates$seg_end[i] <- focal_candidates$seg_start[i] + 5000000
    }
  }
  focal_candidates <- focal_candidates[!((abs(focal_candidates$`seg_mpcf$n.probes`)) < focal_min_probes),]
  focal_candidates$mean <- rowMeans(seg_data_samples)
  focal_candidates$col <- "amplification"
  focal_candid_values <- focal_candidates[,-(1:5)]
  focalcount <- c(1:length(focal_candidates$seg_start))
  for (i in focalcount) {
    focal_candidates$count[i] <- max(0, (sum(focal_candid_values[i,] >= thresh)) )
    if (focal_candidates$mean[i] < 0)
    {
      focal_candidates$count[i] <- focal_candidates$count[i] * -1
      focal_candidates$col[i] <- "deletion"
    }
  }
  
  #broad candidates
  broad_candidates <- df_abberation_candidates
  broad_candidates$mean <- rowMeans(seg_data_samples)
  broad_candidates <- broad_candidates[!((abs(broad_candidates$`seg_mpcf$n.probes`)) < broad_min_probes),]
  broad_candid_values <- abs(broad_candidates[,-(1:5)])
  broadcount <- c(1:length(broad_candid_values$mean))
  broad_candidates$col <- "amplification"
  for (i in broadcount) {
    broad_candidates$count[i] <- max(0, (sum(broad_candid_values[i,] >= thresh)) )
    if (broad_candidates$mean[i] < 0)
    {
      broad_candidates$count[i] <- broad_candidates$count[i] * -1
      broad_candidates$col[i] <- "deletion"
    }
  }
  if (nrow(focal_candidates[(focal_candidates$count!=0),]) == 0) {
    candidates <- broad_candidates[(broad_candidates$count!=0),]
  } else if (nrow(broad_candidates[(broad_candidates$count!=0),]) == 0) {
    candidates <- focal_candidates[(focal_candidates$count!=0),]
  } else {
    candidates <- rbind(focal_candidates[(focal_candidates$count!=0),], broad_candidates[(broad_candidates$count!=0),])
  }
  sample_no <- ncol(seg_data_samples)
  b <- c(-1, -0.5, 0, 0.5, 1)
  y_axis_break <-c(-sample_no,-(sample_no/2) + -(sample_no/4), -(sample_no/2), -(sample_no/4),0, sample_no/4, sample_no/2,(sample_no/2) + (sample_no/4), sample_no)
  y_axis_label <- c("100","75", "50","25", "0","25", "50","75", "100" )
  
  cumFreq <- ggplot2::ggplot() + 
    ggplot2::geom_rect(data = candidates,  aes(xmin = seg_start, xmax = seg_end, ymax = count, ymin = 0 , fill = col))  +
    ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") +
    ggplot2::geom_vline(xintercept = genome_centr, linetype="dotted", colour = "grey") +
    ggplot2::geom_hline( yintercept = 0, colour ="darkgrey") +
    ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::guides(x = guide_axis(angle = 40)) +
    ggplot2::scale_y_continuous(name = "Occurrence over all samples [%]", breaks = y_axis_break, labels = y_axis_label, limits = c(-sample_no, sample_no)) +
    ggplot2::scale_fill_discrete(guide="none") +
    ggplot2::scale_fill_manual(values = c("amplification" = colour.amplification,"deletion" = colour.loss), guide="none") +
    ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
    ggplot2::guides(x = guide_axis(n.dodge = 2))
}