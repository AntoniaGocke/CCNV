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
cumFreq <- function(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, thresh, array_type){
  if (array_type == "mouse") {
    #Genome Data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    
    #defining threshholds
    broad_min_probes <- 300
    focal_min_probes <- 2
    
    #x-axis
    axis_break <- genome_chr
    axis_label <- c(1:length(genome_chr))
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
    
    seg_data_samples <- seg_mpcf[,-(which(names(seg_mpcf) %in% c("chrom",  "start.pos", "end.pos", "n.probes")))]
    df_seg_val <- as.data.frame(cbind(seg_start, seg_end,  rowMeans(seg_data_samples)))
    candidates <-  as.data.frame(cbind(seg_start, seg_end, seg_mpcf$n.probes, seg_data_samples))
    
    
    for (i in 1:length(candidates$seg_end)) {
      if ((candidates$seg_end[i] - candidates$seg_start[i]) <5000000) {
        candidates$seg_end[i] <- candidates$seg_end[i] + 5000000
      }
    }
    
    candidates$mean <- rowMeans(seg_data_samples)
    candidates <- candidates[!((abs(candidates$`seg_mpcf$n.probes`)) < min_probes),]
    candid_values <- candidates[,-(which(names(candidates) %in% c("seg_start", "seg_end", "seg_mpcf$n.probes", "mean")))]
    candidates$col <- "amplification"
    candidates$count <- 0
    
    amplifications <- candidates[1,]
    deletions <- amplifications
    count <- c(1:length(candidates$seg_start))
    
    for (i in count) {
      ampl <- max(0, (sum(candid_values[i,] >= thresh)) )
      del <- max(0, (sum(candid_values[i,] <= -thresh)) )
      if (0 < del)
      { candidates$count[i] <- -del
      candidates$col[i] <- "deletion"
      deletions <- rbind(deletions, candidates[i,])
      }
      if (0 < ampl)
      {
        candidates$count[i] <- ampl
        candidates$col[i] <- "amplification"
        amplifications <- rbind(amplifications, candidates[i,])
      }
      
    }
    
    amplifications <- amplifications[-1,]
    deletions <- deletions[-1,]
    
    if(is.null(amplifications) ) {
      candidates <- deletions
    } else if (is.null(deletions)) {
      candidates <- amplifications
    } else {
      candidates <- rbind(amplifications, deletions)
    }
    sample_no <- ncol(seg_data_samples)
    b <- c(-1, -0.5, 0, 0.5, 1)
    y_axis_break <-c(-sample_no,-(sample_no/2) + -(sample_no/4), -(sample_no/2), -(sample_no/4),0, sample_no/4, sample_no/2,(sample_no/2) + (sample_no/4), sample_no)
    y_axis_label <- c("100","75", "50","25", "0","25", "50","75", "100" )
    
    cumFreq <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = candidates,  aes(xmin = seg_start, xmax = seg_end, ymax = count, ymin = 0 , fill = col))  +
      ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") +
      ggplot2::geom_hline( yintercept = 0, colour ="darkgrey") +
      ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
      ggplot2::theme_classic(base_size = 15) +
      ggplot2::guides(x = guide_axis(angle = 40)) +
      ggplot2::scale_y_continuous(name = "Occurence [%]", breaks = y_axis_break, labels = y_axis_label, limits = c(-sample_no, sample_no)) +
      ggplot2::scale_fill_discrete(guide="none") +
      ggplot2::scale_fill_manual(values = c("amplification" = colour.amplification,"deletion" = colour.loss), guide="none") +
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
      ggplot2::guides(x = guide_axis(n.dodge = 2))
  }
  else {
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    chr_centr <- mSetsAnno$anno_targets@genome$pq
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    genome_centr <- chr_centr + addChr
    
    #defining threshholds
    min_probes <- 2
    
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
    
    seg_data_samples <- seg_mpcf[,-(which(names(seg_mpcf) %in% c("chrom", "start.pos", "end.pos", "n.probes")))]
    
    df_seg_val <- as.data.frame(cbind(seg_start, seg_end,  rowMeans(seg_data_samples)))
    
    candidates <-  as.data.frame(cbind(seg_start, seg_end, seg_mpcf$n.probes, seg_data_samples))
    
    
    for (i in 1:length(candidates$seg_end)) {
      if ((candidates$seg_end[i] - candidates$seg_start[i]) <5000000) {
        candidates$seg_end[i] <- candidates$seg_end[i] + 5000000
      }
    }
    
    candidates$mean <- rowMeans(seg_data_samples)
    candidates <- candidates[!((abs(candidates$`seg_mpcf$n.probes`)) < min_probes),]
    candid_values <- candidates[,-(which(names(candidates) %in% c("seg_start", "seg_end", "seg_mpcf$n.probes", "mean")))]
    candidates$col <- "amplification"
    candidates$count <- 0
    
    amplifications <- candidates[1,]
    deletions <- amplifications
    count <- c(1:length(candidates$seg_start))
    
    for (i in count) {
      ampl <- max(0, (sum(candid_values[i,] >= thresh)) )
      del <- max(0, (sum(candid_values[i,] <= -thresh)) )
      if (0 < del)
      { candidates$count[i] <- -del
      candidates$col[i] <- "deletion"
      deletions <- rbind(deletions, candidates[i,])
      }
      if (0 < ampl)
      {
        candidates$count[i] <- ampl
        candidates$col[i] <- "amplification"
        amplifications <- rbind(amplifications, candidates[i,])
      }
      
    }
    
    amplifications <- amplifications[-1,]
    deletions <- deletions[-1,]
    
    if(is.null(amplifications) ) {
      candidates <- deletions
    } else if (is.null(deletions)) {
      candidates <- amplifications
    } else {
      candidates <- rbind(amplifications, deletions)
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
      ggplot2::scale_y_continuous(name = "Occurence [%]", breaks = y_axis_break, labels = y_axis_label, limits = c(-sample_no, sample_no)) +
      ggplot2::scale_fill_discrete(guide="none") +
      ggplot2::scale_fill_manual(values = c("amplification" = colour.amplification,"deletion" = colour.loss), guide="none") +
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
      ggplot2::guides(x = guide_axis(n.dodge = 2))
  }
  return(cumFreq)
}
