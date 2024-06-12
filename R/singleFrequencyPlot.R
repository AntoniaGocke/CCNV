#' Generates the single sample segmentation aberration frequency plot
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param segmentation_data a dataframe of the segmentation results of the single sample segmentation.
#' @param target_ratios a dataframe of the intesity of all bins for each sample.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param thresh a float specifying the threshold when a segment is called as an aberration
#'
#' @return returns the frequency plot of the samples that were segmented together using the single sample segmentation algorithm
singleFrequencyPlot <- function(mSetsAnno, segmentation_data, colour.amplification, colour.loss, thresh, array_type){
  if (array_type == "mouse") {
    #genome data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    
    axis_break <- genome_chr
    axis_label <- c(1:length(genome_chr))
    
    cnFreq_seg <- as.data.frame(GenVisR::cnFreq(segmentation_data, CN_low_cutoff = -thresh,  CN_high_cutoff = thresh, out = "data", genome = "mm10"))
    cnFreq_seg$data.chromosome <- gsub("chr", "", cnFreq_seg$data.chromosome)
    cnFreq_seg$data.chromosome <- as.numeric(cnFreq_seg$data.chromosome)
    
    data.start = c()
    data.end = c()
    counterCnFreq <- 1:length(cnFreq_seg$data.chromosome)
    for(i in counterCnFreq){
      data.start[i] <- cnFreq_seg$data.start[i] + addChr[cnFreq_seg$data.chromosome[i]]
      data.end[i] <- cnFreq_seg$data.end[i] + addChr[cnFreq_seg$data.chromosome[i]]
    }
    
    cnFreq_seg$data.start <- as.numeric(data.start)
    cnFreq_seg$data.end <- as.numeric(data.end)
    cnFreq_seg$data.lossProportion <- cnFreq_seg$data.lossProportion*-1
    
    #duplicate rows/segment with gains and loss
    df.expanded <- cnFreq_seg[rep(row.names(cnFreq_seg), ifelse(cnFreq_seg$data.gainProportion != 0 & cnFreq_seg$data.lossProportion != 0 ,2 ,1)),]
    
    #now every row is one rectangle of the final plot, extra rows have .1 in their name
    df.expanded$new_rownames <- rownames(df.expanded)
    
    #define the heights of these rectangles
    df.expanded$Y_MIN <- ifelse(df.expanded$data.gainProportion == 0 | df.expanded$data.lossProportion == 0, df.expanded$data.lossProportion, ifelse(grepl("\\.1", df.expanded$new_rownames), df.expanded$data.lossProportion, 0))
    
    df.expanded$Y_MAX <- ifelse(df.expanded$data.gainProportion == 0 | df.expanded$data.lossProportion == 0, df.expanded$data.gainProportion, ifelse(grepl("\\.1", df.expanded$new_rownames), 0, df.expanded$data.gainProportion))
    
    #decribe segment properties
    df.expanded$type <- ifelse(df.expanded$data.gainProportion == 0 & df.expanded$data.lossProportion == 0, "neutral", ifelse(df.expanded$Y_MAX == 0, "loss", "gain"))
    
    #separate into two dataframes for filtering
    gains <- df.expanded[which(df.expanded$type == "gain"),]
    loss <- df.expanded[which(df.expanded$type == "loss"),]
    
    gains_filtered_for_plot <- as.data.frame(cbind(gains$data.start, gains$data.end, gains$data.gainProportion, gains$type))
    loss_filtered_for_plot <- as.data.frame(cbind(loss$data.start, loss$data.end, loss$data.lossProportion, loss$type))
    
    #names(gains_filtered_for_plot) <- c("start", "end", "count", "type")
    names(loss_filtered_for_plot) <- names(gains_filtered_for_plot)
    
    to_plot <- as.data.frame(rbind(gains_filtered_for_plot, loss_filtered_for_plot ))
    names(to_plot) <- c("start", "end", "count", "type")
    
    y_axis_break <-c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
    y_axis_label <- c("100","75", "50","25", "0","25", "50","75", "100" )
    
    to_plot$start <- as.numeric(to_plot$start)
    to_plot$end <- as.numeric(to_plot$end)
    for(i in 1:length(to_plot$end)){
      if ((to_plot$end[i] - to_plot$start[i]) <5000000) {
        to_plot$end[i] <- to_plot$start[i] + 5000000
      }
    }
    
    to_plot$count <- as.numeric(to_plot$count)
    
    singleFreqPlot <- ggplot2::ggplot() + geom_rect(data = to_plot,  aes(xmin = start, xmax = end, ymax = count, ymin = 0 , fill = type))  +
      ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
      ggplot2::geom_hline(yintercept = 0, colour ="darkgrey") + 
      ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
      ggplot2::theme_classic(base_size = 12) + 
      ggplot2::guides(x = guide_axis(angle = 40)) +
      ggplot2::scale_y_continuous(name = "Occurence [%]", breaks = y_axis_break, labels = y_axis_label, limits = c(-1, 1))+
      ggplot2::scale_fill_discrete(guide="none") + 
      ggplot2::scale_fill_manual(values = c("gain" = colour.amplification,"loss" = colour.loss), guide="none") + 
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
      ggplot2::guides(x = guide_axis(n.dodge = 2))
  }
else
  {
    #genome data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    chr_centr <- mSetsAnno$anno_targets@genome$pq
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    genome_centr <- chr_centr + addChr
    
    axis_break <- genome_centr
    axis_label <- c(1:22)
    
    cnFreq_seg <- as.data.frame(GenVisR::cnFreq(segmentation_data, CN_low_cutoff = -thresh,  CN_high_cutoff = thresh, out = "data"))
    cnFreq_seg$data.chromosome <- gsub("chr", "", cnFreq_seg$data.chromosome)
    cnFreq_seg$data.chromosome <- as.numeric(cnFreq_seg$data.chromosome)
    
    data.start = c()
    data.end = c()
    counterCnFreq <- 1:length(cnFreq_seg$data.chromosome)
    for(i in counterCnFreq){
      data.start[i] <- cnFreq_seg$data.start[i] + addChr[cnFreq_seg$data.chromosome[i]]
      data.end[i] <- cnFreq_seg$data.end[i] + addChr[cnFreq_seg$data.chromosome[i]]
    }
    
    cnFreq_seg$data.start <- as.numeric(data.start)
    cnFreq_seg$data.end <- as.numeric(data.end)
    cnFreq_seg$data.lossProportion <- cnFreq_seg$data.lossProportion*-1
    
    #duplicate rows/segment with gains and loss
    df.expanded <- cnFreq_seg[rep(row.names(cnFreq_seg), ifelse(cnFreq_seg$data.gainProportion != 0 & cnFreq_seg$data.lossProportion != 0 ,2 ,1)),]
    
    #now every row is one rectangle of the final plot, extra rows have .1 in their name
    df.expanded$new_rownames <- rownames(df.expanded)
    
    #define the heights of these rectangles
    df.expanded$Y_MIN <- ifelse(df.expanded$data.gainProportion == 0 | df.expanded$data.lossProportion == 0, df.expanded$data.lossProportion, ifelse(grepl("\\.1", df.expanded$new_rownames), df.expanded$data.lossProportion, 0))
    
    df.expanded$Y_MAX <- ifelse(df.expanded$data.gainProportion == 0 | df.expanded$data.lossProportion == 0, df.expanded$data.gainProportion, ifelse(grepl("\\.1", df.expanded$new_rownames), 0, df.expanded$data.gainProportion))
    
    #decribe segment properties
    df.expanded$type <- ifelse(df.expanded$data.gainProportion == 0 & df.expanded$data.lossProportion == 0, "neutral", ifelse(df.expanded$Y_MAX == 0, "loss", "gain"))
    
    #separate into two dataframes for filtering
    gains <- df.expanded[which(df.expanded$type == "gain"),]
    loss <- df.expanded[which(df.expanded$type == "loss"),]
    
    gains_filtered_for_plot <- as.data.frame(cbind(gains$data.start, gains$data.end, gains$data.gainProportion, gains$type))
    loss_filtered_for_plot <- as.data.frame(cbind(loss$data.start, loss$data.end, loss$data.lossProportion, loss$type))
    
    #names(gains_filtered_for_plot) <- c("start", "end", "count", "type")
    names(loss_filtered_for_plot) <- names(gains_filtered_for_plot)
    
    to_plot <- as.data.frame(rbind(gains_filtered_for_plot, loss_filtered_for_plot ))
    names(to_plot) <- c("start", "end", "count", "type")
    
    y_axis_break <-c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
    y_axis_label <- c("100","75", "50","25", "0","25", "50","75", "100" )
    
    to_plot$start <- as.numeric(to_plot$start)
    to_plot$end <- as.numeric(to_plot$end)
    for(i in 1:length(to_plot$end)){
      if ((to_plot$end[i] - to_plot$start[i]) <5000000) {
        to_plot$end[i] <- to_plot$start[i] + 5000000
      }
    }
    
    to_plot$count <- as.numeric(to_plot$count)
    
    singleFreqPlot <- ggplot2::ggplot() + geom_rect(data = to_plot,  aes(xmin = start, xmax = end, ymax = count, ymin = 0 , fill = type))  +
      ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
      ggplot2::geom_vline(xintercept = genome_centr, linetype="dotted", colour = "grey") + 
      ggplot2::geom_hline(yintercept = 0, colour ="darkgrey") + 
      ggplot2::scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
      ggplot2::theme_classic(base_size = 12) + 
      ggplot2::guides(x = guide_axis(angle = 40)) +
      ggplot2::scale_y_continuous(name = "Occurence [%]", breaks = y_axis_break, labels = y_axis_label, limits = c(-1, 1))+
      ggplot2::scale_fill_discrete(guide="none") + 
      ggplot2::scale_fill_manual(values = c("gain" = colour.amplification,"loss" = colour.loss), guide="none") + 
      ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
      ggplot2::guides(x = guide_axis(n.dodge = 2))
  }
  
  return(singleFreqPlot)
}