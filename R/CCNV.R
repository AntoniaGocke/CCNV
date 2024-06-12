#' Generates the cumulative copy number plot
#' @param mSetsAnno A list of the RGSet of the target data, the control data and the annotation data
#' @param seg_mpcf a dataframe of the segmentation results of the multi sample segmentation.
#' @param target_ratios a dataframe of the intesity of all bins for each sample.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#'
#' @return returns the plot of the cumulative copy number variation plot

CCNV <- function(mSetsAnno, seg_mpcf, target_ratios, ArrayType, colour.amplification, colour.loss, detail.regions, array_type){
  annotate = TRUE
  if(is.null(detail.regions)){
    annotate = FALSE
  }
  
  if (array_type == "mouse") {
    bin_data_samples <- target_ratios[,-(which(names(target_ratios) %in% c("Chrom", "Median.bp")))]
    
    #here anno longer due to na.omit
    bin_data_all <- as.data.frame(cbind(as.numeric(mSetsAnno$anno_targets@bins@ranges@start), as.numeric(rowMeans(bin_data_samples)), target_ratios$Chrom))
    colnames(bin_data_all)<- c("bin_start", "bin_value", "names")
    
    #genome data
    chr_name <- mSetsAnno$anno_targets@genome$chr
    chr_size <- mSetsAnno$anno_targets@genome$size
    genome_chr <- cumsum(as.numeric(chr_size))
    
    #just for getting pq - lines
    addChr <- c(0,genome_chr)
    addChr <- addChr[-length(addChr)]
    
    #calculate position on genome
    start = c()
    foreach(i=1:length(bin_data_all$bin_start)) %do%
      {
        start[i] <- bin_data_all$bin_start[i] + addChr[bin_data_all$names[i]]
      }
    
    # bin data consisting out of the starting point of the bin and the value
    final_bin_data <- as.data.frame(cbind(as.numeric(start), bin_data_all$bin_value))
    final_bin_data$bin_nr <- rownames(final_bin_data)
    
    #segment data
    seg_data_samples <- seg_mpcf[,-(which(names(seg_mpcf) %in% c("chrom", "arm", "start.pos", "end.pos", "n.probes")))]
    #mean seg value
    segSta <- as.data.frame(cbind(seg_mpcf$start.pos, rowMeans(seg_data_samples), as.numeric(seg_mpcf$chrom)))
    segEnd <- as.data.frame(cbind(seg_mpcf$end.pos, rowMeans(seg_data_samples), as.numeric(seg_mpcf$chrom)))
    names(segEnd) <- names(segSta)
    segment_data <- as.data.frame(rbind(segSta, segEnd))
    colnames(segment_data) <- c("seg_lim", "seg_val", "chr")
    
    #calculate position on genome
    seg_start = c()
    foreach(i=1:length(segment_data$chr)) %do%
      {
        seg_start[i] <- segment_data$seg_lim[i] + addChr[segment_data$chr[i]]
      }
    
    final_seg_data <- as.data.frame(cbind(as.numeric(seg_start), as.numeric(segment_data$seg_val)))
    final_seg_data <-final_seg_data[order(final_seg_data$V1),]
    
    #x-axis
    axis_break <- genome_chr
    axis_label <- c(1:length(genome_chr))
    b <- c(-1, -0.5, 0, 0.5, 1)
    
    
    if (annotate == TRUE) {
      
      print("Annotation function is not yet accessible for mouse data")
      
      cumCNV <- ggplot2::ggplot() +
        ggplot2::geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
        ggplot2::ylim(-2.5, 2.5) + 
        ggplot2::scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
        ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
        ggplot2::geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
        ggplot2::labs(y = "Intensity") +
        ggplot2::scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
        ggplot2::theme_classic(base_size = 15) + 
        ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
        ggplot2::guides(x = guide_axis(n.dodge = 2))
    }
    
    if (annotate == FALSE){
      cumCNV <- ggplot2::ggplot() +
        ggplot2::geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
        ggplot2::ylim(-2.5, 2.5) + 
        ggplot2::scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
        ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
        ggplot2::geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
        ggplot2::labs(y = "Intensity") +
        ggplot2::scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
        ggplot2::theme_classic(base_size = 15) + 
        ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
        ggplot2::guides(x = guide_axis(n.dodge = 2))
    }
  }
  else {
    bin_data_samples <- target_ratios[,-(which(names(target_ratios) %in% c("Chrom", "Median.bp")))]
    bin_data_all <- as.data.frame(cbind(as.numeric(mSetsAnno$anno_targets@bins@ranges@start), as.numeric(rowMeans(bin_data_samples)), target_ratios$Chrom))
    colnames(bin_data_all)<- c("bin_start", "bin_value", "names")
    
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
    start = c()
    foreach(i=1:length(bin_data_all$bin_start)) %do%
      {
        start[i] <- bin_data_all$bin_start[i] + addChr[bin_data_all$names[i]]
      }
    
    # bin data consisting out of the starting point of the bin and the value
    final_bin_data <- as.data.frame(cbind(as.numeric(start), bin_data_all$bin_value))
    final_bin_data$bin_nr <- rownames(final_bin_data)
    
    #segment data
    seg_data_samples <- seg_mpcf[,-(which(names(seg_mpcf) %in% c("chrom", "arm", "start.pos", "end.pos", "n.probes")))]
    #mean seg value
    segSta <- as.data.frame(cbind(seg_mpcf$start.pos, rowMeans(seg_data_samples), as.numeric(seg_mpcf$chrom)))
    segEnd <- as.data.frame(cbind(seg_mpcf$end.pos, rowMeans(seg_data_samples), as.numeric(seg_mpcf$chrom)))
    names(segEnd) <- names(segSta)
    segment_data <- as.data.frame(rbind(segSta, segEnd))
    colnames(segment_data) <- c("seg_lim", "seg_val", "chr")
    
    #calculate position on genome
    seg_start = c()
    foreach(i=1:length(segment_data$chr)) %do%
      {
        seg_start[i] <- segment_data$seg_lim[i] + addChr[segment_data$chr[i]]
      }
    
    final_seg_data <- as.data.frame(cbind(as.numeric(seg_start), as.numeric(segment_data$seg_val)))
    final_seg_data <-final_seg_data[order(final_seg_data$V1),]
    
    #x-axis
    axis_break <- genome_centr
    axis_label <- c(1:22)
    b <- c(-1, -0.5, 0, 0.5, 1)
    
    
    if (annotate == TRUE) {
      if (ArrayType != "450k"){
        anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      }
      if (ArrayType == "EPIC"){
        anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      }
      if (ArrayType == "EPICv2"){
        print("Annotation function is not yet accessible for EPICv2 Arrays")
        cumCNV <- ggplot2::ggplot() +
          ggplot2::geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
          ggplot2::ylim(-2.5, 2.5) + 
          ggplot2::scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
          ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
          ggplot2::geom_vline(xintercept = genome_centr,linetype="dotted", colour = "grey")  + 
          ggplot2::geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
          ggplot2::labs(y = "Intensity") +
          ggplot2::scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
          ggplot2::theme_classic(base_size = 15) + 
          ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
          ggplot2::guides(x = guide_axis(n.dodge = 2))
      }
      
      #detail data -> here used to be anno _EPIC
      gene_regions <- as.data.frame(cbind(anno@listData$UCSC_RefGene_Name, anno@rownames, anno@listData$chr, anno@listData$pos))
      colnames(gene_regions) <- c("genes", "CpG_sites", "chr", "pos")
      gene_regions$genes <- gsub(";"," ",gene_regions$genes)
      gene_regions$chr <- gsub("chr", "", gene_regions$chr)
      df_gene_regions <- gene_regions[!(gene_regions$chr=="X" | gene_regions$chr=="Y"),]
      df_gene_regions$chr <- as.numeric(df_gene_regions$chr)
      df_gene_regions$pos <- as.numeric(df_gene_regions$pos)
      
      #calculate position on genome
      cpg_pos = c()
      counter3 <- 1:length(df_gene_regions$chr)
      for(i in counter3){
        cpg_pos[i] <- df_gene_regions$pos[i] + addChr[df_gene_regions$chr[i]]
      }
      df_gene_regions$pos <- as.numeric(cpg_pos)
      genes <- as.data.frame(detail.regions)
      
      final_bin_data <-  na.omit(final_bin_data)
      df_gene_regions$bin <- final_bin_data$bin_nr[findInterval(df_gene_regions$pos, final_bin_data$V1)]
      df_gene_regions <- as.data.frame(na.omit(as.matrix(df_gene_regions)))
      df_gene_regions <- tidyr::separate(df_gene_regions, genes, into = c("genes1", "genes2"), sep = " (?=[^ ]+$)")
      
      #put together the genes and the corresponding bin data
      geneBins <- as.data.frame(dplyr::left_join(genes, df_gene_regions, by = dplyr::join_by("detail.regions" == "genes1")))
      
      #remove duplicates and only keep the bin with the highest number of CpG sites present
      geneBins <- geneBins %>% dplyr::group_by(geneBins$bins) %>% dplyr::count(geneBins$detail.regions, geneBins$bin)
      colnames(geneBins) <- c("genes", "bins", "n")
      geneBins <- geneBins[order(geneBins$n, decreasing  = TRUE),]
      geneBins <- na.omit(geneBins[-which(duplicated(geneBins$genes)),])
      
      geneBins <- dplyr::left_join(geneBins, final_bin_data, by= dplyr::join_by("bins" == "bin_nr"))
      geneBins <- geneBins[!((abs(geneBins$V2)) <0.15),]
      
      cumCNV <- ggplot2::ggplot() +
        ggplot2::geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
        ggplot2::ylim(-2.5, 2.5) + 
        ggplot2::scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
        ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
        ggplot2::geom_vline(xintercept = genome_centr,linetype="dotted", colour = "grey")  + 
        ggplot2::geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
        ggplot2::labs( y = "Intensity") + 
        ggplot2::geom_point(data = geneBins, aes(x = V1, y = V2)) + 
        ggrepel::geom_label_repel(data = geneBins, aes(x = V1, y = V2), label = geneBins$genes, label.size = 0.2, nudge_x = 0.1, nudge_y = 0.1,) + 
        ggplot2::scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
        ggplot2::theme_classic(base_size = 15) + 
        ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
        ggplot2::guides(x = guide_axis(n.dodge = 2))
    }
    
    if (annotate == FALSE){
      cumCNV <- ggplot2::ggplot() +
        ggplot2::geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
        ggplot2::ylim(-2.5, 2.5) + 
        ggplot2::scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
        ggplot2::geom_vline(xintercept = genome_chr, colour = "grey") + 
        ggplot2::geom_vline(xintercept = genome_centr,linetype="dotted", colour = "grey")  + 
        ggplot2::geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
        ggplot2::labs(y = "Intensity") +
        ggplot2::scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
        ggplot2::theme_classic(base_size = 15) + 
        ggplot2::coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
        ggplot2::guides(x = guide_axis(n.dodge = 2))
    }
  }
  
  return(cumCNV)
}