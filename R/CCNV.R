CCNV <- function(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, detail.regions){
  annotate = TRUE
  if(is.null(detail.regions)){
    annotate = FALSE
  }
  
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
    if (ArrayType != "EPIC"){
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    }
    if (ArrayType == "EPIC"){
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
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
    genes <- as.data.frame(detail_regions$name)
    
    
    final_bin_data <-  na.omit(final_bin_data)
    df_gene_regions$bin <- final_bin_data$bin_nr[findInterval(df_gene_regions$pos, final_bin_data$V1)]
    df_gene_regions <- as.data.frame(na.omit(as.matrix(df_gene_regions)))
    df_gene_regions <- separate(df_gene_regions, genes, into = c("genes1", "genes2"), sep = " (?=[^ ]+$)")
    
    #put together the genes and the corresponding bin data
    geneBins <- left_join(genes, df_gene_regions ,join_by("detail_regions$name" == "genes1"))
    
    #remove duplicates and only keep the bin with the highest number of CpG sites present
    geneBins <- geneBins %>% group_by(geneBins$`detail_regions$name`) %>% dplyr::count(geneBins$`detail_regions$name`, geneBins$bin)
    colnames(geneBins) <- c("genes", "bins", "n")
    geneBins <- geneBins[order(geneBins$n, decreasing  = TRUE),]
    geneBins <- na.omit(geneBins[-which(duplicated(geneBins$genes)),])
    
    geneBins <- left_join(geneBins, final_bin_data, join_by("bins" == "bin_nr"))
    geneBins <- geneBins[!((abs(geneBins$V2)) <0.15),]
    
    cumCNV <- ggplot() +
      geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
      ylim(-2, 2) + 
      scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
      geom_vline(xintercept = genome_chr, colour = "grey") + 
      geom_vline(xintercept = genome_centr,linetype="dotted", colour = "grey")  + 
      geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
      labs( y = "Intensity") + 
      geom_point(data = geneBins, aes(x = geneBins$V1, y = geneBins$V2)) + 
      geom_label_repel(data = geneBins, aes(x = geneBins$V1, y = geneBins$V2), label = geneBins$genes, label.size = 0.2, nudge_x = 0.1, nudge_y = 0.1,) + 
      scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
      theme_classic(base_size = 15) + 
      coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
      guides(x = guide_axis(n.dodge = 2))
  }
  
  if (annotate == FALSE){
    cumCNV <- ggplot() +
      geom_point(final_bin_data, mapping = aes(x=V1, y=V2, colour=V2))+ 
      ylim(-2, 2) + 
      scale_color_gradient2(low = colour.loss ,  mid = "white",  high = colour.amplification,  midpoint = 0, breaks = b, guide = "colourbar" , space = "Lab", name = "Intensity") + 
      geom_vline(xintercept = genome_chr, colour = "grey") + 
      geom_vline(xintercept = genome_centr,linetype="dotted", colour = "grey")  + 
      geom_path(data = final_seg_data, aes(x = V1, y = V2)) + 
      labs(y = "Intensity") +
      scale_x_continuous(limits = c(10001, 2881033286), name = "Chromosome", breaks = axis_break, labels = axis_label ) +
      theme_classic(base_size = 15) + 
      coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
      guides(x = guide_axis(n.dodge = 2))
  }
  return(cumCNV)
}