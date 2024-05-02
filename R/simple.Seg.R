simple.Seg <- function(mSetsAnno, ArrayType, thresh, colour.amplification, colour.loss){
    
    start.time <- Sys.time()
    #load and bin each sample in conumee
    x <- CNV.segment(CNV.bin(CNV.fit(query = target_mset_loaded, ref = control_mset_loaded, anno_targets)))
    foreach(i=1:ncol(mSetsAnno$target_mset_loaded@intensity)) %do%
        {
            if (i == 1) {
                x <- CNV.segment(CNV.bin(CNV.fit(query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets)))
                segmentation_data <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i])))
            }
            else {
                x <- CNV.segment(CNV.bin(CNV.fit(query = mSetsAnno$target_mset_loaded[names(mSetsAnno$target_mset_loaded[i])], ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets)))
                target_segmentation <- as.data.frame(cbind(x@seg$summary$chrom, x@seg$summary$loc.start, x@seg$summary$loc.end, x@seg$summary$seg.mean, names(mSetsAnno$target_mset_loaded[i])))
                names(target_segmentation) <- names(segmentation_data)
                segmentation_data <- rbind(segmentation_data, target_segmentation)
            }
        }
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    names(segmentation_data) <- c("chromosome", "start", "end","segmean", "sample")
    
    segmentation_data$chromosome <- gsub("chr","",segmentation_data$chromosome)
    segmentation_data$chromosome <- as.numeric(segmentation_data$chromosome)
    sample_no <- as.numeric(length(unique(segmentation_data$sample)))
    segmentation_data$start <- as.numeric(segmentation_data$start)
    segmentation_data$end <- as.numeric(segmentation_data$end)
    segmentation_data$segmean <- as.numeric(segmentation_data$segmean)
    segmentation_data2 <- as.data.frame(segmentation_data)
    
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
            seg_end[i] <- segmentation_data$end[i] + addChr[segmentation_data$chromosome[i]] + 10000000
        }
    
    segmentation_data$start <- seg_start
    segmentation_data$end <- seg_end
    
    
    axis_break <- genome_centr
    axis_label <- c(1:22)
    
    segcount <- c(1:length(segmentation_data$segmean))
    segmentation_data$col <- "amplification"
    for (i in segcount) {
        if (segmentation_data$segmean[i] < 0)
        {
            segmentation_data$col[i] <- "deletion"
        }
    }
    
    overlayPlot <- ggplot() + 
        geom_rect(data = segmentation_data,  aes(xmin = start, xmax = end, ymax = segmean, ymin = 0 , fill = col))  +
        geom_vline(xintercept = genome_chr, colour = "grey") +
        ylim(-2, 2) + 
        geom_vline(xintercept = genome_centr, linetype="dotted", colour = "grey") +
        geom_hline( yintercept = 0, colour ="darkgrey") +
        scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
        theme_classic(base_size = 15) +
        guides(x = guide_axis(angle = 40)) +
        scale_fill_discrete(guide="none") +
        coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) +
        guides(x = guide_axis(n.dodge = 2))+
        scale_fill_manual(values = alpha(c("amplification" = colour.amplification,"deletion" = colour.loss), 1/sample_no), guide="none")
    
    cnFreq_seg <- as.data.frame(cnFreq(segmentation_data2, CN_low_cutoff = -thresh,  CN_high_cutoff = thresh, out = "data"))
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
    
    b <- c(-1, -0.5, 0, 0.5, 1)
    y_axis_break <-c(-sample_no,-(sample_no/2) + -(sample_no/4), -(sample_no/2), -(sample_no/4),0, sample_no/4, sample_no/2,(sample_no/2) + (sample_no/4), sample_no)
    y_axis_label <- c("100","75", "50","25", "0","25", "50","75", "100" )
    
    to_plot$start <- as.numeric(to_plot$start)
    to_plot$end <- as.numeric(to_plot$end)
    to_plot$end <- to_plot$end + 5000000
    to_plot$count <- as.numeric(to_plot$count)
    
    singleFreqPlot <- ggplot() + geom_rect(data = to_plot,  aes(xmin = start, xmax = end, ymax = count, ymin = 0 , fill = type))  +
        geom_vline(xintercept = genome_chr, colour = "grey") + 
        geom_vline(xintercept = genome_centr, linetype="dotted", colour = "grey") + 
        geom_hline(yintercept = 0, colour ="darkgrey") + 
        scale_x_continuous(name = "Chromosome", breaks = axis_break, labels = axis_label, limits = c(0, 2881033286)) +
        theme_classic(base_size = 12) + 
        guides(x = guide_axis(angle = 40)) +
        scale_y_continuous(name = "Intensity", breaks = y_axis_break, labels = y_axis_label, limits = c(-sample_no, sample_no))+
        scale_fill_discrete(guide="none") + 
        scale_fill_manual(values = c("gain" = colour.amplification,"loss" = colour.loss), guide="none") + 
        coord_cartesian(xlim = c(50001, 2881033286), expand = FALSE) + 
        guides(x = guide_axis(n.dodge = 2))
    
    #draw plots
    print(overlayPlot)
    print(singleFreqPlot)
}