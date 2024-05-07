multiSampleSeg <- function(mSetsAnno, ArrayType, thresh, colour.amplification, colour.loss, detail.regions){
  
  #load and bin each sample in conumee
  x <- CNV.bin(CNV.fit(query = mSetsAnno$target_mset_loaded, ref = mSetsAnno$control_mset_loaded, mSetsAnno$anno_targets))
  target_ratios <- cbind(x@anno@bins@ranges@NAMES, x@anno@bins@ranges@start, as.data.frame(x@bin$ratio))
  names(target_ratios) <- c("Chrom", "Median.bp", names(mSetsAnno$target_mset_loaded@intensity))
  
  target_ratios$Chrom <- gsub("chr","",target_ratios$Chrom)
  target_ratios$Chrom <- substr(target_ratios$Chrom,1,nchar(target_ratios$Chrom)-4)
  target_ratios$Chrom <- sub("\\-", "", target_ratios$Chrom)
  
  target_ratios <- na.omit(as.data.frame(sapply(target_ratios, as.numeric)))
  
  #################### Segmentation #################################
  seg_mpcf <- updatempcf(target_ratios, gamma = 5)
  
  cumCNV <- CCNV(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, detail.regions)
  cumFreq <- cumFreq(mSetsAnno, seg_mpcf, target_ratios, colour.amplification, colour.loss, thresh)
  
  
  #draw plots
  print(cumCNV)
  print(cumFreq)
  
}