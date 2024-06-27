#' Takes the weighted mean intensity percommom segment for each sample or chromosome and determines whether a chromosome 
#' has am amplification (1), a loss (-1) or no change (0). This is loaded into a dataframe and returned. A similarity clustering 
#' is performed usind the euclidean distance between the samples to reorder them, based on their similarity. Sample that show similar
#' segment boudaries (from the combined segmentation) can be processed segment wise (chrom = FALSE) or chromosome wise (chrom = TRUE).
#' Samples segmented separately will be processed on the choromosomal level.
#' @param segmentation_data dataframe that has been returned from the function cumul.CNV
#' @param thresh a threshold to determine from when on a segment is considered an aberration
#' @param chrom boolean, whether the aberration should be calculated on the segment level (FALSE) or on the chromosomal level (TRUE)
#'
#' @return returns the frequency plot of the samples that were segmented togethher using the multi sample segmentation algorithm
#' @export
#'
get.chromAberrations <- function(segmentation_data, thresh = 0.2, chrom = FALSE){
  require(dplyr)
  if(all(c("chrom", "start.pos", "end.pos", "n.probes") %in% names(segmentation_data)))
  {
    if (chrom == TRUE) {
      segmentation_data$n.bp <- segmentation_data$end.pos - segmentation_data$start.pos
      segmentation_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("start.pos", "end.pos", "n.probes"))]
      sample_names <- names(segmentation_data_filtered)
      sample_names <- sample_names[-1]
      chromosomal_weights <- segmentation_data_filtered %>% group_by(chrom) %>%  mutate(weight = n.bp / sum(n.bp),) 
      
      chromosomal_average_multi <- chromosomal_weights %>% group_by(chrom) %>% summarise(across(sample_names, ~weighted.mean(., w = weight)))
      chromosomal_average_multi <- chromosomal_average_multi[,-which(names(chromosomal_average_multi) %in% c("n.bp", "chrom")) ,drop=FALSE]
      
      chromosomal_average_multi[chromosomal_average_multi > 0 & chromosomal_average_multi <= thresh ] <- 1
      chromosomal_average_multi[chromosomal_average_multi <= -thresh ] <- -1
      chromosomal_average_multi[chromosomal_average_multi != 1 & chromosomal_average_multi != -1 ] <- 0
      
      euclideanDist <- hclust(dist(t(chromosomal_average_multi)))
      ordereuclideanDist <- euclideanDist$order
      chromosomal_average_multi <- chromosomal_average_multi[,ordereuclideanDist]
      output <- cbind(chromosome = c(1:length(t(chromosomal_average_multi)[1,])), chromosomal_average_multi)
      
    } else {
      segmentation_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("chrom", "start.pos", "end.pos", "n.probes"))]
      segment_average_multi <- segmentation_data_filtered
      segment_average_multi[segment_average_multi > 0 & segment_average_multi <= thresh ] <- 1
      segment_average_multi[segment_average_multi <= -thresh ] <- -1
      segment_average_multi[segment_average_multi != 1 & segment_average_multi != -1 ] <- 0
      
      euclideanDist <- hclust(dist(t(segment_average_multi)))
      ordereuclideanDist <- euclideanDist$order
      segment_average_multi <- segment_average_multi[,ordereuclideanDist]
      output <- as.data.frame(cbind(segmentation_data[,c(1:3)], segment_average_multi))
    }
  }
  else if (all(c("chromosome", "start", "end", "segmean", "sample", "noise") %in% names(segmentation_data)))
    {
    segmentation_data$n.bp <- segmentation_data$end - segmentation_data$start
    single_data_filtered <- segmentation_data[,-which(names(segmentation_data) %in% c("start", "end", "noise"))]
    
    single_data_weights <- single_data_filtered %>% group_by(sample) %>% group_by(chromosome, .add = TRUE) %>% mutate(weight = n.bp / sum(n.bp),) 
    
    chromosomal_average_single <- single_data_weights %>% group_by(sample, chromosome) %>% summarise(segmean_weighted_mean = weighted.mean(segmean, weight, na.rm = TRUE))
    require(tidyverse)
    chromosomal_wide <- chromosomal_average_single %>% pivot_wider(names_from = sample, values_from = segmean_weighted_mean)
    chromosomal_average_single <- chromosomal_wide[,-1]
    
    chromosomal_average_single[chromosomal_average_single > 0 & chromosomal_average_single <= thresh ] <- 1
    chromosomal_average_single[chromosomal_average_single <= -thresh ] <- -1
    chromosomal_average_single[chromosomal_average_single != 1 & chromosomal_average_single != -1 ] <- 0
    
    euclideanDist <- hclust(dist(t(chromosomal_average_single)))
    ordereuclideanDist <- euclideanDist$order
    chromosomal_average_single <- chromosomal_average_single[,ordereuclideanDist]
    output <- cbind(chromosome = c(1:length(t(chromosomal_average_single)[1,])), chromosomal_average_single)
  } else {
    stop("Please verify that your input is a dataframe from the cumul.CNV function or is ordered accordingly.")
  }
  return(output)
}