#' binning target data and publicly available control data and mapping them to 
#' the genome and generating annotation file for bins   
#'
#' @param target_rgset RGchannelset of data to be binned
#' @param ArrayType character representing a that Array Type used for aquiring
#' data
#' @return A list containing three entires, the mapped MethylSet of the target 
#' data, the mapped MethylSet of the control Data and the annotation file of 
#' the bins.

sampleBinContr <- function(target_rgset, ArrayType) {
  #generate bins
  anno_targets <- conumee::CNV.create_anno(bin_minprobes = 15, bin_minsize = 50000, bin_maxsize = 5000000, array_type = ArrayType)
  # Illumina normalisation
  target_mset <- preprocessIllumina(target_rgset)
  target_mset_mapped <- mapToGenome(target_mset)
  target_mset_loaded <- conumee::CNV.load(target_mset)
  
  #load controls based on ArrayType
  if (ArrayType == "overlap" || ArrayType == "450k"){
    control_mset <- minifData::MsetEx
  }
  if (ArrayType == "EPIC"){
    control_mset <- minfiDataEPIC::MsetEPIC
  }
  control_mset_loaded <- conumee::CNV.load(control_mset)
  
  # find overlapping probes between arraydata and annotations
  anno_targets@probes <- subsetByOverlaps(anno_targets@probes, granges(target_mset_mapped))
  
  output <- list("target_mset_loaded" = target_mset_loaded, "control_mset_loaded" = control_mset_loaded, "anno_targets" = anno_targets)
  return(output)
  
}