#' binning target data and publicly available control data and mapping them to
#' the genome and generating annotation file for bins
#'
#' @param target_rgset RGchannelSet of data to be binned
#' @param ArrayType character representing a that Array Type used for acquiring
#' data
#' @return A list containing three objects, the mapped MethylSet of the target
#' data, the mapped MethylSet of the control Data and the annotation file of
#' the bins.
sampleBinContr <- function(target_rgset, ArrayType, controls) {
    #generate bins with some good default values
    anno_targets <- conumee::CNV.create_anno(array_type = ArrayType)
    # Illumina normalisation
    target_mset <- minfi::preprocessIllumina(target_rgset)
    target_mset_mapped <- minfi::mapToGenome(target_mset)
    target_mset_loaded <- conumee::CNV.load(target_mset) 
    
    if (is.null(controls)) {
      #load controls based on ArrayType
      if (ArrayType == "overlap" || ArrayType == "450k") {
        control_mset <- minfiData::MsetEx
      }
      if (ArrayType == "EPIC") {
        control_mset <- minfiDataEPIC::MsetEPIC
      }
    } else {
      control_mset <- controls
    }
    
    control_mset_loaded <- conumee::CNV.load(control_mset)
    
    # find overlapping probes between arraydata and annotations
    anno_targets@probes <-
        IRanges::subsetByOverlaps(anno_targets@probes, granges(target_mset_mapped))
    
    output <-
        list(
            "target_mset_loaded" = target_mset_loaded,
            "control_mset_loaded" = control_mset_loaded,
            "anno_targets" = anno_targets
        )
    return(output)
}