#' binning target data and publicly available control data and mapping them to
#' the genome and generating annotation file for bins using conumee2.0
#'
#' @param target_rgset RGchannelSet of data to be binned
#' @param ArrayType character representing a that Array Type used for acquiring
#' data
#' @return A list containing three objects, the mapped MethylSet of the target
#' data, the mapped MethylSet of the control Data and the annotation file of
#' the bins.
sampleBinContr2 <- function(target_rgset, ArrayType, controls) {
  #generate bins with some good default values
  
  
  if (ArrayType == "mouse") {
    anno_targets <- conumee2.0::CNV.create_anno(array_type = ArrayType, chrXY = TRUE)
    target_mset <- conumee2.0::CNV.import(target_rgset$array_type, target_rgset$directory, target_rgset$sample_sheet)
  
  }
  else if (ArrayType == "EPICv2") {
    anno_targets <- conumee2.0::CNV.create_anno(array_type = ArrayType)
    target_mset <- conumee2.0::CNV.import(target_rgset$array_type, target_rgset$directory, target_rgset$sample_sheet)
    
  }
    else {
      # Illumina normalisation
      anno_targets <- conumee2.0::CNV.create_anno(array_type = ArrayType)
      target_mset <- minfi::preprocessIllumina(target_rgset)
    }
  target_mset_loaded <- conumee2.0::CNV.load(target_mset) 
  
  
  if (is.null(controls))  {
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
  
  control_mset_loaded <- conumee2.0::CNV.load(control_mset)
  
  output <-
    list(
      "target_mset_loaded" = target_mset_loaded,
      "control_mset_loaded" = control_mset_loaded,
      "anno_targets" = anno_targets
    )
  return(output)
}