#'Determines the array type of the user input dataframe
#'
#'@param dataFiles Dataframe with a column labelled "ArrayType" containing the arraytype (450K/EPIC/EPIC2) for each sample.
#'@return A string with either "450K" (only 450K samples), "EPIC" (only EPIC samples), "combined" (both 450K and EPIC samples) and EPIC2 (only EPIC2)
get.ArrayType <- function(dataFiles) {
    ArrayType <- NULL
    types = unique(dataFiles$ArrayType)
    if ((all(c("EPIC", "450K") %in% types)) &
        (length(types) == 2)) {
        ArrayType <- "combined"
    } else if (("450K" %in% types) & (length(types) == 1)) {
        ArrayType <- "450K"
    } else if (("EPIC" %in% types) & (length(types) == 1)) {
        ArrayType <- "EPIC"
    } else if (("EPIC2" %in% types) & (length(types) == 1)) {
        ArrayType <- "EPIC2"
    }
    return(ArrayType)
}

#' Determines the conumee version from the array type
#'
#' @param ArrayType The ArrayType as returned by get.ArrayType
#'
#' @return The conumee version (either 1 or 2)
get.ConumeeVersion <- function(ArrayType) {
    v <- 1
    if (ArrayType %in% c("EPIC2")) {
        v <- 2
    }
    return(v)
}

#' Reads the specified methylation arrays into an RGSet. Note that any combination of 450K and EPIC arrays will be coerced into an RGSet of 450K type.
#' @param dataFiles Dataframe with a column batch name as requested by minfi for reading in experiments.
#' @param ArrayTye A string (either "450K", "EPIC", "combined" or "EPIC2")
#'
#' @return A list of the RGSet of the target data, the control data and the annotation data
read.RGSet <- function(dataFiles, ArrayType) {
    stopifnot(
        "Only 450K, EPIC and combined are permitted as ArrayType at the moment. EPIC2 is still missing" =
            (ArrayType %in% c("450K", "EPIC", "combined"))
    )
    types = unique(dataFiles$ArrayType)
    
    #read data transform to RGChannelSet
    if (ArrayType == "combined") {
        # separate file names
        data_EPIC <-
            dataFiles[which(dataFiles$ArrayType == "EPIC"),]
        data_450k <-
            dataFiles[which(dataFiles$ArrayType == "450k"),]
        # separately prepare 450k and 850k and then combine them
        rgset_EPIC <-
            minfi::read.metharray.exp(targets = data_EPIC, force = TRUE)
        rgset_450k <-
            minfi::read.metharray.exp(targets = data_450k, force = TRUE)
        # combine into array with 450K cpg sites
        target_rgset <-
            minfi::combineArrays(rgset_EPIC, rgset_450k, outType = "IlluminaHumanMethylation450k")
    }
    if (ArrayType == "450K") {
        rgset_450k <-
            minfi::read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_450k
    }
    if (ArrayType == "EPIC") {
        rgset_EPIC <-
            minfi::read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_EPIC
        
    }
    
    return(target_rgset)
}

#' Segments and visualizes DNA methylation data and generative cumulative plots
#' @param dataFiles A dataframe with columns ArrayType and Basename
#' @param segmentationMode specifying the segmentation mode. Allowed values are single/multi/all.
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param gamma A positive integer (>0) indicating the threshold for the multi-sample segmentation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param conumee.version The conumee version to use.
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
segment.Plot <-
    function(target_rgset,
             array_type,
             segmentationMode,
             gamma,
             thresh,
             colour.amplification,
             colour.loss,
             detail.regions,
             conumee.version) {
      if(conumee.version == 1) {
        require(conumee)
        mSetsAnno <-  sampleBinContr(target_rgset, array_type)
      } else {
        require(conumee2.0)
        mSetsAnno <-  sampleBinContr(target_rgset, array_type)
      }
        
        if (segmentationMode == "single") {
            if (conumee.version == 1) {
                singleSampleSeg(mSetsAnno,
                                thresh,
                                colour.amplification,
                                colour.loss)
            } else{
                singleSampleSeg2(mSetsAnno,
                                 thresh,
                                 colour.amplification,
                                 colour.loss)
            }
            
        } else if (segmentationMode == "multi") {
            multiSampleSeg(mSetsAnno,
                           thresh,
                           colour.amplification,
                           colour.loss,
                           detail.regions)
        } else if (segmentationMode == "all") {
            multiSampleSeg(mSetsAnno,
                           thresh,
                           colour.amplification,
                           colour.loss,
                           detail.regions)
            if (conumee.version == 1) {
                singleSampleSeg(mSetsAnno,
                                thresh,
                                colour.amplification,
                                colour.loss)
            } else{
                singleSampleSeg2(mSetsAnno,
                                 thresh,
                                 colour.amplification,
                                 colour.loss)
            }
            
        }
    }

#' compute the cumulative CNV plots
#' @param dataFiles A dataframe with columns ArrayType and Basename
#' @param segmentationMode specifying the segmentation mode. Allowed values are single/multi/all.
#' @param thresh A positive float (>=0) indicating the threshold for an abberation.
#' @param gamma A positive integer (>0) indicating the threshold for the multi-sample segmentation.
#' @param colour.amplification Colour for amplification
#' @param colour.loss Colour for loss
#' @param detail.regions Either NULL or a vector of gene names.
#' @param conumee.version The version of conumee to use (either 1 or 2). 1 is incompatible with mouse or EPICv2 arrays. NULL will set the version heuristically to 1 for 450K, EPIC and to 2 for Mouse and EPICv2
#'
#' @return Nothing. Will print the figures to the default plotting terminal.
#' @export
cum.CNV <-
    function(dataFiles,
             segmentationMode = "all",
             thresh = 0.2,
             gamma = 5,
             colour.amplification = "red3",
             colour.loss = "blue4",
             detail.regions = NULL,
             conumee.version = NULL) {
        # check user input
        stopifnot("dataFiles must be a dataframe" = typeof(dataFiles) == "list")
        stopifnot("Basename must be a column of input dataframe dataFiles" =
                      ("Basename" %in% colnames(dataFiles)))
        stopifnot(
            "ArrayType must be a column of input dataframe dataFiles" = ("ArrayType" %in% colnames(dataFiles))
        )
        stopifnot(
            "Parameter segmentationMode must be one of single/multi/all" = segmentationMode %in% c("multi", "all", "single")
        )
        stopifnot("Parameter thresh must be a float >=0" = (thresh >= 0) &&
                      (typeof(thresh) == "double"))
        stopifnot("Parameter gamma must be an integer >0" = (gamma > 0) &&
                      (typeof(gamma) == "double"))
        stopifnot(
            "Parameter colour.amplification must be a string" = typeof(colour.amplification) ==
                "character"
        )
        stopifnot("Parameter colour.loss must be a string" = typeof(colour.loss) ==
                      "character")
        stopifnot(
            "Parameter colour.loss must be a string" = (detail.regions == NULL) &
                (typeof(detail.regions) == "list")
        )
        stopifnot(
            "Parameter detail.regions must either be NULL or a vector of strings" =
                (
                    is.null(detail.regions) | typeof(detail.regions) == "character"
                )
        )
        stopifnot(
            "Parameter conumee.version must either be NULL, 1 or 2" = (conumee.version %in% c(1, 2) ||
                                                                           is.null(conumee.version))
        )
        
        # determine type of input files
        array_type <- get.ArrayType(dataFiles)
        # determine the conumee version (if not specified by the user)
        if (is.null(conumee.version)) {
            conumee.version <- get.ConumeeVersion(array_type)
        }
        # read in RGSet
        target_rgset <- read.RGSet(dataFiles, array_type)
        
        # segment and plot
        segment.Plot(
            target_rgset,
            array_type,
            segmentationMode,
            gamma,
            thresh,
            colour.amplification,
            colour.loss,
            detail.regions,
            conumee.version
        )
        
    }