get.ArrayType <- function(dataFiles) {
    ArrayType <- NULL
    if ((all(c("EPIC", "450K") %in% types)) &
        (length(types) == 2)) {
        ArrayType <- "combined"
    } else if (("450K" %in% types) & (length(types) == 1)) {
        ArrayType <- "450K"
    } else if (("EPIC" %in% types) & (length(types) == 1)) {
        ArrayType <- "EPIC"
    }
    return(ArrayType)
}

read.RGSet <- function(dataFiles, ArrayType) {
    types = unique(dataFiles$ArrayType)
    
    #read data transform to RGChannelSet
    if (ArrayType == "combined") {
        # separate file names
        data_EPIC <-
            dataFiles[which(dataFiles$ArrayType == "EPIC"), ]
        data_450k <-
            dataFiles[which(dataFiles$ArrayType == "450k"), ]
        # separately prepare 450k and 850k and then combine them
        rgset_EPIC <-
            read.metharray.exp(targets = data_EPIC, force = TRUE)
        rgset_450k <-
            read.metharray.exp(targets = data_450k, force = TRUE)
        # combine into array with 450K cpg sites
        target_rgset <-
            combineArrays(rgset_EPIC, rgset_450k, outType = "IlluminaHumanMethylation450k")
    }
    if (ArrayType == "450K") {
        rgset_450k <- read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_450k
    }
    if (ArrayType == "EPIC") {
        rgset_EPIC <- read.metharray.exp(targets = dataFiles, force = TRUE)
        target_rgset <- rgset_EPIC
        
    }
    return_list = list("ArrayType" = ArrayType, "RGSet" = target_rgset)
    return(return_list)
}

segment.Plot <-
    function(target_rgset,
             array_type,
             thresh,
             colour.amplification,
             colour.loss,
             detail.regions) {
        mSetsAnno <-  sampleBinContr(target_rgset, array_type)
        
        if (segmentationMode == "single") {
            simple.Seg(mSetsAnno,
                            array_type,
                            thresh,
                            colour.amplification,
                            colour.loss)
        } else if (segmentationMode == "multi") {
            multiSampleSeg(
                mSetsAnno,
                array_type,
                thresh,
                colour.amplification,
                colour.loss,
                detail.regions
            )
        } else if (segmentationMode == "all") {
            multiSampleSeg(
                mSetsAnno,
                array_type,
                thresh,
                colour.amplification,
                colour.loss,
                detail.regions
            )
            simple.Seg(mSetsAnno,
                            ArrayType,
                            thresh,
                            colour.amplification,
                            colour.loss)
            
        }
    }

cum.CNV <-
    function(dataFiles,
             segmentationMode = "all",
             thresh = 0.2,
             gamma = 5,
             colour.amplification = "red3",
             colour.loss = "blue4",
             detail.regions = NULL) {
        # check user input
        stopifnot("dataFiles must be a dataframe" = typeof(dataFiles) == "list")
        stopifnot("Basename must be a column of input dataframe dataFiles" =
                      ("Basename" %in% colnames(dataFiles)))
        stopifnot(
            "ArrayType must be a column of input dataframe dataFiles" = ("ArrayType" %in% colnames(dataFiles))
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
        
        # determine type of input files
        array_type <- get.ArrayType(dataFiles)
        # read in RGSet
        target_rgset <- read.RGSet(dataFiles, array_type)
        # segment and plot
        segment.Plot(
            target_rgset,
            array_type,
            thresh,
            colour.amplification,
            colour.loss,
            detail.regions
        )
        
    }