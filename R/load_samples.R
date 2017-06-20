#' Load sample metadata file
#'
#' Loads a tab-delimited file containing the sample metadata.  There should be
#' at least two columns - one called basename (containing the basename of the
#' samples) and another containing the contrasts of interest.  Another column
#' containing an additional factor of interest which should be taken into
#' account may also be included.
#'
#' @param sample.path Character string.  Full path to sample metadata file.
#' @param contrast.column Chacter string.  Column in sample metadata file
#'   that contains the contrasts of interest.
#' @param block.column Character string.  Column in sample metadata file
#'   that contains blocking factors.
#' @param contrast.levels Character vector.  The contrast levels listed
#'   in the order they should be compared as present in the contrast
#'   column of the metadata file (e.g. c("WT", "KO")).
#'
#' @return A dataframe of the sample metadata.
#'
#' @export

load_samples <- function(sample.path = NULL,
                         contrast.column = NULL,
                         block.column = NULL,
                         contrast.levels) {

    # Read in the sample metadata file.
    sampleData <- read.table(sample.path,
                           header = T,
                           sep = "\t",
                           strip.white = TRUE)

    # Check the sample metadata file for errors.
    if (!("basename" %in% names(sampleData))) {
        stop(cat("Error: the column 'basename' is missing
             from the metadata file\n\n"))
    }
    if (!(contrast.column %in% names(sampleData))) {
        stop(cat("Error: the specified 'contrast' column is missing
             from the metadata file"))
    }
    if (!is.null(block.column) && !(block.column %in% names(sampleData))) {
        stop(cat("Error: the specified 'block' column is missing
             from the metadata file"))
    }
    sampleData[, contrast.column] <- as.factor(sampleData[, contrast.column])
    if (!all(contrast.levels %in%
             as.character(sampleData[, contrast.column]))) {
        stop(cat("Error: one or more of the levels are missing
             from the 'contrast' column"))
    }
    if (min(table(sampleData[, contrast.column])) < 2) {
        stop(cat("Error: one (or more) of the contrasts has no replicates;
             analysis can not proceed"))
    }

    # Relevel the contrasts.
    sampleData[, contrast.column] <- factor(sampleData[, contrast.column],
                                          levels = contrast.levels)

    # Make sure the blocking factor column consists of factors.
    if (!is.null(block.column) &&
        is.numeric(sampleData[, block.column])) {
        sampleData[, block.column] <-
        factor(sampleData[, block.column])
    }

    # Order the metadata by contrast levels.
    sampleData <- sampleData[order(sampleData[, contrast.column]),]

    # Show the sample metadata on the screen.
    cat("Sample metadata:\n")
    print(sampleData)

    # Rename the column containing the contrasts to "Contrast".
    names(sampleData)[names(sampleData) == contrast.column] <- "Contrast"

    # Convert the rownames to the basename
    rownames(sampleData) <- sampleData$basename
    return(sampleData)
}
