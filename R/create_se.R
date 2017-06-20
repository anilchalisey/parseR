#' Create Summarized Experiment
#'
#' @param samples Sample metadata object created by \code{load_samples()}.
#' @param counts Count matrix where row names are features and column names are
#'   samples.  The column names must match the basenames in the samples object.
#' @param experimentTitle Title to be stored in the summarizedExperiment
#'   metadata.
#'
#' @return
#' A SummarizedExperiment Object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#'
#' @export


create_se <- function(samples = NULL,
                      counts = NULL,
                      experimentTitle = NULL) {

  if (ncol(counts) != nrow(samples)) {
    stop("The number of samples in the count matrix and the sample metadata
         do not match")
  }
  if (!(all(colnames(counts) %in% samples$basename))) {
    stop("The names of the samples in the count matrix do not match the sample
         basenames")
  }

  counts <- counts[, as.vector(samples$basename)]

  se <-
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = samples,
    metadata = S4Vectors::SimpleList(experimentTitle = experimentTitle)
    )

  rownames(se) <- rownames(counts)
  return(se)
}
