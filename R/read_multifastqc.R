#' Read a a directory of FastQC data files
#' @description Read a directory of FastQC data files into R
#' @param fqc.dir Path to the directory of FASTQC zip files to be imported.
#'
#' @return Returns an object of class fqRes_multi.  Each object within the list
#'  is a fqRes object containing the FastQC data for a specific sample.
#'
#' @examples
#' \dontrun{
#' read_multifastqc("qcfile.zip")
#' }
#'
#' @export

read_multifastqc <- function(fqc.dir) {
  fqc.files <- list.files(path = fqc.dir,
                          pattern = "*fastqc.zip$",
                          full.names = TRUE)

  if (length(fqc.files) == 0) {
    stop("Can't find any fastqc files in the directory specified.")
  }

  res <- lapply(fqc.files, read_fastqc)
  names(res) <- gsub("_fasqc.zip", "", basename(fqc.files))
  res <- structure(res, class = c(class(res), "fqcRes_multi"))
}
