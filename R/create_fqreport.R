#' Create a QC report
#' @description Create an HTML file containing FastQC reports of one or multiple
#'  files. Inputs can be either a directory containing multiple FastQC zipped
#'  files or a single FastQC file.
#'
#' @param fqc Either a single zipped FASTQC file or a directory containing such
#'   files.
#' @param outdir The directory to which to write out the results.  The default
#'   is the current working directory.
#' @param preview Logical indicating whether to open the HTML file after it is
#'   created.
#' @param experiment Brief experiment description
#' @param author Author name
#'
#' @importFrom rmarkdown render
#' @import tools
#'
#' @return An HTML file (or files) containing FASTQC reports.
#' @export

create_fqreport <- function(fqc = NULL,
                            outdir = ".",
                            experiment = NULL,
                            author = NULL,
                            preview = TRUE) {

  if (dir.exists(fqc)) {
    qc.files <- list.files(path = fqc, pattern = "*fastqc.zip$",
                           full.names = TRUE)
    if (length(qc.files) == 0) {
      stop("Can't find any FASTQC files in this location")
    }
  }
  if (tools::file_ext(fqc) == "zip") {
    qc.files <- fqc
    if (length(qc.files) == 0) {
      stop("This file does not end in \"zip\".  Is this a zipped FASTQC file?")
    }
  }

  .create_dir(outdir)

  if (length(qc.files) == 1) {
    template <- system.file("extdata", "single_sample.Rmd", package = "parseR")
    qc.files <- normalizePath(qc.files)
    output <- gsub("_fastqc.zip", ".html", basename(qc.files))
    output <- file.path(normalizePath(outdir), output)
    rmarkdown::render(template,
                      output_file = output,
                      params = list(qc.files = qc.files,
                                    experiment = experiment,
                                    author = author))
  } else {
    template <- system.file("extdata", "multi_sample.Rmd", package = "parseR")
    qc.files <- normalizePath(qc.files)
    output <- file.path(normalizePath(outdir), "multiQCreport.html")
    rmarkdown::render(template,
                      output_file = output,
                      params = list(qc.files = qc.files,
                                    experiment = experiment,
                                    author = author))
    }
  if (isTRUE(preview)) utils::browseURL(output)
}

