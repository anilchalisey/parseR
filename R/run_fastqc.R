#' Run FastQC
#' @description Run the FastQC tool
#' @param fq.files Vector of the full names and paths of the fastq files.
#' @param out.dir Path to the directory to which to write the results.  If NULL,
#'   which is the default, a directory named "FASTQC" is created in the current
#'   working directory.
#' @param threads The number of threads to use.  The default is 1.
#' @param fastqc The path to the fastqc executable.  On Unix systems, if the
#'   executable is in $PATH, then it may be left as the default ("fastqc").
#'   If it is not in $PATH, then the absolute path should be given.  If using
#'   the WSL on Windows 10, then the path must be the absolute path in WSL.
#'
#' @return A directory containing the QC reports
#'
#' @examples
#' \dontrun{
#' fq.files <- list.files(path = "data-raw", pattern = "*.fastq$")
#' fastqc(fq.files, out.dir = "FASTQ", threads = 2, fastqc = "~/bin/fastqc")
#' }
#'
#' @export

run_fastqc <- function(fq.files = NULL,
                       out.dir  = NULL,
                       threads  = 1,
                       fastqc   = "fastqc") {

  if (.check_cmd(fastqc) == "Not Found") {
    stop("the fastqc path is invalid or it is not installed correctly.")
  }

  fq <- paste(fq.files, collapse = " ")
  if (is.null(out.dir)) out.dir <- "FASTQC"
  .create_dir(out.dir)

  fastqc.run <- sprintf('%s %s --outdir=%s --threads=%s',
                        fastqc, fq, out.dir, threads)
  .cmd(fastqc.run)
}

