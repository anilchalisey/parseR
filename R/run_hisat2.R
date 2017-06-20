#' Wrapper script to run HISAT2.
#'
#' @description Script to align reads to a reference genome using hisat2.  This
#'   requires an existing index which may be created using hisat2 itself.
#'   Commonly used genome indices may also be downloaded from the HISAT2
#'   homepage.
#'
#' @param hisat2 Path to hisat2 (if using WSL, then this should be the
#'   full path on the linux subsystem)
#' @param idx The basename of the index for the reference genome. The basename
#' is the name of any of the index files up to but not including the final
#' .1.ht2, etc.
#' @param mate1 Comma-separated list of files containing mate 1s (filename
#' usually includes _1)
#' @param mate2 Comma-separated list of files containing mate 2s (filename
#' usually includes _2). Sequences specified with this option must correspond
#' file-for-file and read-for-read with those specified in .
#' @param fastq Logical indicating if reads are FASTQ files.
#' @param fasta Logical indicating if reads are FASTA files.
#' @param softClipPenalty Sets the maximum (MX) and minimum (MN) penalties for
#' soft-clipping per base, both integers.  Must be given in the format "MX,MN".
#' @param noSoftClip Logical indicating whether to disallow soft-clipping.
#' @param noSplice Logical indicating whether to switch off spliced alignment,
#' e.g., for DNA-seq analysis.
#' @param knownSplice Path to text file containing known splice sites.
#' @param strand Specify strand-specific information.  Default is unstranded.
#' @param tmo Logical indicating whether to report only those reads aligning to
#' known transcripts.
#' @param maxAlign Integer indicating the maximum number of distinct primary
#' alignments to search for each read.
#' @param secondary Logical indicating whether to report secondary alignments.
#' @param minInsert The minimum fragment length for valid paired-end alignments.
#' This option is valid only with noSplice = TRUE.
#' @param maxInsert The maximum fragment length for valid paired-end alignments.
#' This option is valid only with noSplice = TRUE.
#' @param nomixed By default, when hisat2 cannot find a concordant or discordant
#' alignment for a pair, it then tries to find alignments for the individual
#' mates. If TRUE, this option disables that behavior.
#' @param nodiscordant By default, hisat2 looks for discordant alignments if it
#' cannot find any concordant alignments.  If true, this option disables that
#' behavior.
#' @param threads Integer value indicating number of parallel threads to use.
#' @param rgid Character string, to which the read group ID is set.
#' @param quiet If TRUE, print nothing except alignments and serious errors.
#' @param non_deterministic When set to TRUE, HISAT2 re-initializes its
#' pseudo-random generator for each read using the current time.
#'
#' @return
#' Alignment file in SAM format
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_hisat2(hisat2 = "hisat2", idx = "../prana/data-raw/index/UCSC.hg19",
#' mate1 = "../prana/data-raw/seqFiles/HB1_sample_1.fastq.gz",
#' mate2 = "../prana/data-raw/seqFiles/HB1_sample_2.fastq.gz",
#' fastq = TRUE, fasta = FALSE, softClipPenalty = NULL, noSoftClip = FALSE,
#' noSplice = FALSE, knownSplice = NULL, strand = NULL, tmo = FALSE,
#' maxAlign = NULL, secondary = FALSE, minInsert = NULL, maxInsert = NULL,
#' nomixed = FALSE, nodiscordant = FALSE,
#' threads = (parallel::detectCores() - 1), rgid = NULL, quiet = FALSE,
#' non_deterministic = TRUE)
#' }

run_hisat2 <- function(
  # Path to hisat2
  hisat2 = "hisat2",
  # Mandatory
  idx = NULL,
  mate1 = NULL,
  mate2 = NULL,
  # input options
  fastq = TRUE,
  fasta = FALSE,
  softClipPenalty = NULL,
  noSoftClip = FALSE,
  # spliced alignment options
  noSplice = FALSE, # If TRUE - DNA-seq; if FALSE - RNA-seq
  knownSplice = NULL,
  strand = NULL,
  tmo = FALSE,
  # reporting options
  maxAlign = NULL,
  secondary = FALSE,
  # paired end options
  minInsert = NULL,
  maxInsert = NULL,
  nomixed = FALSE,
  nodiscordant = FALSE,
  # other
  threads = 1,
  rgid = NULL,
  quiet = FALSE,
  non_deterministic = FALSE) {

  if ((isTRUE(fastq)) && (isTRUE(fasta))) {
    stop("fastq = TRUE and fasta = TRUE cannot both be true.")
  }
  if ((!isTRUE(fastq)) && (!isTRUE(fasta))) {
    stop("fastq = FALSE and fasta = FALSE cannot both be true.")
  }

  bnames <- .remove_ext(mate1)
  bnames <- sub("_1$", "", bnames)
  outfile <- paste0(bnames, ".sam")
  logfile <- paste0(bnames, ".log")

  if (is.null(mate2)) {
    paired <- FALSE
    maxInsert <- NULL
    minInsert <- NULL
    nomixed <- FALSE
    nodiscordant <- FALSE
    run.hisat <- sprintf("%s -p %s -x %s -1 %s",
                         hisat2, threads, idx, mate1)
  } else {
    paired <- TRUE
    run.hisat <- sprintf("%s -p %s -x %s -1 %s -2 %s",
                         hisat2, threads, idx, mate1, mate2)
  }

  if (is.null(hisat2) | is.null(mate1) | is.null(mate2) | is.null(idx)) {
    stop("The options hisat2, mate1, mate2, and idx are mandatory")
  }

  if (isTRUE(fastq)) {run.hisat <-
    sprintf("%s -q", run.hisat)}
  if (isTRUE(fasta)) {run.hisat <-
    sprintf("%s -f", run.hisat)}
  if (!is.null(softClipPenalty)) {run.hisat <-
    sprintf("%s --sp %s", run.hisat, softClipPenalty)}
  if (isTRUE(noSoftClip)) {run.hisat <-
    sprintf("%s --no-softclip", run.hisat)}
  if (isTRUE(noSplice)) {run.hisat <-
    sprintf("%s --no-spliced-alignment", run.hisat)}
  if (!is.null(knownSplice)) {run.hisat <-
    sprintf("%s --known-splicesite-infile", run.hisat, knownSplice)}
  if (!is.null(strand)) {run.hisat <-
    sprintf("%s --rna-strandedness", run.hisat, strand)}
  if (isTRUE(tmo)) {run.hisat <-
    sprintf("%s --tmo", run.hisat)}
  if (!is.null(maxAlign)) {run.hisat <-
    sprintf("%s -k %s", run.hisat, maxAlign)}
  if (isTRUE(secondary)) {run.hisat <-
    sprintf("%s --secondary", run.hisat)}
  if (!is.null(maxInsert)) {run.hisat <-
    sprintf("%s -X %s", run.hisat, maxInsert)}
  if (!is.null(minInsert)) {run.hisat <-
    sprintf("%s -I %s", run.hisat, minInsert)}
  if (isTRUE(nomixed)) {run.hisat <-
    sprintf("%s --no-mixed", run.hisat)}
  if (isTRUE(nodiscordant)) {run.hisat <-
    sprintf("%s --no-discordant", run.hisat)}
  if (!is.null(rgid)) {run.hisat <-
    sprintf("%s -rg-id %s", run.hisat, rgid)}
  if (isTRUE(quiet)) {run.hisat <-
    sprintf("%s --quiet", run.hisat)}
  if (isTRUE(non_deterministic)) {run.hisat <-
    sprintf("%s --non-deterministic", run.hisat)}
  run.hisat <- sprintf("%s -S %s", run.hisat, outfile, logfile)
  run.hisat <- paste(run.hisat, "2>>", logfile)

  # For some reason, hisat2 throws an error if the bash window is not open
  # in the background on Windows 10

  if (.Platform$OS.type == "windows") shell("start cmd /k bash.exe")

  for (i in seq_along(run.hisat)) {
    .cmd(run.hisat[i])
  }
  return(outfile)
}
