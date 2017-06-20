#' Read a single FastQC data file
#' @description Read FastQC data into R
#' @param fqc Path to the FastQC result zip file to be imported.
#'
#' @return Returns an object of class fqRes containing the data from each of
#'  the FastQC modules.
#'
#' @examples
#' \dontrun{
#' read_fastqc("qcfile.zip")
#' }
#'
#' @export

read_fastqc <- function(fqc) {

  zf <- tools::file_path_sans_ext(fqc)
  utils::unzip(fqc, exdir = paste0(zf, "/.."), overwrite = TRUE)
  on.exit(unlink(zf, recursive = TRUE, force = TRUE))

  fd <- readLines(file.path(zf, "fastqc_data.txt"))
  sd <- readLines(file.path(zf, "summary.txt"))

  fd <- gsub("#", "", fd)
  fd <- fd[!grepl(">>END_MODULE", fd)]
  fd <- split(fd, cumsum(grepl("^>>", fd)))[-1]
  names(fd) <- sapply(fd, function(x) {
    gsub("^>>", "", gsub("\t.*", "", gsub(" ", "_", x[1])))
  })
  names(fd) <- tolower(names(fd))
  fd <- lapply(fd, function(x) {
    if (length(x) == 1) return(data.frame())
    x <- strsplit(x[-1], split = "\t")
    tab <- as.data.frame(do.call(rbind, x[-1]), stringsAsFactors = FALSE)
    for (i in 1:ncol(tab)) {
      if (!any(is.na(suppressWarnings(as.numeric(tab[,i]))))) {
        tab[,i] <- as.numeric(tab[,i])
      }
    }
    colnames(tab) <- x[[1]]
    tab
  })
  pct_dup <- 100 - as.numeric(colnames(fd$sequence_duplication_levels)[2])
  pct_dup <- round(pct_dup, digits = 2)
  fd$pct_duplication <- pct_dup

  colnames(fd$sequence_duplication_levels) <-
    unname(as.vector(fd$sequence_duplication_levels[1, ]))
  fd$sequence_duplication_levels <- fd$sequence_duplication_levels[-1, ]
  fd$sequence_duplication_levels[, 2] <-
    as.numeric(fd$sequence_duplication_levels[, 2])
  fd$sequence_duplication_levels[, 3] <-
    as.numeric(fd$sequence_duplication_levels[, 3])

  sd <- gsub("#", "", sd)
  sd <- strsplit(sd, split = "\t")
  sd <- as.data.frame(do.call(rbind, sd), stringsAsFactors = FALSE)
  colnames(sd) <- c("status", "module", "sample")

  ad <- vector("list", 14)
  for (i in 1:13) ad[[i]] <- fd[[i]]
  ad[[14]] <- sd
  names(ad) <- c(names(fd), "Summary")

  ad <- structure(ad, class = c(class(ad), "fqcRes"))
}
