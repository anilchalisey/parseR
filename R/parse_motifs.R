#' Parse HOMER motif files
#'
#' @description \code{homer_known} is a function to parse the 'known' motifs
#'   identified by HOMER.
#'
#' @param motif_dir Directory in which the motif files created by HOMER
#'                  are located
#'
#' @return
#' An R object with the known motifs in an easily readable format.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' motifdir <- system.file("extdata", "motifs.tar.gz", package = "dealr")
#' extract <- sprintf("tar -xvzf %s", motifdir)
#' system(extract)
#' motifs <- list.files("motifs")
#' hk <- homer_known("motifs")
#' hd <- homer_denovo("motifs")
#' compress <- sprintf("tar -cvzf %s %s", motifs, motifdir)
#' system(compress)
#' unlink("motifs", recursive = TRUE)
#' }

homer_known <- function(motif_dir = "motifs") {

  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("stringr")

  # Read the known results file produced by HOMER
  known_motifs <- read.table(file.path(motif_dir, "knownResults.txt"),
                             sep = "\t", header = T, comment.char = "")

  # Parse the table to a suitable format
  colnames(known_motifs) <- c("Name", "Consensus", "Pval",
                              "Log-pval", "Qval", "Tcount",
                              "T", "Bcount", "B")

  known_motifs <- known_motifs %>%
    dplyr::mutate(T = gsub(pattern = "%", replacement = "", x = T) %>%
                    as.numeric %>%
                    magrittr::divide_by(100),
                  B = gsub(pattern = "%", replacement = "", x = B) %>%
                    as.numeric %>%
                    magrittr::divide_by(100)) %>%
    dplyr::select(-Tcount, -Bcount)

  known_motifs <- known_motifs %>%
    dplyr::mutate(
      Name = gsub(pattern = ")/.*", replacement = "", x = Name)) %>%
    tidyr::separate(
      ., Name, into = c("Name", "Family"), sep = "\\(", extra = "drop")

  return(known_motifs)
}

#' Parse HOMER motif files
#'
#' \code{homer_denovo} is a function to parse the 'denovo' motifs identified by
#'   HOMER.
#'
#' @return
#' An R object with the \emph{de novo} motifs in an easily readable format.
#'
#' @rdname homer_known
#'
#' @export

homer_denovo <- function(motif_dir = "motifs") {

  ### Folder motifs
  # Files
  motif_fnames <- list.files(file.path(motif_dir, "homerResults"),
                             full.names = TRUE)
  motif_fnames <- motif_fnames[grepl(pattern = ".motif", x = motif_fnames,
                                     fixed = TRUE)]
  motif_fnames <- motif_fnames[!grepl(pattern = "similar", x = motif_fnames,
                                      fixed = T)]
  motif_fnames <- motif_fnames[!grepl(pattern = "RV", x = motif_fnames,
                                      fixed = T)]

  # Separate PWMS
  homer_PWMs <- lapply(motif_fnames, read.table, skip = 1,
                       col.names = c("A", "C", "G", "T"))

  # Reformat info
  folder_motifs <- lapply(motif_fnames, read.table, nrows = 1, sep = "\t") %>%
    Reduce(rbind, .) %>%
    as.data.frame()
  colnames(folder_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval",
                               "Placeholder", "Occurence")

  # Split names
  folder_motifs <- tidyr::separate(folder_motifs, Name,
                                   into = c("Name", "Guess"),
                                   sep = ",", extra = "drop")

  # Rename PWMs
  names(homer_PWMs) <- folder_motifs$Name

  ### File motifs
  # Read from files starting with >
  heap <- readLines(file.path(motif_dir, "homerMotifs.all.motifs"))
  heap <- heap[grep(pattern = ">", x = heap)]

  # Reformat info
  file_motifs <- heap %>%
    stringr::str_split(pattern = "\t") %>%
    as.data.frame() %>% t() %>% as.data.frame()

  colnames(file_motifs) <- c("Consensus", "Name", "Log-odds", "Log-pval",
                             "Placeholder", "Occurence", "Statistics")

  rownames(file_motifs) <- NULL

  file_motifs <- file_motifs %>%
    dplyr::mutate(`Log-odds` = as.numeric(as.character(`Log-odds`)),
                  `Log-pval` = as.numeric(as.character(`Log-pval`)),
                  Placeholder = as.integer(as.character(Placeholder)))

  ### Merge info
  # Merge frames
  homer_motifs <- suppressWarnings(
    dplyr::left_join(folder_motifs, file_motifs,
                     by = c("Consensus", "Name", "Log-odds",
                            "Log-pval", "Placeholder", "Occurence")))
  # Split columns
  homer_motifs <- tidyr::separate(homer_motifs, Occurence,
                                  into = c("T", "B", "P"),
                                  sep = ",", extra = "drop")
  homer_motifs <- tidyr::separate(homer_motifs, Statistics,
                                  into = c("Tpos", "Tstd", "Bpos", "Bstd",
                                           "StrandBias", "Multiplicity"),
                                  sep = ",", extra = "drop")

  # Reformat
  homer_motifs <- homer_motifs %>%
    dplyr::mutate(Consensus = gsub(pattern = ">", replacement = "",
                                   x = Consensus),
                  Guess = gsub(pattern = "BestGuess:", replacement = "",
                               x = Guess) %>%
                    gsub(pattern = "/.*", replacement = "", .),
                  T = stringr::str_split(string = T,
                                         pattern = ":|\\(|\\)|%") %>%
                    sapply(function(x) as.numeric(x[3]) / 100),
                  B = stringr::str_split(string = B,
                                         pattern = ":|\\(|\\)|%") %>%
                    sapply(function(x) as.numeric(x[3]) / 100),
                  P = gsub(pattern = "P:", replacement = "", x = P) %>%
                    as.numeric,
                  Tpos  =  gsub(pattern = "Tpos:",
                                replacement = "", x = Tpos) %>%
                    as.numeric,
                  Tstd = gsub(pattern = "Tstd:",
                              replacement = "", x = Tstd) %>%
                    as.numeric,
                  Bpos = gsub(pattern = "Bpos:",
                              replacement = "", x = Bpos) %>%
                    as.numeric,
                  Bstd = gsub(pattern = "Bstd:",
                              replacement = "", x = Bstd) %>%
                    as.numeric,
                  StrandBias = gsub(pattern = "StrandBias:",
                                    replacement = "", x = StrandBias) %>%
                    as.numeric,
                  Multiplicity = gsub(pattern = "Multiplicity:",
                                      replacement = "", x = Multiplicity) %>%
                    as.numeric) %>%
    dplyr::mutate(Orientation = ifelse(test = grepl(pattern = "RV",
                                                    x = motif_fnames),
                                       yes = "reverse", no = "forward"),
                  Length = sapply(homer_PWMs, nrow))

  ### Trim
  final_motifs <- !is.na(homer_motifs$Tpos)
  homer_motifs <- subset(homer_motifs, final_motifs, select = -Placeholder)
  homer_PWMs <- homer_PWMs[final_motifs]

  list(homer_motifs = homer_motifs, homer_PWMs = homer_PWMs)
}
