###############
# GENERAL FUNCTIONS
###############

# .remove_ext: remove file path and extension from file name
.remove_ext <- function(filename) {
  pattern1 <- "\\.(fastq|fq|fasta|fa|sam|bam|bed|bedgraph|bg|bigwig|bw)"
  pattern2 <- "\\.(txt|tar|bz2|bz|csv).*"
  reduced <- tools::file_path_sans_ext(basename(filename))
  reduced <- sub(pattern1, "", reduced)
  (reduced <- sub(pattern2, "", reduced))
}

# .get_ext: get extension from a file name
.get_ext <- function(filename) {
  if (!is.character(filename)) {
    stop("the file name should be of type character()")
  }
  strip.file.frags <- function(x) {
    file.segs <- strsplit(x, ".", fixed = T)[[1]]
    lss <- length(file.segs)
    if (lss > 1) {
      out <- paste(file.segs[lss])
    } else {
      out <- ""
    }
    return(out)
  }
  return(sapply(filename, strip.file.frags))
}

# .remove_trail: remove trailing white spaces (or other character)
.remove_trail <- function(str, before = TRUE, after = TRUE, char = " ") {
  if (!is.character(str)) {
    warning("not a character() type")
    return(str)
  }
  ch <- substr(paste(char)[1], 1, 1)
  kk <- (length(str))
  if (kk < 1) return(str)
  for (cc in 1:kk) {
    if (isTRUE(before)) {
      while (substr(str[cc], 1, 1) == ch) {
        if (nchar(str[cc]) > 1) {
          str[cc] <- substr(str[cc], 2, nchar(str[cc]))
        } else {
          str[cc] <- gsub(ch, "", str[cc])
        }
      }
    }
    if (after) {
      while (substr(str[cc], nchar(str[cc]), nchar(str[cc])) == ch) {
        if (nchar(str[cc]) > 1) {
          str[cc] <- substr(str[cc], 1, nchar(str[cc]) - 1)
        } else {
          str[cc] <- gsub(ch, "", str[cc])
        }
      }
    }
  }
  return(str)
}

# .cmd: function to run commands either directly on bash or via WSL
.cmd <- function(cmd, intern = FALSE) {
  if (.Platform$OS.type != "windows") {
    system(command = cmd, intern = intern)
  } else {
    shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
  }
}

# .check_cmd: function to check whether a command exists and is callable
.check_cmd <- function(x) {
  test <- sprintf(
    "type -P %s &>//dev//null && echo 'Found' || echo 'Not Found'",
    x)
  .cmd(test, intern = TRUE)
}

# .create_dir: create a directory
.create_dir <- function(dir.name) {
  dir.create(dir.name, showWarnings = FALSE, recursive = TRUE)
}



# .isin and .isnotin: filtering one list by another list
.isin <- function(theList, toMatch){
  matches <- unique(grep(paste(toMatch,collapse = "|"),
                         theList, value = TRUE))
  return(matches)
}

.isnotin <- function(theList, toMatch){
  return(setdiff(theList, .isin(theList,toMatch)))
}

# .extract_numbers: extract numbers from a list
.extract_numbers <- function(string) {
  x <- unlist(regmatches(string, gregexpr('\\(?[0-9,.]+', string)))
  x <- as.numeric(gsub('\\(', '', gsub(',', '', x)))
}

###############
# PLOTTING FUNCTIONS
###############

# .theme_Publication: theme for ggplot2 figures
.theme_Publication <- function() {
  if (Sys.info()['sysname'] == "Windows") {
    windowsFonts(Helvetica = windowsFont("TT Helvetica"))
  }
  (ggthemes::theme_foundation(base_size   = 14,
                              base_family = "Helvetica") +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(face   = "bold",
                                            size   = ggplot2::rel(1.2),
                                            hjust  = 0.5),
        text                 = ggplot2::element_text(size   = 12),
        panel.background     = ggplot2::element_rect(colour = NA),
        plot.background      = ggplot2::element_rect(colour = NA),
        panel.border         = ggplot2::element_rect(colour = NA),
        axis.title           = ggplot2::element_text(face   = "bold",
                                                     size   = ggplot2::rel(1)),
        axis.title.y         = ggplot2::element_text(angle  = 90,
                                                     vjust  = 2),
        axis.title.x         = ggplot2::element_text(vjust  = -0.2),
        axis.text            = ggplot2::element_text(),
        axis.line            = ggplot2::element_line(colour = "black"),
        axis.ticks           = ggplot2::element_line(),
        panel.grid.major     = ggplot2::element_line(colour = "#f0f0f0"),
        panel.grid.minor     = ggplot2::element_blank(),
        legend.key           = ggplot2::element_rect(colour = NA),
        legend.position      = "bottom",
        legend.direction     = "horizontal",
        legend.key.size      = ggplot2::unit(0.7, "lines"),
        legend.spacing       = ggplot2::unit(0, "cm"),
        legend.title         = ggplot2::element_text(face   = "italic"),
        plot.margin          = ggplot2::unit(c(10, 5, 5, 5), "mm"),
        strip.background     = ggplot2::element_rect(colour = "#f0f0f0",
                                                     fill   = "#f0f0f0"),
        strip.text           = ggplot2::element_text(face   = "bold"))
    )
}

# .grid_arrange_shared_legend: share a single legend between multiple ggplots
.grid_arrange_shared_legend <- function(...,
                                        nrow = 1,
                                        ncol = length(list(...)),
                                        position = c("bottom", "right")) {
  requireNamespace("ggplot2")
  requireNamespace("gridExtra")
  requireNamespace("grid")

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = gridExtra::arrangeGrob(
                       do.call(gridExtra::arrangeGrob, gl),
                       legend, ncol = 1,
                       heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = gridExtra::arrangeGrob(
                       do.call(gridExtra::arrangeGrob, gl),
                       legend, ncol = 2,
                       widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
}

###############
# QC MODULE FUNCTIONS
###############

# .get_status: identify whether FASTQC modules were PASS, FAIL or WARN
.get_status <- function(fqcRes, module) {
  fqcRes$Summary[grep(module, fqcRes$Summary$module, ignore.case = TRUE), 1]
}

# .valid_modules: determine whether the modules chosen are valid
.valid_modules <- function(modules = "all"){
  allowed.modules <- c("Summary",
                       "Basic_Statistics",
                       "Per_base_sequence_quality",
                       "Per_sequence_quality_scores",
                       "Per_base_sequence_content",
                       "Per_sequence_GC_content",
                       "Per_base_N_content",
                       "Sequence_Duplication_Levels",
                       "Overrepresented_sequences",
                       "Adapter_Content",
                       "Kmer_content")

  # Modules
  if ( "all" %in% modules)
    modules <- allowed.modules
  else{
    # partial matching of module names
    modules <- .remove_trail(modules)
    modules <- gsub(" ", "_", modules)
    modules <- grep(pattern = paste(modules, collapse = "|"),
                    allowed.modules,
                    ignore.case = TRUE, value = TRUE)
    modules <- unique(modules)
    if (length(modules) == 0) {
      stop("Incorect module names provided. Allowed values include: \n\n",
           paste(allowed.modules, collapse = "\n- "))
    }
  }
  modules
}

###############
# DIFF ANALYSIS MODULE FUNCTIONS
###############

# .create_folders: create directory tree for diff analysis results
.create_folders <- function(modules = "LED") {
  mainfolder <- "DE_results"
  mods <- unlist(strsplit(modules, ""))
  subfolders <- file.path(mainfolder, "Common_results")
  if ("L" %in% mods) subfolders <-
    c(subfolders, file.path(mainfolder, "limma_results"))
  if ("E" %in% mods) subfolders <-
    c(subfolders, file.path(mainfolder, "edger_results"))
  if ("D" %in% mods) subfolders <-
    c(subfolders, file.path(mainfolder, "deseq2_results"))
  lapply(subfolders, function(x) {
    dir.create(x, showWarnings = FALSE, recursive = TRUE)
  })
}

# .extract_common: identify genes determined as differentially expressed by all
# the chosen programs
.extract_common <- function(results = results, p.value = 0.05) {
  DE <- results$diff_genes

  for (i in seq_along(DE)) {
    if (names(DE)[i] != "limma") {
      DE[[i]] <- lapply(DE[[i]], function(x) x <- x[x$FDR <= p.value, ])
    }
    if (names(DE)[i] == "limma") {
      DE[[i]] <- lapply(DE[[i]], function(x) x <- x[x$`adj.P.Val` <= p.value, ])
    }
  }

  cDE <- do.call(Map, c(list, DE))

  if (length(names(DE)) == 2) {
    if (all(names(DE) == c("limma", "edger"))) {
      cDE <- lapply(cDE, function(x) {
        Reduce(function(...) merge(..., by = "genes",
                                   all = FALSE), x)
      })
      sn <- (length(cDE[[1]]) - 12) / 2
      for (i in seq_along(cDE)) {
        cDE[[i]] <- cDE[[i]][c(1:7, (8 + sn):length(cDE[[i]]))]
        names(cDE[[i]]) <- c("genes",
                             paste0(names(cDE[[i]])[2:7], ".limma"),
                             paste0(names(cDE[[i]])[8:12], ".edger"),
                             gsub(".y", "",
                                  names(cDE[[i]])[13:length(names(cDE[[i]]))]))
      }
    }
    if (all(names(DE) == c("limma", "deseq2"))) {
      cDE <- lapply(cDE, function(x) {
        Reduce(function(...) merge(..., by = "genes",
                                   all = FALSE), x)
      })
      sn <- (length(cDE[[1]]) - 13) / 2
      for (i in seq_along(cDE)) {
        cDE[[i]] <- cDE[[i]][c(1:7, (8 + sn):length(cDE[[i]]))]
        names(cDE[[i]]) <- c("genes",
                             paste0(names(cDE[[i]])[2:7], ".limma"),
                             paste0(names(cDE[[i]])[8:13], ".deseq2"),
                             gsub(".y", "",
                                  names(cDE[[i]])[14:length(names(cDE[[i]]))]))
      }
    }
    if (all(names(DE) == c("edger", "deseq2"))) {
      cDE <- lapply(cDE, function(x) {
        Reduce(function(...) merge(..., by = "genes",
                                   all = FALSE), x)
      })
      sn <- (length(cDE[[1]]) - 12) / 2
      for (i in seq_along(cDE)) {
        cDE[[i]] <- cDE[[i]][c(1:7, (8 + sn):length(cDE[[i]]))]
        names(cDE[[i]]) <- c("genes",
                             paste0(names(cDE[[i]])[2:7], ".edger"),
                             paste0(names(cDE[[i]])[8:13], ".deseq2"),
                             gsub(".y", "",
                                  names(cDE[[i]])[14:length(names(cDE[[i]]))]))
      }
    }
  }
  if (length(names(DE)) == 3) {
    cDE <- lapply(cDE, function(x) {
      Reduce(function(...) merge(..., by = "genes",
                                 all = FALSE), x)
    })
    sn <- (length(cDE[[1]]) - 18) / 3
    for (i in seq_along(cDE)) {
      cDE[[i]] <- cDE[[i]][c(1:7, (8 + sn):(12 + sn),
                             (13 + (2 * sn)):length(cDE[[i]]))]
      names(cDE[[i]]) <- c("genes",
                           paste0(names(cDE[[i]])[2:7], ".limma"),
                           paste0(names(cDE[[i]])[8:12], ".edger"),
                           paste0(names(cDE[[i]])[13:18], ".deseq2"),
                           names(cDE[[i]])[19:length(names(cDE[[i]]))])
    }
  }
  cDE <- lapply(cDE, function(x) {
    x <- x[, which(grepl("genes|logFC|adj.P.Val|FDR", names(x)))]
    names(x) <- gsub("\\.x\\.", "\\.", names(x))
    names(x) <- gsub("\\.y\\.", "\\.", names(x))
    x$logFC.mean <- rowMeans(x[, which(grepl("logFC", names(x)))])
    return(x)
  })
  return(cDE)
}

# .venn2: create a venn diagram from two overlapping lists
.venn2 <- function(set1, set2, names, title) {
  stopifnot(length(names) == 2)

  #Form universe as union of sets
  universe <- sort(unique(c(set1, set2)))
  Counts <- matrix(0, nrow = length(universe), ncol = 2)
  colnames(Counts) <- names
  for (i in 1:length(universe)) {
    Counts[i, 1] <- universe[i] %in% set1
    Counts[i, 2] <- universe[i] %in% set2
  }
  par(mar = c(1, 1, 2, 1), oma = c(1, 1, 2, 1))
  limma::vennDiagram(vennCounts(Counts))
  mtext(title,outer = T, line = 1)
  tab <- data.frame(universe, Counts, stringsAsFactors = FALSE)
  colnames(tab) <- c("genes", names)
  return(tab)
}

# .venn3: create a venn diagram from three overlapping lists
.venn3 <- function(set1, set2, set3, names, title) {
  stopifnot( length(names) == 3)

  # Form universe as union of all three sets
  universe <- sort(unique(c(set1, set2, set3)))
  Counts <- matrix(0, nrow = length(universe), ncol = 3)
  colnames(Counts) <- names

  for (i in 1:length(universe)) {
    Counts[i, 1] <- universe[i] %in% set1
    Counts[i, 2] <- universe[i] %in% set2
    Counts[i, 3] <- universe[i] %in% set3
  }
  par(mar = c(1, 1, 2, 1), oma = c(1, 1, 2, 1))
  limma::vennDiagram(vennCounts(Counts))
  mtext(title, outer = T, line = 1)
  tab <- data.frame(universe, Counts, stringsAsFactors = FALSE)
  colnames(tab) <- c("genes", names)
  return(tab)
}

