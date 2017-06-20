#' Design and contrast matrices
#'
#' \code{create_design} - Function to generate design matrix for edgeR
#' and limma+voom.
#'
#' @param se A SummarizedExperiment object containing the sample metadata
#'   information.
#' @param block.column The column in `samples` containing information
#'                     regarding block or batch effect.
#'
#' @return
#' \code{create_design} returns a design matrix.
#'
#' @rdname create_contrast
#'
#' @importFrom limma makeContrasts
#'
#' @export
#'
#' @examples
#' \dontrun{
#' design <- create_design(se, block.column = "replicate")
#' contrast <- create_contrast(design = design)
#' }

create_design <- function(se = NULL,
                          block.column = NULL) {

  dform <- formula(paste0("~0+", "Contrast",
                           ifelse(!is.null(block.column),
                                  paste0("+", block.column), "")))
  design <- model.matrix(dform, data = SummarizedExperiment::colData(se))
  print(as.data.frame.array(design))
  return(design)
}

#' Design and contrast matrices
#'
#' \code{create_contrast} - Function to generate contrast matrix for edgeR,
#' limma+voom, and DESEq2 analyses.
#'
#' @param design Design matrix produced by \code{create_design}.
#'
#' @return
#' \code{create_contrast} returns a contrast matrix.
#'
#' @rdname create_contrast
#'
#' @export

create_contrast <- function(design = NULL) {
  contrastColumns <- factor(colnames(design)[grep("Contrast",
                                                  colnames(design))],
                            levels = colnames(design)[grep("Contrast",
                                                           colnames(design))])
  x <- combn(contrastColumns, 2, simplify = F)
  y <- as.list(unlist(lapply(x, function(x) paste0(x[2], "-", x[1]))))
  contr.matrix <- limma::makeContrasts(contrasts = y, levels = colnames(design))
  print(contr.matrix)
  return(contr.matrix)
}
