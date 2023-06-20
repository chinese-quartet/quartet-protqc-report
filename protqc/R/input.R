#' Input files
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @import stats
#' @import utils
#' @importFrom data.table fread
#' @importFrom dplyr %>%
#' @importFrom dplyr rename_with
#' @export

input_data <- function(exp_path, meta_path) {

  # Load data ------------------------------------------------
  expr <- fread(exp_path)
  meta <- fread(meta_path)

  expr <- as.data.frame(expr)
  meta <- meta %>%
    rename_with(tolower)

  # Check data format ----------------------------------------
  # Compatible with version 0.1.x.
  # Caution: The data format is changed from version 0.2.0.
  # metadata file: name -> library, sample
  # data file: rowname -> Type, Feature, xxx

  if(length(which(duplicated(colnames(expr))))) {
    stop('Duplicated column names in data.')
    print(1)
  }

  col_check1 <- c("name", "library", "sample") %in% colnames(meta)
  if ((col_check1[1] | col_check1[2]) & col_check1[3]) {

    meta_final <- meta %>% select(any_of(c("name", "library", "sample")))
    colnames(meta_final) <- c("library", "sample")

    if(length(which(duplicated(meta_final$library)))) {

      stop('Duplicated column names in metadata.')
      print(1)
    }
  } else {

    stop('The columns named "library" and "sample" are required in metadata.')
    print(1)
  }

  col_check2 <- c("Type", "Feature") %in% colnames(expr)
  if (col_check2[1] & col_check2[2]) {
    col_check3 <- colnames(expr)[3:ncol(expr)] %in% meta_final$library
    if (length(which(!col_check3)) == 0) {
      expr_pro <- expr[grep("Gene|Protein|protein", expr[, 1]), 2:ncol(expr)]
      expr_pep <- expr[grep("Peptide|peptide", expr[, 1]), 2:ncol(expr)]
      if (nrow(expr_pep) == 0) {
        expr_pep <- NULL
      }
    } else {
      stop("The column names does not correspond to input metadata.")
    }
  } else {
    message("The first column of your expression data is used as features.")
    col_check3 <- colnames(expr)[2:ncol(expr)] %in% meta_final$library
    if (length(which(!col_check3)) == 0) {
      expr_pro <- expr
      expr_pep <- NULL
      colnames(expr_pro) <- c("Feature", colnames(expr_pro)[2:ncol(expr_pro)])
    } else {
      stop("The column names does not correspond to input metadata.")
    }
  }

  # Check if data at peptide levels provided -----------------
  if (!is.null(expr_pep)) {
    data_list <- list(
      "expdata_proteinLevel" = expr_pro,
      "expdata_peptideLevel" = expr_pep,
      "metadata" = meta_final
    )

  } else {
    message("You only input data at protein levels.")
    data_list <- list(
      "expdata_proteinLevel" = expr_pro,
      "metadata" = meta_final
    )

  }

  return(data_list)
}
