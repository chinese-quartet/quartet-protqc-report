#' Input files
#'
#' @param pro_path A file path of the expression table file (at protein level)
#' @param pep_path A file path of the expression table file (at peptide level)
#' @param meta_path A file path of the metadata file
#' @export

input_data <- function(pro_path, meta_path, pep_path = NULL) {
  expr_pro <- read.csv(pro_path)
  meta <- read.csv(meta_path)

  if(is.null(pep_path)) {
    data_list <- list(
      'expdata_proteinLevel' = expr_pro,
      'metadata' = meta
    )
  }else {
    expr_pep <- read.csv(pep_path)
    data_list <- list(
      'expdata_proteinLevel' = expr_pro,
      'expdata_peptideLevel' = expr_pep,
      'metadata' = meta
    )
  }

  return(data_list)
}
