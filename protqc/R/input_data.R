#' Input files
#'
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @export

input_data <- function(exp_path, meta_path) {
  expr <- read.csv(exp_path)
  meta <- read.csv(meta_path)
  expr_pro <- expr[grep("Gene|Protein|protein", expr[, 1]), 2:ncol(expr)]
  expr_pep <- expr[grep("Peptide|peptide", expr[, 1]), 2:ncol(expr)]

  if(nrow(expr_pep) != 0 & nrow(expr_pro) != 0) {
    data_list <- list(
      'expdata_proteinLevel' = expr_pro,
      'expdata_peptideLevel' = expr_pep,
      'metadata' = meta
    )
  }else if(nrow(expr_pep) == 0 & nrow(expr_pro) != 0){
    data_list <- list(
      'expdata_proteinLevel' = expr_pro,
      'metadata' = meta
    )
  }else {
    print("Please check your input format.")
  }

  return(data_list)
}
