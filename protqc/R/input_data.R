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
    # Compatible with version 0.1.x.
    # Caution: The data format is changed from version 0.2.0.
    # metadata file: name -> library, sample
    # data file: rowname -> Feature, xxx
    expr_colnames <- colnames(expr)
    colnames(expr) <- c("Feature", expr_colnames[2:length(expr_colnames)])
    if ("name" %in% colnames(meta)) {
      library <- meta$name
    } else {
      library <- meta$library
    }
    metadata <- data.frame(library=library, sample=meta$sample)
    data_list <- list(
      'expdata_proteinLevel' = expr,
      'metadata' = metadata
    )
  }

  return(data_list)
}
