#' Generating a table of conclusion
#'
#' @param expr_dt_path A file path of the expression table file
#' @param meta_dt_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr %>%
#' @importFrom ggthemes theme_few
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export

table_conclusion <- function(expr_dt_path,meta_dt_path,output_dir){
  output_signoise_db <- plot_pca(expr_dt_path,meta_dt_path,output_dir)
  output_cor_value <- plot_corr(expr_dt_path,meta_dt_path,output_dir)
  output_conclusion <- rbind(output_signoise_db,output_cor_value)

  output_dir_final <- file.path(output_dir,'conclusion_table.tsv')
  write.table(output_conclusion,output_dir_final,sep = '\t',row.names = F)
}
