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

  ref_snrcorr_dir <- file.path(system.file(package = "ProtQC"), "data/ref_snrcorr.rds")
  snrcorr <- readRDS(ref_snrcorr_dir)
  snrcorr_new <- rbind(snrcorr,data.frame(batch='Lot2_test',SNR=output_conclusion$Value[1],COR=output_conclusion$Value[2]))

  snr_max <- max(snrcorr_new$SNR)
  snr_min <- min(snrcorr_new$SNR)
  snrcorr_new$SNR_normalized <- sapply(snrcorr_new$SNR,function(x){1+(x-snr_min)*99/(snr_max-snr_min)})

  cor_max <- max(snrcorr_new$COR)
  cor_min <- min(snrcorr_new$COR)
  snrcorr_new$COR_normalized <- sapply(snrcorr_new$COR,function(x){1+(x-cor_min)*99/(cor_max-cor_min)})

  snrcorr_new$Total <- apply(snrcorr_new[,4:5],1,mean)
  total_quantiles <- quantile(snrcorr_new$Total)
  snrcorr_new$Class <- sapply(snrcorr_new$Total,function(x){
    if(between(x,total_quantiles[1],total_quantiles[2])){return("low")
      }else if(between(x,total_quantiles[2],total_quantiles[3])){return("mid-low")
        }else if(between(x,total_quantiles[3],total_quantiles[4])){return("mid-high")
          }else if(between(x,total_quantiles[4],total_quantiles[5])) return("high")
  })

  output_dir_rank <- file.path(output_dir,'rank_table.tsv')
  output_dir_final <- file.path(output_dir,'conclusion_table.tsv')
  write.table(snrcorr_new,output_dir_rank,sep = '\t',row.names = F)
  write.table(output_conclusion,output_dir_final,sep = '\t',row.names = F)
}
