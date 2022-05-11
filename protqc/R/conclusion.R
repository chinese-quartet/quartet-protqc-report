#' Generating a table of conclusion
#'
#' @param pro_path A file path of the expression table file (at protein level)
#' @param pep_path A file path of the expression table file (at peptide level)
#' @param meta_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @import data.table
#' @import ggplot2
#' @importFrom psych geometric.mean
#' @export

table_conclusion <- function(pro_path, meta_path,
                             output_dir = NULL, pep_path = NULL){
  ref_qc_dir <- file.path(
    system.file(package = "protqc"), "data/historical_qc.rds")
  ref_qc_norm_dir <- file.path(
    system.file(package = "protqc"), "data/historical_qc_norm.rds")
  ref_qc <- readRDS(ref_qc_dir)
  ref_qc_norm <- readRDS(ref_qc_norm_dir)

  data_list <- input_data(pro_path,meta_path,pep_path)
  pro_data <- data_list$expdata_proteinLevel
  meta <- data_list$metadata

  pro_info <- qc_info(pro_data,meta)
  snr_results <- qc_snr(pro_data, meta, output_dir)
  snr_value <- snr_results$SNR

  if(!is.null(pep_path)) {
    pep_data <- data_list$expdata_peptideLevel
    cor_results <- qc_cor(pep_data, meta, output_dir)
    cor_value <- cor_results$COR
  }else {
    cor_value <- NA
  }

  output_table <- rbind(pro_info, data.table(
    "Quality Metrics" = c(
      'Signal-to-Noise Ratio (SNR)',
      'Relative Correlation with Reference Datasets (RC)'),
    "Value" =c(snr_value, cor_value)
  ))

  output_class <- c()
  output_norm <- c()
  output_rank <- c()
  output_ms <- c()
  metrics <- output_table$`Quality Metrics`
  for(m in metrics) {
    x <- output_table$Value[output_table$`Quality Metrics` %in% m]
    x_ref_norm <- as.numeric(ref_qc_norm[, colnames(ref_qc_norm) %in% m])
    x_ref <- ref_qc[, colnames(ref_qc) %in% m]

    x_ref_mean <- round(mean(x_ref, na.rm = T),3)
    x_ref_sd <- round(sd(x_ref, na.rm = T),3)
    x_ref_ms <- paste(x_ref_mean,' ± ',x_ref_sd,sep = '')

    x_max <- max(x_ref, na.rm = T)
    x_min <- min(x_ref, na.rm = T)
    if(m %in% c('Coefficient of variantion (CV, %)','Missing percentage (%)')){
      x_norm <- qc_linear_norm(x, x_min, x_max, decreasing = T)
    }else {
      x_norm <- qc_linear_norm(x, x_min, x_max)
    }
    x_all <- c(x_norm,x_ref_norm[!is.na(x_ref_norm)])
    x_pos <- floor(rank(-x_all)[1])
    x_rank <- c(paste(x_pos,'/', length(x_all),sep = ''))

    x_ref_perc <- quantile(x_ref_norm, c(0, 0.2, 0.5, 0.8, 1), na.rm = T)
    if(between(x_norm,x_ref_perc[1],x_ref_perc[2])) {x_class <- 'Bad'
    }else if(between(x_norm,x_ref_perc[2],x_ref_perc[3])) {x_class <- "Fair"
    }else if(between(x_norm,x_ref_perc[3],x_ref_perc[4])) {x_class <- "Good"
    }else if(between(x_norm,x_ref_perc[4],x_ref_perc[5])) x_class <- "Great"

    output_ms <- c(output_ms, x_ref_ms)
    output_norm <- c(output_norm, x_norm)
    output_rank <- c(output_rank, x_rank)
    output_class <- c(output_class, x_class)
  }

  his_ref <- as.numeric(ref_qc_norm$Total)
  total_value <- round(geometric.mean(as.numeric(output_norm)), 3)
  total_norm <- qc_linear_norm(total_value, min(his_ref), max(his_ref))

  output_table <- rbind(output_table, data.table(
    "Quality Metrics" = "Total Score",
    "Value" = total_norm
  ))

  ref_qc_norm_new <- rbind(
    ref_qc_norm, c(
      'QUERIED DATA',
      output_norm[c(3, 4, 2, 1, 6, 5)],
      total_value,total_norm
    )
  )

  his_mean <- round(mean(his_ref, na.rm = T),3)
  his_sd <- round(sd(his_ref, na.rm = T),3)
  his_ms <- paste(his_mean,' ± ',his_sd,sep = '')

  total_ref <- as.numeric(ref_qc_norm_new$Total_norm)
  total_pos <- floor(rank(-total_ref)[nrow(ref_qc_norm_new)])
  total_rank <- c(paste(total_pos,'/', length(total_ref),sep = ''))
  total_ref_perc <- quantile(total_ref, c(0, 0.2, 0.5, 0.8, 1))
  if(between(total_norm,total_ref_perc[1],total_ref_perc[2])) {
    total_class <- 'Bad'
  }else if(between(total_value,total_ref_perc[2],total_ref_perc[3])) {
    total_class <- "Fair"
  }else if(between(total_value,total_ref_perc[3],total_ref_perc[4])) {
    total_class <- "Good"
  }else if(between(total_value,total_ref_perc[4],total_ref_perc[5])) {
    total_class <- "Great"
  }

  output_table_final <- cbind(
    output_table,
    data.table(
      "Historical Value (mean ± SD)" = c(output_ms, his_ms),
      "Rank" = c(output_rank, total_rank),
      "Performance" = c(output_class,total_class)
    )
  )

  total_cf <- paste((1-round(total_pos/length(total_ref),4))*100,'%',sep = '')
  output_cutoff <- data.table(
    'Cut-off' = c('0%', '20%', '50%', '80%', '100%', total_cf),
    'Percentile' = c(total_ref_perc, total_norm)
  )

  if(!is.null(output_dir)){
    output_dir_cutoff <- file.path(output_dir, 'cutoff_table.tsv')
    output_dir_rank <- file.path(output_dir,'rank_table.tsv')
    output_dir_final <- file.path(output_dir,'conclusion_table.tsv')
    write.table(ref_qc_norm_new,output_dir_rank,sep = '\t',row.names = F)
    write.table(output_table_final,output_dir_final,sep = '\t',row.names = F)
    write.table(output_cutoff, output_dir_cutoff, sep = '\t',row.names = F)
  }

  return(list(results = ref_qc_norm_new, conclusion = output_table_final))
}
