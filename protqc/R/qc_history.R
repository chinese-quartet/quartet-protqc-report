#' Quality assessment of historical data sets
#'
#' @param historical_pro_path A file path of data (at protein level)
#' @param historical_pep_path A file path of data (at peptide level)
#' @param meta_dt_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @import data.table
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom psych corr.test
#' @importFrom psych geometric.mean
#' @export

qc_history <- function(historical_pro_path, historical_pep_path,
                       historical_meta_path){
  all_pro <- readRDS(historical_pro_path)
  all_pep <- readRDS(historical_pep_path)
  all_meta <- readRDS(historical_meta_path)

  df_history <- c()
  batches <- unique(all_meta$batch)
  for(b in batches) {
    meta <- all_meta[all_meta$batch %in% b, ]
    pro_data <- all_pro[, colnames(all_pro) %in% c(meta$library, 'Gene')]
    pro_data <- pro_data[apply(
      pro_data, 1, function(x) length(which(is.na(x))) < ncol(pro_data) - 1), ]
    pro_info <- qc_info(pro_data,meta)
    snr_results <- qc_snr(pro_data,meta)
    snr_value <- snr_results$SNR

    if(grepl('Lot2', b)) {
      pep_data <- all_pep[, colnames(all_pep) %in% c(meta$library, 'Sequence')]
      pep_data <- pep_data[apply(
        pep_data, 1, function(x) length(which(is.na(x))) < ncol(pep_data) - 1),]
      cor_results <- qc_cor(pep_data,meta)
      cor_value <- cor_results$COR
    }else {
      cor_value <- NA
    }

    df_per_batch <- rbind(pro_info, data.table(
      "Quality Metrics" = c(
        'Signal-to-Noise Ratio (SNR)',
        'Relative Correlation with Reference Datasets (RC)'),
      "Value" =c(snr_value, cor_value)
    ))
    df_per_batch$Batch <-b
    df_history <- rbind(df_history, df_per_batch)
  }

  df_long <- melt(df_history, id = colnames(df_history)[c(1,3)])
  df_wide <- dcast(df_long, Batch ~ `Quality Metrics`)

  df_wide_norm <- c()
  metrics <- colnames(df_wide)[2:ncol(df_wide)]
  for(m in metrics) {
    x <- df_wide[, colnames(df_wide) %in% m]
    x_max <- max(x, na.rm = T)
    x_min <- min(x, na.rm = T)
    if(m %in% c('Coefficient of variantion (CV, %)','Missing percentage (%)')){
      x_norm <- qc_linear_norm(x, x_min, x_max, decreasing = T)
    }else {
      x_norm <- qc_linear_norm(x, x_min, x_max)
    }
    df_wide_norm <- cbind(df_wide_norm, x_norm)
  }

  df_wide_norm <- cbind(df_wide$Batch, df_wide_norm,
                        apply(df_wide_norm, 1, function(x) {
                          round(geometric.mean(as.numeric(x)), 3)})
  )
  colnames(df_wide_norm) <- c(colnames(df_wide), 'Total')
  df_wide_norm <- as.data.frame(df_wide_norm)

  total_all <- as.numeric(df_wide_norm$Total)
  total_max <- max(total_all, na.rm = T)
  total_min <- min(total_all, na.rm = T)
  df_wide_norm$Total_norm <- qc_linear_norm(total_all, total_min, total_max)

  return(list(
    'historical_qc_statisctics' = df_wide,
    'historical_qc_norm' = df_wide_norm)
  )

}
