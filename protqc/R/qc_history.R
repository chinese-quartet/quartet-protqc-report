#' Quality assessment of historical data sets
#' @import stats
#' @import utils
#' @importFrom data.table data.table
#' @importFrom stringi stri_escape_unicode
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom psych corr.test
#' @importFrom psych geometric.mean
#' @export

qc_history <- function() {
  # Load data ------------------------------------------
  load(system.file("data/historical_meta.rda", package = "protqc"))
  load(system.file("data/historical_meta.rda", package = "protqc"))
  load(system.file("data/historical_meta.rda", package = "protqc"))

  # Global variables -----------------------------------
  all_pro <- historical_data_genesymbols
  all_pep <- historical_data_peptides
  all_meta <- historical_meta

  # Run the QC pipeline: per batch ---------------------
  df_history <- c()
  batches <- unique(all_meta$batch)
  for (b in batches) {
    meta <- all_meta[all_meta$batch %in% b, ]
    pro_dt <- all_pro[, colnames(all_pro) %in% c(meta$library, "Gene")]
    sample_num <- ncol(pro_dt) - 1
    pro_dt <- pro_dt[apply(X = pro_dt,
                           MARGIN = 1,
                           function(x) length(which(is.na(x))) < sample_num), ]
    if (grepl("Lot2", b)) {
      pep_dt <- all_pep[, colnames(all_pep) %in% c(meta$library, "Sequence")]
      sample_num <- ncol(pep_dt) - 1
      pep_dt <- pep_dt[apply(X = pep_dt,
                             MARGIN = 1,
                             function(x) length(which(is.na(x))) < sample_num), ]
    } else {
      pep_dt <- NULL
    }

    allmetrics_results <- qc_allmetrics(pro_dt, meta, pep_dt)
    output_table <- allmetrics_results$output_table
    df_per_batch <- data.table("Batch" = b, output_table)
    df_history <- rbind(df_history, df_per_batch)
  }

  # Converting: long to wide ---------------------------
  df_long <- melt(df_history, id = colnames(df_history)[c(1, 2)])
  df_wide <- dcast(df_long, Batch ~ `Quality Metrics`)

  # Normalizing: 1 ~ 10 --------------------------------
  df_wide_norm <- c()
  df_meansd <- c()
  metrics <- colnames(df_wide)[2:ncol(df_wide)]
  for(m in metrics) {
    x <- df_wide[, colnames(df_wide) %in% m]
    x_mean <- round(mean(x, na.rm = T), 3)
    x_sd <- round(sd(x, na.rm = T), 3)
    x_ms <- paste(x_mean, " \\u00b1 ", x_sd, sep = "")

    x_max <- max(x, na.rm = T)
    x_min <- min(x, na.rm = T)
    if (m %in% c("Coefficient of variantion (CV, %)",
                 "Missing percentage (%)")) {
      x_norm <- qc_linear_norm(x, x_min, x_max, decreasing = T)
    }else {
      x_norm <- qc_linear_norm(x, x_min, x_max)
    }
    df_wide_norm <- cbind(df_wide_norm, x_norm)

    m_ms <- data.table("Quality Metrics" = m,
                       "Historical Value (mean \\u00b1 SD)" = x_ms)
    df_meansd <- rbind(df_meansd, m_ms)
  }

  total_norm <- apply(X = df_wide_norm,
                      MARGIN = 1,
                      function(x) round(geometric.mean(as.numeric(x)), 3))
  df_wide_norm <- cbind(df_wide$Batch, df_wide_norm, total_norm)

  colnames(df_wide_norm) <- c(colnames(df_wide), "Total")
  df_wide_norm <- as.data.frame(df_wide_norm)

  total_all <- as.numeric(df_wide_norm$Total)
  total_max <- max(total_all, na.rm = T)
  total_min <- min(total_all, na.rm = T)
  df_wide_norm$Total_norm <- qc_linear_norm(total_all, total_min, total_max)

  his_ref_norm <- as.numeric(df_wide_norm$Total_norm)
  his_mean <- round(mean(his_ref_norm, na.rm = T), 3)
  his_sd <- round(sd(his_ref_norm, na.rm = T), 3)
  his_ms <- paste(his_mean, " \\u00b1 ", his_sd, sep = "")
  df_meansd_final <- rbind(
    df_meansd,
    data.table("Quality Metrics" = "Total Score",
               "Historical Value (mean \\u00b1 SD)" = his_ms))

  # Output ------------------------------------------
  historical_qc_stat <- df_meansd_final
  historical_qc <- df_wide
  historical_qc_norm <- df_wide_norm

  output_list <- list(statistics = historical_qc_stat,
                      raw = historical_qc,
                      normalized = historical_qc_norm)

  return(output_list)
}
