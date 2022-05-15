#' Linear normalization: 1~10
#' @param x A value to be normalized
#' @param x_max Queried maximum in historical values
#' @param x_min Queried minimum in historical values
#' @param decreasing if True, a decreased linear normalization will be applied.
#' @export
qc_linear_norm <- function(x, x_min, x_max, decreasing=F) {
  if (length(x) == 1) {
    if (x <= x_min) {
      x <- x_min
    } else if (x >= x_max) {
      x <- x_max
    }
  }
  if (decreasing == F) {
    x_norm <- sapply(x, function(a) 1 + (a - x_min) * 9 / (x_max - x_min))
  } else {
    x_norm <- sapply(x, function(a) 10 - (a - x_min) * 9 / (x_max - x_min))
  }

  return(x_norm <- round(x_norm, digits = 3))

}

#' Label the performance: bad, fair, good, great
#' @param x The value of the metric to be labeled
#' @param x_ref The reference values of the metric
#' @param cutoff The cut-offs of the labeling
#' @import stats
#' @importFrom dplyr between
#' @export

qc_performance <- function(x, x_ref, cutoff=c(0, 0.2, 0.5, 0.8, 1)) {
  ref_perc <- quantile(x_ref, cutoff, na.rm = T)
  if (between(x, ref_perc[1], ref_perc[2])) {
    x_class <- "Bad"
  } else if (between(x, ref_perc[2], ref_perc[3])) {
    x_class <- "Fair"
  } else if (between(x, ref_perc[3], ref_perc[4])) {
    x_class <- "Good"
  } else if (between(x, ref_perc[4], ref_perc[5])) {
    x_class <- "Great"
  }

  return(x_class)

}

#' Rank in historical performances
#' @param x The value of the metric to be labeled
#' @param x_hist The reference values of the metric
#' @export
qc_rank <- function(x, x_hist) {
  x_all <- c(x, x_hist[!is.na(x_hist)])
  x_pos <- floor(rank(- x_all)[1])
  x_rank <- c(paste(x_pos, "/", length(x_all), sep = ""))

  return(x_rank)

}

#' Calculating all QC metrics
#' @param pro_dt A expression table file (at protein level)
#' @param pep_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory for results
#' @param plot if True, a plot will be output.
#' @export

qc_allmetrics <- function(pro_dt, meta_dt, pep_dt=NULL,
                          output_dir=NULL, plot=FALSE) {
  # Basic information ---------------------------
  pro_info <- qc_info(pro_dt, meta_dt)

  # SNR -----------------------------------------
  snr_results <- qc_snr(pro_dt, meta_dt, output_dir, plot)
  snr_value <- snr_results$SNR

  # RC ------------------------------------------
  if (!is.null(pep_dt)) {
    cor_results <- qc_cor(pep_dt, meta_dt, output_dir, plot)
    cor_value <- cor_results$COR
  }else {
    cor_results <- NULL
    cor_value <- NA
  }

  # QC results ----------------------------------
  metrics <- c(
    "Number of features",
    "Missing percentage (%)",
    "Absolute Correlation",
    "Coefficient of variantion (CV, %)",
    "Signal-to-Noise Ratio (SNR)",
    "Relative Correlation with Reference Datasets (RC)"
  )
  qc_values <- c(pro_info, snr_value, cor_value)

  # Output --------------------------------------
  output_table <- data.table(
    "Quality Metrics" = metrics,
    "Value" = c(pro_info, snr_value, cor_value)
  )

  all_results <- list(
    snr_results = snr_results,
    cor_results = cor_results,
    output_table = output_table
  )
  return(all_results)
}


#' Calculating: Total Score
#' @param allmetrics_dt the output table from qc_allmetrics
#' @param ref_qc historical data set
#' @param ref_qc_norm historical data set
#' @param ref_qc_stat historical data set
#' @param normalized if True, the qc values will be linearly normalized to 1~10.
#' @export

qc_total <- function(allmetrics_dt,
                     ref_qc, ref_qc_norm, ref_qc_stat, normalized=T) {
  # Normalize & Rank: All metrics ----------------------
  output_class <- c()
  output_norm <- c()
  output_rank <- c()
  metrics <- allmetrics_dt$`Quality Metrics`
  for (m in metrics) {
    x <- allmetrics_dt$Value[allmetrics_dt$`Quality Metrics` %in% m]
    x_ref_norm <- as.numeric(ref_qc_norm[, colnames(ref_qc_norm) %in% m])
    x_ref <- ref_qc[, colnames(ref_qc) %in% m]

    if (!is.na(x)) {
      x_max <- max(x_ref, na.rm = T)
      x_min <- min(x_ref, na.rm = T)
      if (m %in% c("Coefficient of variantion (CV, %)",
                   "Missing percentage (%)")) {
        x_norm <- qc_linear_norm(x, x_min, x_max, decreasing = T)
      } else {
        x_norm <- qc_linear_norm(x, x_min, x_max)
      }
      x_rank <- qc_rank(x_norm, x_ref_norm)
      x_class <- qc_performance(x_norm, x_ref_norm)

    } else {
      x_norm <- NA
      x_rank <- NA
      x_class <- NA
    }

    output_norm <- c(output_norm, x_norm)
    output_rank <- c(output_rank, x_rank)
    output_class <- c(output_class, x_class)
  }

  # Normalize & Rank: Total score ----------------------
  total_ref <- as.numeric(ref_qc_norm$Total)
  total_value <- round(geometric.mean(as.numeric(output_norm)), 3)
  total_norm <- qc_linear_norm(total_value, min(total_ref), max(total_ref))

  total_ref_norm <- as.numeric(ref_qc_norm$Total_norm)
  total_rank <- qc_rank(total_norm, total_ref_norm)
  total_c <- qc_performance(total_norm, total_ref_norm)

  # Output ---------------------------------------------
  allmetrics <- c(metrics, "Total", "Total_norm")
  allnorm_dt <- data.table(
    "Quality Metrics" = allmetrics,
    "Value" = c(output_norm, total_value, total_norm)
  )
  allmetrics_dt <- rbind(
    allmetrics_dt,
    data.table(
      "Quality Metrics" = "Total Score",
      "Value" = total_norm
    )
  )
  output_table <- merge(
    allmetrics_dt,
    ref_qc_stat,
    by = "Quality Metrics"
  )
  output_table <- cbind(
    output_table[c(4, 3, 1, 2, 6, 5, 7), ],
    data.table(
      "Rank" = c(output_rank, total_rank),
      "Performance" = c(output_class, total_c)
    )
  )
  output_list <- list(
    Normalized = allnorm_dt,
    Raw = output_table
  )
  return(output_list)
}
