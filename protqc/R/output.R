#' Generating a table of conclusion
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @import stats
#' @import utils
#' @importFrom data.table data.table
#' @importFrom psych geometric.mean
#' @export

qc_conclusion <- function(exp_path, meta_path, output_dir = NULL, plot=FALSE) {
  # Load historical QC results -------------------------
  load(system.file("data/historical_qc.rda", package = "protqc"))
  load(system.file("data/historical_qc_norm.rda", package = "protqc"))
  load(system.file("data/historical_qc_stat.rda", package = "protqc"))
  ref_qc <- historical_qc
  ref_qc_norm <- historical_qc_norm
  ref_qc_stat <- historical_qc_stat

  # Load the input data --------------------------------
  data_list <- input_data(exp_path, meta_path)
  pro_data <- data_list$expdata_proteinLevel
  meta <- data_list$metadata
  if (length(data_list) == 3) {
    pep_data <- data_list$expdata_peptideLevel
  }else {
    pep_data <- NULL
  }

  # Run the QC pipelines --------------------------------
  allmetrics_results <- qc_allmetrics(pro_data, meta, pep_data, output_dir)
  allmetrics_dt <- allmetrics_results$output_table
  output_list <- qc_total(allmetrics_dt, ref_qc, ref_qc_norm, ref_qc_stat)
  output_table <- output_list$Raw

  # Cut-off --------------------------------------------
  total_norm <- output_table$Value[nrow(output_table)]
  total_ref_norm <- c(total_norm, as.numeric(ref_qc_norm$Total_norm))
  his_ref_norm <- as.numeric(ref_qc_norm$Total_norm)
  if (total_norm == 10) {
    total_cf <- "100%"
  } else if (total_norm %in% his_ref_norm) {
    total_num <- length(his_ref_norm)
    total_pos <- floor(rank(- his_ref_norm)[his_ref_norm == total_norm])
    total_perc <- (total_num - total_pos) / total_num
    total_cf <- paste((round(total_perc, 4)) * 100, "%", sep = "")
  } else {
    total_num <- length(total_ref_norm)
    total_pos <- floor(rank(- total_ref_norm)[1])
    total_perc <- (total_num - total_pos) / total_num
    total_cf <- paste((round(total_perc, 4)) * 100, "%", sep = "")
  }
  output_cutoff <- data.table(
    "Cut-off" = c("0%", "20%", "50%", "80%", "100%", total_cf),
    "Percentile" = c(quantile(his_ref_norm, na.rm = T), total_norm)
  )

  # Save & Output --------------------------------------
  output_table2 <- output_list$Normalized
  output_table_o <- output_table2[order(output_table2$`Quality Metrics`), ]
  ref_qc_norm_new <- rbind(ref_qc_norm, c("QUERIED DATA", output_table_o$Value))

  if (!is.null(output_dir)) {
    output_dir_cutoff <- file.path(output_dir, "cutoff_table.tsv")
    output_dir_rank <- file.path(output_dir, "rank_table.tsv")
    output_dir_final <- file.path(output_dir, "conclusion_table.tsv")
    write.table(ref_qc_norm_new, output_dir_rank, sep = "\t", row.names = F)
    write.table(output_table, output_dir_final, sep = "\t", row.names = F)
    write.table(output_cutoff, output_dir_cutoff, sep = "\t", row.names = F)
  }

  final_list <- list(
    results = allmetrics_results,
    conclusion = output_table
  )

  return(final_list)
}
