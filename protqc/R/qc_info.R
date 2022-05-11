#' Statistics for basic information
#'
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @importFrom psych corr.test
#' @import data.table
#' @export

qc_info <- function(expr_dt, meta_dt) {

  d <- expr_dt
  d[d == 0] <- NA
  stat_num <- length(unique(d[, 1]))

  d_missings <- apply(d, 2, function(x) length(which(is.na(x))))
  stat_missing <- round(median(d_missings)*100 / stat_num, 3)

  d_cortest <- corr.test(d[, 2:ncol(d)], method = 'pearson', adjust="fdr")
  d_cormtx <- d_cortest$r
  d_cormtx[d_cortest$p > 0.05] <- 0
  samples <- table(meta_dt$sample)
  samples <- samples[samples > 1]
  m <- meta_dt[, colnames(meta_dt) %in% c('library', 'sample')]
  d_cordf <- melt(d_cormtx)
  d_cordf <- d_cordf[d_cordf$Var2 != d_cordf$Var1, ]
  d_cordf <- merge(d_cordf, m, by.x = 'Var1', by.y = 'library')
  d_cordf <- merge(d_cordf, m, by.x = 'Var2', by.y = 'library')
  d_cordf <- d_cordf[d_cordf$sample.x == d_cordf$sample.y, ]
  stat_repcor <- round(median(d_cordf$value), 3)

  d_long <- melt(d)
  d_long$value2 <- 2^(d_long$value)
  d_long <- na.omit(d_long)
  colnames(d_long) <- c('feature', "variable", "value", "value2")
  d_cv <- aggregate(
    value2 ~ feature, data = d_long,
    FUN = function(x) sd(x)/mean(x))
  stat_cv <- round(median(d_cv$value2, na.rm = T) * 100, 3)

  stat_all <- data.table(
    "Quality Metrics" = c(
      'Number of features',
      'Missing percentage (%)',
      'Absolute Correlation',
      'Coefficient of variantion (CV, %)'
    ),
    "Value" = c(
      stat_num,
      stat_missing,
      stat_repcor,
      stat_cv)
  )

  return(stat_all)
}
